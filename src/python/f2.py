#!/usr/bin/env python3
# Copyright (c) 2018, Panos Theodoropoulos

import argparse
import time
import math
import random
import os
import os.path
import shlex
import subprocess
import networkx as nx
from networkx.algorithms import bipartite
import sympy as sp
import scipy.optimize
import logging as log
import numpy as np
import pickle
from collections import defaultdict
import resource
import sys
# import lapjv
import lap
import cProfile
import re

import operator as op
from functools import reduce
from itertools import combinations
#should be faster than scipy, but crashing on import...
# from pymatgen.optimization import linear_assignment
from scipy.optimize import linear_sum_assignment
VERBOSE = False

def do_cprofile(func):
    '''
    https://zapier.com/engineering/profiling-python-boss/
    '''
    def profiled_func(*args, **kwargs):
        profile = cProfile.Profile()
        try:
            profile.enable()
            result = func(*args, **kwargs)
            profile.disable()
            return result
        finally:
            profile.print_stats()
    return profiled_func


class Config:
    def __init__(self):
        self.formula_filename = None
        self.has_ind_header = None
        self.var_set = None  # list of variables used in the constraints
        self.working_dir = None
        self.keep_cnf = False
        self.skip_sharpsat = False
        self.random_seed = None
        self.mode = None
        self.default_error_prob = 0.05
        self.total_max_time = math.nan
        # self.f2_program_git_version = ''  # DEVTODO: Autocomplete
        # self.cmsat_git_version = ''  # DEVTODO: Autocomplete
        self.cmsat_major_version = 5

        # Preferred time allocated for each solver invocation, actual may be
        # larger (if iterations are few) or smaller (if total time is being
        # exhausted)
        self.alg1_loop_timeout = None
        self.degn = math.nan  # Average variable degree of LDPC code
        self.alg1_maxsol = 4 #this is s, the solution cutoff, in algorithm 1 SPARSE-COUNT in the paper "Adaptive Hashing for Model Counting"
        # self.alg1_maxsol/2 = 128
        self.alg1_threshold_coeff = 2
        assert(self.alg1_threshold_coeff == self.alg1_maxsol/2)
        self.alg1_n_iter = math.nan
        self.alg2_n_iter = math.nan
        self.alg2_loop_timeout = math.nan
        self.alg2_maxsol = 999999999
        self.alg2_error_prob = math.nan
        self.f2_alg_confidence_radius = math.nan
        self.f2_alg_error_prob = math.nan
        self.args = None
        self.cmsat_exe = None
        self.sharpsat_exe = None
        self.extra_configs = None

    def setup(self, args):
        self.extra_configs = args.extra_configs
        global formula_lines
        self.mode = args.mode
        if not os.path.isfile(args.input_filename):
            log.error('File {} does not exist.'.format(args.input_filename))
            print('File {} does not exist.'.format(args.input_filename))
            exit(101)
        else:
            self.formula_filename = args.input_filename
            with open(args.input_filename, 'r', errors='replace') as f:
                formula_lines = f.readlines()
                # print(formula_lines)
                self.var_set, self.has_ind_header = \
                    extract_support_vars()
        print("in setup, args.working_dir:", args.working_dir)                        
        if args.working_dir == '':
            conf.working_dir = os.getcwd()
            print("conf.working_dir:", conf.working_dir)
        elif not os.path.isdir(args.working_dir):
            log.error('Working directory does not exist.')
            exit(102)
        else:
            conf.working_dir = args.working_dir
            os.chdir(conf.working_dir)
        print("in setup, conf.working_dir:", conf.working_dir)                                    
        self.keep_cnf = args.keep_cnf
        # self.keep_cnf = False
        self.skip_sharpsat = args.skip_sharpsat
        if args.random_seed >= 0:
            self.random_seed = args.random_seed
            random.seed(args.random_seed)
            np.random.seed(args.random_seed)
        else:
            random.seed()
            rs = random.randint(0, 2**31)
            random.seed(rs)
            np.random.seed(rs)            
            self.random_seed = rs
        if (not math.isnan(args.lower_bound)) and\
                (not math.isnan(args.upper_bound)):
            if args.lower_bound >= math.isnan(args.upper_bound):
                log.error('lower_bound must be less than upper_bound')
                exit(103)
        self.degn = args.var_degree
        self.cmsat_major_version = args.cmsat_version
        self.total_max_time = args.max_time
        if args.cmsat_exe is not None:
            self.cmsat_exe = os.path.join('.', args.cmsat_exe)
            if not os.path.isfile(self.cmsat_exe):
                log.error('CMSat executable does not exist')
                exit(104)
        else:
            f2_dir = os.path.dirname(os.path.realpath(__file__))
            cms = os.path.join(f2_dir, 'cryptominisat5')
            if not os.path.isfile(cms):
                log.error('CMSat executable does not exist in F2 dir')
                exit(105)
            else:
                self.cmsat_exe = cms
        if args.sharpsat_exe is not None:
            self.sharpsat_exe = os.path.join('.', args.sharpsat_exe)
            if not os.path.isfile(self.sharpsat_exe):
                log.error('sharpsat executable does not exist')
                exit(108)
        else:
            pass
            # f2_dir = os.path.dirname(os.path.realpath(__file__))
            # shs = os.path.join(f2_dir, 'sharpSAT')
            # if not os.path.isfile(shs):
            #     log.error('sharpSAT executable does not exist in F2 dir')
            #     exit(109)
            # else:
            #     self.sharpsat_exe = shs
        if math.isnan(args.error_prob):
            error_probability = conf.default_error_prob
        else:
            error_probability = args.error_prob
        if args.mode == 'lb':
            lb_mode_sanity_check(args)
            self.setup_alg1(args, error_probability)
        elif args.mode == 'ub':
            ub_mode_sanity_check(args)
            if not math.isnan(args.lower_bound):
                self.setup_alg2(args, error_probability)
            else:  # Lower bound not given
                if not math.isnan(args.lb_n_iter):
                    self.setup_alg1(args, math.nan)
                    self.setup_alg2(args, error_probability)
                else:
                    self.setup_alg1(args, error_probability / 2)
                    self.setup_alg2(args, error_probability / 2)
        elif args.mode == 'appr':
            num_of_phases = 1  # In how many parts to split the total error pr.
            if not math.isnan(args.lower_bound):
                num_of_phases += 1
            if not math.isnan(args.upper_bound):
                num_of_phases += 1
            er_prob = error_probability / num_of_phases
            if math.isnan(args.lower_bound):
                self.setup_alg1(args, er_prob)
            if math.isnan(args.upper_bound):
                self.setup_alg2(args, er_prob)
            self.setup_f2_alg(args, er_prob)
        else:
            assert False

    def setup_alg2(self, args, er_prob):
        """
        Setup the configuration for algorithm 2.
        If number of iterations is given, it is used. Otherwise it will be
        calculated from the available 'share' of error probability given as
        second parameter inside the algorithm
        :param args:
        :param er_prob: Error probability for the invocation of algorithm
        """
        if not math.isnan(args.ub_n_iter):
            self.alg2_n_iter = args.ub_n_iter
        else:
            self.alg2_error_prob = er_prob
        self.alg2_loop_timeout = args.max_time / min(5, self.alg2_n_iter)

    def setup_alg1(self, args, er_prob):
        """
        Setup the configuration for algorithm 1.
        If number of iterations is given, it is used. Otherwise it is
        calculated from the available 'share' of error probability given as
        second parameter.
        :param args:
        :param er_prob: Error probability for the invocation of algorithm
        """
        if not math.isnan(args.lb_n_iter):
            self.alg1_n_iter = args.lb_n_iter
        else:
            # self.alg1_n_iter = int(math.ceil(-8 * math.log(er_prob)))
            n_var = len(conf.var_set)
            self.alg1_n_iter = int(math.ceil(8 * math.log(math.ceil(math.log(n_var,2))/er_prob)))
        self.alg1_loop_timeout = args.max_time / min(5, self.alg1_n_iter)

    def setup_f2_alg(self, args, er_prob):
        self.f2_alg_confidence_radius = args.confidence_radius
        self.f2_alg_error_prob = er_prob


# === Globals ===
conf = Config()  # Global configuration of all algorithms
formula_lines = []  # The lines of the original cnf file
start_time = time.process_time()
aug_counter = 0  # Counter for augmented formulas filenames
subprocess_cumulative_time = 0  # Accumulate the time spent in subprocesses
# ===============


def lb_mode_sanity_check(args):
    if (not math.isnan(args.lb_n_iter) and
            not math.isnan(args.error_prob)):
        log.error('In \'lb\' mode --error_prob and --lb_n_iter '
                  'cannot be both specified.')
        exit(107)
    if math.isnan(args.lb_n_iter) and math.isnan(args.error_prob):
        args.error_prob = 0.05


def ub_mode_sanity_check(args):
    if (not math.isnan(args.ub_n_iter) and
            not math.isnan(args.error_prob)):
        log.error('In \'ub\' mode --error_prob and --ub_n_iter '
                  'cannot be both specified.')
        exit(106)
    if math.isnan(args.ub_n_iter) and math.isnan(args.error_prob):
        args.error_prob = 0.05


def check_time():
    # print('@'*80)
    # print("start_time", start_time)
    this_process_time = time.process_time() - start_time
    rt = conf.total_max_time - this_process_time - subprocess_cumulative_time
    if rt <= -1:
        print("check_time, rt =", rt)
        raise TotalTimeoutException()

def check_time_jdk():
    # print('@'*80)
    # print("start_time", start_time)
    global SAT_SOLVER_TIME_BETTER_PRECISION
    rt = conf.total_max_time - SAT_SOLVER_TIME_BETTER_PRECISION
    print('check_time_jdk, remaining time:', rt)
    if rt <= 0:
        raise TotalTimeoutException()


def remaining_time():
    """
    Calculate the remaining time allowed for the program to run
    """
    return remaining_time_jdk()
    this_process_time = time.process_time() - start_time
    rt = conf.total_max_time - this_process_time - subprocess_cumulative_time
    # return max(rt, 1)
    return rt

def remaining_time_jdk():
    """
    Calculate the remaining time allowed for the program to run
    """
    global SAT_SOLVER_TIME_BETTER_PRECISION
    rt = conf.total_max_time - SAT_SOLVER_TIME_BETTER_PRECISION
    assert(rt > 0)
    return max(rt, 0)

def consumed_time():
    """
    Total time spent by this process the the subprocesses it spawned
    :return: time in seconds
    """
    this_process_time = time.process_time() - start_time
    return this_process_time + subprocess_cumulative_time


class SolverBoundsExceededException(Exception):
    """
    A solver timeout occurred or the number of solutions exceeded a bound
    in a case that it is fatal for the algorithm.
    """
    pass


class TotalTimeoutException(Exception):
    """
    The total time allowed running time has expired.
    """
    pass


def main():
    global conf
    log.basicConfig(level=log.INFO)
    parser = argparse.ArgumentParser(
        prog='F2',
        usage='%(prog)s [options] formula_file.cnf',
        description=''
        'Probabilistically approximate number of models of a cnf formula')
    parser.add_argument('input_filename', type=str,
                        help='The formula in extended DIMACS form')
    parser.add_argument('--working-dir', type=str,
                        help='Where to create auxiliary files', default='')
    parser.add_argument('--keep-cnf', action='store_true',
                        help='Keep generated auxiliary formulas?')
    parser.add_argument('--random-seed', type=int,
                        help='Initialize the random generator', default=42)
    parser.add_argument('--var-degree', type=int,
                        help='Average variable degree of LDPC XORs',
                        default=12)
    parser.add_argument('--cmsat-exe', type=str,
                        help='The location of a cryptominisat executable',
                        default=None)
    parser.add_argument('--cmsat-version', type=int,
                        help='The major version of cryptominisat executable. '
                             '(Allowed versions 2 or 5, default 5)',
                        default=5)
    parser.add_argument('--sharpsat-exe', type=str,
                        help='The location of a SharpSAT executable',
                        default=None)
    parser.add_argument('--mode', type=str,
                        help='Mode of operation. Allowed values: lb, ub, appr',
                        choices=['lb', 'ub', 'appr'],
                        default='appr')
    parser.add_argument('--lower-bound', type=int,
                        help='Binary logarithm of lower bound (modes: ub, '
                             'appr)', default=math.nan)
    parser.add_argument('--upper-bound', type=int,
                        help='Binary logarithm of upper bound (only appr '
                             'mode)', default=math.nan)
    parser.add_argument('--error-prob', type=float,
                        help='Probability of error', default=math.nan)
    parser.add_argument('--confidence-radius', type=float,
                        help='Tolerance of error as a ratio of real value',
                        default=0.2)
    parser.add_argument('--max-time', type=int,
                        help='Maximum system running time (in seconds)',
                        default=3600)
    parser.add_argument('--lb-n-iter', type=int,
                        help='Override num of iterations in algorithm 1',
                        default=math.nan)
    parser.add_argument('--ub-n-iter', type=int,
                        help='Override num of iterations in algorithm 2',
                        default=math.nan)
    parser.add_argument('--skip-sharpsat', action='store_true',
                        help='Skip the SharpSAT invocation.')
    args = parser.parse_args()
    print("type(args):", type(args))
    print("(args):", args)
    conf.setup(args)
    run(args)


def run(args):
    try:
        if not conf.skip_sharpsat:
            log.info('Running SharpSAT for quick test...')
            tout, n, ct = sharpsat_count(conf.formula_filename, 2)
            if not tout:
                if n > 1:
                    print('F2: Exact count (by SharpSAT) is {} '
                          '(t={:.2f})'.format(n, ct))
                    print('F2Sharp:{}:{:.2f}:{:.2f}'.format(
                        conf.formula_filename, math.log2(n), ct))
                    if conf.has_ind_header:
                        print_sharpsat_warning()
                    exit(0)
                elif n == 0:
                    log.info(
                        'F2: The formula is UNSAT (as reported by sharpSAT)')
                    print('F2: UNSAT')
                    print('F2Sharp:{}:{:.2f}:{:.2f}'.format(
                        conf.formula_filename, 0, ct))
                    exit(0)
                else:  # n == 1
                    log.info('SharpSAT returned 1 which means that the formula'
                             ' is either UNSAT or has 1 solution.')
                    if is_unsat():
                        print('F2: UNSAT')
                        exit(0)
            else:
                log.info('SharpSAT timed out. Proceeding...')
        if conf.mode == 'lb':
            lb_init = prepare_lower_bound()
            lb = find_lower_bound(lb_init)
            pt = consumed_time()
            print('F2: Lower bound is 2^{} (processor time: {:.2f})'
                  ''.format(lb, pt))
        elif conf.mode == 'ub':
            if math.isnan(args.lower_bound):
                lb_init = prepare_lower_bound()
                lb = find_lower_bound(lb_init)
                ct1 = consumed_time() - 2  # Remove the sharpsat time
            else:
                lb = args.lower_bound
                ct1 = 0
            ub_a, ub_b = find_upper_bound(lb)
            ct2 = consumed_time()
            print('F2: Lower bound is 2^{} (t={:.2f}) and upper bound '
                  'is {:.2f} * 2^{:.2f} (total processor time: {:.2f})'
                  ''.format(lb, ct1, ub_a, ub_b, ct2))
            print('F2UB:{}:{}:{:.2f}:{:.2f}:{:.2f}'.format(
                  conf.formula_filename, lb, ct1,
                  math.log2(ub_a) + ub_b, ct2 - ct1))
        elif conf.mode == 'appr':
            if math.isnan(args.lower_bound):
                lb_init = prepare_lower_bound()
                lb = find_lower_bound(lb_init)
            else:
                lb = args.lower_bound
            if math.isnan(args.upper_bound):
                ub_a, ub_b = find_upper_bound(lb)
                ub = math.log2(ub_a) + ub_b
            else:
                ub = args.upper_bound
            if lb > ub:
                log.error('F2: Lower bound estimate cannot be greater than'
                          ' upper')
                exit(202)
            try:
                a, b = F2_algorithm(lb, ub, conf.f2_alg_confidence_radius,
                                    conf.f2_alg_error_prob)
                if a != -1:
                    pt = consumed_time()
                    print('F2: Z = {} * 2^{}  (processor time: {:.2f})'
                          ''.format(a, b, pt))
                else:
                    print('F2: Approximation Algorithm Failed!')
                    exit(203)
            except SolverBoundsExceededException:
                print('F2: Approximation Algorithm Failed!')
                exit(203)
        else:
            assert False, 'Unknown mode: {}'.format(conf.mode)
    except TotalTimeoutException:
        log.error('Total Timeout occured!')
        print('F2: The allowed running time expired before finish. Exiting.')
        exit(204)

    global SAT_SOLVER_TIME
    print("SAT_SOLVER_TIME:", SAT_SOLVER_TIME)

    global SAT_SOLVER_TIME_BETTER_PRECISION
    print("SAT_SOLVER_TIME_BETTER_PRECISION:", SAT_SOLVER_TIME_BETTER_PRECISION)

    exit(0)


def prepare_lower_bound():
    few_solutions = 2**4
    if len(conf.var_set) <= 10:  # Take care of really "small" cases
        t, b, n, ctime = sat_count(conf.formula_filename,
                                   conf.total_max_time - 2, few_solutions)
        if t:
            raise TotalTimeoutException
        elif n == 0:
            print('F2: UNSAT')
            exit(0)
        elif n < few_solutions:
            print('F2: The formula has exactly {} models'.format(n))
            exit(0)
        else:
            return math.log2(few_solutions)
    else:
        if is_unsat():
            print('F2: UNSAT')
            exit(0)
        else:
            return 0


def is_unsat():
    t, b, n, ctime = sat_count(conf.formula_filename,
                               conf.total_max_time - 2, 1)
    if t:
        raise TotalTimeoutException
    elif b:
        return False
    else:
        assert (n == 0)
        return True

class Bunch:
    # http://code.activestate.com/recipes/52308-the-simple-but-handy-collector-of-a-bunch-of-named/?in=user-97991
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


def find_lower_bound_call_from_python(problem_name, random_seed, var_degree, method='original', extra_configs=None, time_limit=500000):
    print("problem_name:", problem_name)

    #for create_biregular_adding_variables_orderByMarginals_randomInChunks
    #just keep a running sum of the satisfying solutions which can be divided by the number of solutions to get the marginals
    #benefit: avoid resumming stuff up
    #drawback: can't do something smarter like upweight higher quality samples or check that we got the same sample twice
    global USE_SUM_OF_SATISFYING_SOLUTIONS
    if method == 'bi_regular_order_vars_by_marginals_randomChunks':
        USE_SUM_OF_SATISFYING_SOLUTIONS = True
    else:
        USE_SUM_OF_SATISFYING_SOLUTIONS = False
    global SUM_OF_SATISFYING_SOLUTIONS
    SUM_OF_SATISFYING_SOLUTIONS = None

    global SATISFYING_SOLUTIONS_AS_LIST
    global DOUBLE_MARGINAL
    if method == 'bi_regular_order_vars_by_double_marginals':
        SATISFYING_SOLUTIONS_AS_LIST = True
        DOUBLE_MARGINAL = True
    else:
        SATISFYING_SOLUTIONS_AS_LIST = False
        DOUBLE_MARGINAL = False

    global SATISFYING_SOLUTIONS
    global SATISFYING_SOLUTIONS_ARRAY
    if SATISFYING_SOLUTIONS_AS_LIST:
        SATISFYING_SOLUTIONS = []
        SATISFYING_SOLUTIONS_ARRAY = None
    else:
        SATISFYING_SOLUTIONS = None

    global SAT_SOLVER_TIME
    SAT_SOLVER_TIME = 0  
    global SAT_SOLVER_TIME_BETTER_PRECISION
    SAT_SOLVER_TIME_BETTER_PRECISION = 0  
    global SAT_SOLVER_TIME_BETTER_PRECISION_PARALLEL
    SAT_SOLVER_TIME_BETTER_PRECISION_PARALLEL = 0  

    global CONSTRAINTS_CALLED_COUNT
    CONSTRAINTS_CALLED_COUNT = 0

    global SOLUTION_COUNT
    SOLUTION_COUNT = 0


    args = Bunch()

    args.input_filename = '/atlas/u/jkuck/approxmc/counting2/%s' % problem_name
    # args.input_filename = '/atlas/u/jkuck/low_density_parity_checks/SAT_problems_cnf/%s' % problem_name

    # args.working_dir = '/atlas/u/jkuck/F2/fireworks/augform'
    args.working_dir = os.getcwd()    
    args.keep_cnf = False
    args.random_seed = random_seed
    args.var_degree = var_degree
    args.cmsat_exe = '/atlas/u/jkuck/software/cryptominisat_BIRD/cryptominisat-5.6.6/build/cryptominisat5'
    # args.cmsat_exe = '/atlas/u/jkuck/XORModelCount/SATModelCount/cryptominisat'
    args.cmsat_version = 5
    # args.sharpsat_exe = 
    args.error_prob = .05
    args.max_time = time_limit

    args.mode = 'lb'
    args.skip_sharpsat = True

    args.sharpsat_exe = None
    args.lower_bound = math.nan
    args.upper_bound = math.nan
    args.confidence_radius = None
    args.lb_n_iter = math.nan
    args.ub_n_iter = math.nan
    args.extra_configs = extra_configs

    ##### other setup #####
    formula_lines = []  # The lines of the original cnf file
    global start_time
    start_time = time.process_time()
    global aug_counter
    aug_counter = 0  # Counter for augmented formulas filenames
    global subprocess_cumulative_time
    subprocess_cumulative_time = 0  # Accumulate the time spent in subprocesses    

    global conf
    print("1 conf.working_dir:", conf.working_dir)
    print("1 args.working_dir:", args.working_dir)
    conf = Config()
    conf.setup(args)
    print("2 conf.working_dir:", conf.working_dir)

    lb_init = prepare_lower_bound()
    lb, sat_solver_time, time_out = find_lower_bound(lb_init, method=method, return_time=True)
    print("remaining_time():", remaining_time())

    print("finished")
    print("problem_name:", problem_name, "var_degree:", var_degree, "method:", method)
    print("extra_configs['sum_of_T_solutions']:", extra_configs['sum_of_T_solutions'])
    return lb, sat_solver_time, time_out, SAT_SOLVER_TIME_BETTER_PRECISION_PARALLEL


def sharp_sat_call_from_python(problem_name, time_limit, problem_directory='/atlas/u/jkuck/approxmc/counting2/'):
    global SATISFYING_SOLUTIONS
    SATISFYING_SOLUTIONS = []
    global SAT_SOLVER_TIME
    SAT_SOLVER_TIME = 0  
    global SAT_SOLVER_TIME_BETTER_PRECISION
    SAT_SOLVER_TIME_BETTER_PRECISION = 0  
    global CONSTRAINTS_CALLED_COUNT
    CONSTRAINTS_CALLED_COUNT = 0

    args = Bunch()

    args.input_filename = '%s/%s' % (problem_directory, problem_name) 

    # args.input_filename = '/atlas/u/jkuck/low_density_parity_checks/SAT_problems_cnf/%s' % problem_name

    # args.working_dir = '/atlas/u/jkuck/F2/fireworks/augform'
    args.working_dir = os.getcwd()        
    args.keep_cnf = False
    args.random_seed = False
    args.var_degree = None
    args.cmsat_exe = '/atlas/u/jkuck/software/cryptominisat_BIRD/cryptominisat-5.6.6/build/cryptominisat5'
    # args.cmsat_exe = '/atlas/u/jkuck/XORModelCount/SATModelCount/cryptominisat'
    args.cmsat_version = 5
    # args.sharpsat_exe = 
    args.error_prob = .05
    args.max_time = time_limit

    args.mode = 'lb'
    args.skip_sharpsat = True

    args.sharpsat_exe = '/atlas/u/jkuck/sharpSAT/build/Release/sharpSAT'
    args.lower_bound = math.nan
    args.upper_bound = math.nan
    args.confidence_radius = None
    args.lb_n_iter = math.nan
    args.ub_n_iter = math.nan
    args.extra_configs = None

    ##### other setup #####
    formula_lines = []  # The lines of the original cnf file
    global start_time
    start_time = time.process_time()
    global aug_counter
    aug_counter = 0  # Counter for augmented formulas filenames
    global subprocess_cumulative_time
    subprocess_cumulative_time = 0  # Accumulate the time spent in subprocesses    

    global conf
    conf = Config()
    conf.setup(args)
    print("conf.working_dir:", conf.working_dir)

    #a triplet (time_out, solution_count, ctime) such that:
    #            time_out: [boolean] has the time_limit been reached without finishing?
    #            solution_count: number of solutions (if time_out is False)
    #            ct: The processor time consumed in the subprocess
    input_filename = '%s/%s' % (problem_directory, problem_name) 
    # input_filename = '/atlas/u/jkuck/approxmc/counting2/%s' % problem_name 
    # input_filename = '/atlas/u/jkuck/low_density_parity_checks/SAT_problems_cnf/%s' % problem_name

    print("input_filename:", input_filename)


    t0 = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime
    t0_bad = time.time()
    time_out, solution_count, ct = sharpsat_count(formula=input_filename, time_limit=time_limit)
    t1 = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime
    t1_bad = time.time()
    print("our timing of sharpSAT:", t1-t0)
    print("our (bad?) timing of sharpSAT:", t1_bad-t0_bad)
    sharp_sat_time = t1-t0
    return time_out, solution_count, sharp_sat_time


# @do_cprofile
def find_lower_bound(start=0, method='original', return_time=False):
    """
    Find a lower bound.

    PRECONDITION: The formula should not be UNSAT.

    :return: The binary logarithm of the lower bound

    Inputs:
    - method: (string)
        'original': use biregular constraints
        'iid': for comparing the lower bound, use iid constraints
    - return_time: (bool) return the time calling the SAT solver in addition to the bound 
    Outputs:
    - lb: (float) the lower bound
    - sat_solver_time: (float) the time used by the sat solver
    - time_out: (bool) True -> the method timed out
    """
    #this time matches cryptominisat output well, but has more precision
    global SAT_SOLVER_TIME_BETTER_PRECISION
    SAT_SOLVER_TIME_BETTER_PRECISION = 0

    #this time doesn't match cryptominisat well, but leaving it so code runs for now 
    global SAT_SOLVER_TIME
    SAT_SOLVER_TIME = 0

        
    print()
    print('-'*80)
    log.info('find_lower_bound(): phase 1')
    print("start =", start)
    print(sys.version)
    time1 = SAT_SOLVER_TIME_BETTER_PRECISION
    j = 0
    print("start + 2**j:", start + 2**j)
    (bound_holds, time_out) = algorithm1(int(start + 2**j), 1, sampleRepeats=10, method=method)
    if time_out:
        assert(bound_holds == False)
        return start, SAT_SOLVER_TIME_BETTER_PRECISION, time_out

    while bound_holds:
        j += 1        
        (bound_holds, time_out) = algorithm1(start + 2**j, 1, sampleRepeats=10, method=method)
        if time_out:
            assert(bound_holds == False)
            return start, SAT_SOLVER_TIME_BETTER_PRECISION, time_out

    if j == 0:
        return start, SAT_SOLVER_TIME_BETTER_PRECISION, time_out
    i = start + 2 ** (j - 1)
    up = start + 2 ** j - 1

    print("i:", i, "up:", up)

    time2 = SAT_SOLVER_TIME_BETTER_PRECISION
    print("phase 1 spent", time2-time1, "seconds calling the sat solver")
    print()
    print('-'*80)
    log.info('find_lower_bound(): phase 2  exponent={}'.format(j))
    global SATISFYING_SOLUTIONS
    # SATISFYING_SOLUTIONS = []
    while i < up:
        m = int(math.ceil((up - i)/2.0) + i)
        (bound_holds, time_out) = algorithm1(m, 1, sampleRepeats=10, method=method)
        if time_out:
            assert(bound_holds == False)
            return start, SAT_SOLVER_TIME_BETTER_PRECISION, time_out
        elif bound_holds:
            i = m
        else:
            up = m - 1
    time3 = SAT_SOLVER_TIME_BETTER_PRECISION
    print("phase 2 spent", time3-time2, "seconds calling the sat solver")

    print()
    print('-'*80)
    log.info('find_lower_bound(): phase 3  i={}'.format(i))

    if i > 2:
        j = 1
    else:
        j = 0
    while True:
        #jdk temp edit
        i -= 2 ** j
        j += 1

        if i <= start:
            break
        print("remaining time:", remaining_time())
        # (bound_holds, time_out) = algorithm1(i, 1, sampleRepeats=10, method=method)
        (bound_holds, time_out) = algorithm1(i, conf.alg1_n_iter, sampleRepeats=10, method=method)
        if time_out:
            assert(bound_holds == False)
            return start, SAT_SOLVER_TIME_BETTER_PRECISION, time_out
        elif bound_holds:
            break
        # i -= 2 ** j #jdk temp edit
        # j += 1 #jdk temp edit

    if i > start:
        lower_bound = i
    else:
        lower_bound = start
    log.info('find_lower_bound(): Returning {} (time elapsed: '
             '{:.2f})'.format(lower_bound, time.process_time()))
    time4 = SAT_SOLVER_TIME_BETTER_PRECISION
    print("phase 3 spent", time4-time3, "seconds calling the sat solver")

    if return_time:
        return lower_bound, SAT_SOLVER_TIME_BETTER_PRECISION, False
    # i=35
    # lower_bound=35

    ###############################################################################################################
    assert(False), "should ignore code below!!"

    print()
    print('-'*80)
    log.info('find_lower_bound(): phase 4, with samples i={}'.format(i))

    i = lower_bound + 1
    ret1 = lower_bound  
    (bound_holds, time_out) = algorithm1(i, conf.alg1_n_iter, sampleRepeats=1, method=method)
    if time_out:
        assert(bound_holds == False)
        return start, SAT_SOLVER_TIME_BETTER_PRECISION, time_out

    while bound_holds:
        lower_bound = i
        i += 1
        (bound_holds, time_out) = algorithm1(i, conf.alg1_n_iter, sampleRepeats=1, method=method)
        if time_out:
            assert(bound_holds == False)
            return start, SAT_SOLVER_TIME_BETTER_PRECISION, time_out

    log.info('after phase 4, find_lower_bound(): Returning {} (time elapsed: '
            '{:.2f})'.format(lower_bound, time.process_time()))
    time5 = SAT_SOLVER_TIME_BETTER_PRECISION
    print("phase 4, with samples spent", time5-time4, "seconds calling the sat solver")



    print()
    print('-'*80)
    log.info('find_lower_bound(): phase 4, no samples i={}'.format(i))

    # i = ret1 + 1
    # i = 171
    # while algorithm1(i, conf.alg1_n_iter, sampleRepeats=1):
    (bound_holds, time_out) = algorithm1(i, 240, sampleRepeats=1, method=method)
    if time_out:
        assert(bound_holds == False)
        return start, SAT_SOLVER_TIME_BETTER_PRECISION, time_out

    while bound_holds:
        ret1 = i
        i += 1
        (bound_holds, time_out) = algorithm1(i, 240, sampleRepeats=1, method=method)
        if time_out:
            assert(bound_holds == False)
            return start, SAT_SOLVER_TIME_BETTER_PRECISION, time_out


    log.info('after phase 4, find_lower_bound(): (would) Returning {} (time elapsed: '
             '{:.2f})'.format(ret1, time.process_time()))
    time6 = SAT_SOLVER_TIME_BETTER_PRECISION

    print("phase 4, no samples spent", time6-time5, "seconds calling the sat solver")



    return lower_bound


def find_upper_bound(lower):
    try:
        ub_a, ub_b = algorithm2(lower, conf.alg2_error_prob, conf.alg2_n_iter)
        return ub_a, ub_b
    except SolverBoundsExceededException:
        print('F2: Algorithm 2 Failed!')
        exit(201)


def extract_support_vars():
    """
    Generate the set of vars over which to generate the constraints.
    Extract the independent support or sampling set variables from header lines
    (starting with 'c ind') in the cnf file if they exist.
    Else return the set of all variables of the formula (as implied by the
    'p' header line of the DIMACS cnf file).

    Uses the global var formula_lines which contains the list of lines of
    the cnf file
    :return: A list of variables
    """
    global formula_lines  # The lines of the cnf file
    num_variables = -1
    indset = set()  # The vars of the independent or sampling set, if any
    for l in formula_lines:
        if l.lower().startswith('c ind'):
            for s in l.split(' ')[2:-1]:
                indset.add(int(s))
        if str(l.strip()[0:5]).lower() == 'p cnf':
            fields = l.strip().split(' ')
            num_variables = int(fields[2])
    if num_variables == -1:
        print('F2: Malformed input file. "P CNF" header line is missing!')
        exit(300)
    # if no 'c ind' variables are given, then all variables are presumed to
    # belong to the support
    if len(indset) == 0:
        indset = range(1, num_variables + 1)
    return list(indset), len(indset) > 0


def get_boost(constr):
    """
    Calculate the max boost for a list of constraint numbers
    :param constr: An int or a list of ints representing constraint numbers
    :return: Maximum value of a bound for boost for these constraint numbers
    """
    check_time()
    left_degree = conf.degn
    n = len(conf.var_set)
    if type(constr) == list:
        bs = [boost(left_degree, n, i) for i in constr]
        r = max(bs)
    else:
        assert type(constr) == int
        r = boost(left_degree, n, constr)
    log.info('The Boost is {:.2f}'.format(float(r)))
    return r


def generateConstraints(n_constr, degn, force_long=False, dummy_vars=0, sampleRepeats=1, method='original'):
    '''
    Inputs:
    - method: (string)
        'original': use biregular constraints
        'iid': for comparing the lower bound, use iid constraints

    '''
    if n_constr == 0:
        return list()

    if method == 'original':
        if needLong(n_constr) or force_long:
            constraints = generateLong(n_constr)
            return constraints, -1, -1, False #should really check there are no empty constraints, but unlikely
        else:
           # degn = 2.5
            constraints = generateLDPC(n_constr, dummy_vars, degn)
            return constraints, -1, -1, False

    # if needLong(n_constr) or force_long:
    #     return generateLong(n_constr), False #should really check there are no empty constraints, but unlikely
    # else:
    #     print("f=degn/n_constr:", degn/n_constr)
    #     # degn = 2.5
    #     return generateLDPC_bipartiteGraph_parameterizeF_pickEvenSampleSplit(n_constr, f=degn/n_constr), False


    # return generateLDPC_bipartiteGraph_parameterizeF(n_constr, f=.021), False


    #WAS USING THIS RECENTLY
    #currently, sample a biregular matrix and keep constraints with good marginals.  resample remaining constraints
    #sample 'problem constraints' iid (can also keep with bad marginals)
    elif method in ['bi_regular_marginals_per_constraint', 'bi_regular_order_vars_by_marginals', 'bi_regular_order_vars_by_double_marginals', 'bi_regular_order_vars_by_marginals_randomChunks']:
        if degn/n_constr > .4:
            print("generating long constraint")
            constraints = generateLong(n_constr)
            # constraints = generateShort(n_constr, f=.05)
            # constraints, fail_apriori, max_product_of_marginals = generate_iid(n_constr, dummy_vars, degn, clause_sample_repeats=1, f=.05)
            return constraints, -1, -1, False #should really check there are no empty constraints, but unlikely
        else:
            print("f =", degn/n_constr)
        global SATISFYING_SOLUTIONS   
        global SOLUTION_COUNT 
        if SOLUTION_COUNT < 10:
            constraints, max_product_of_marginals, i_effective = generateLDPC_bipartiteGraph_parameterizeF_pickEvenSampleSplit(n_constr, f=degn/n_constr, repeats=sampleRepeats)
        else:
            constraints, max_product_of_marginals, i_effective = generateLDPC_bipartiteGraph_parameterizeF_pickEvenSampleSplit_perConstraint(n_constr, f=degn/n_constr, method=method, marginal_cutoff=.45)
        return constraints, max_product_of_marginals, i_effective, False

    elif method == 'bi_regular_marginals_joint_constraint':
        constraints, max_product_of_marginals, i_effective = generateLDPC_bipartiteGraph_parameterizeF_pickEvenSampleSplit(n_constr, f=degn/n_constr, repeats=sampleRepeats)
        return constraints, max_product_of_marginals, i_effective, False


    # return generateLDPC_bipartiteGraph_parameterizeF_proportionalMarginals(n_constr, f=.018, approx_marginals=approx_marginals_sat_grid_pbl_0015), False
 
    # constraints, fail_apriori, max_product_of_marginals = generateLDPC_ours(n_constr, dummy_vars, degn, clause_sample_repeats=15)
    # i_effective = -1
    # return constraints, max_product_of_marginals, i_effective, fail_apriori

    elif method == 'iid':
        global CONSTRAINTS_CALLED_COUNT
        CONSTRAINTS_CALLED_COUNT += 1
        # print("#"*80)
        # print("CONSTRAINTS_CALLED_COUNT:", CONSTRAINTS_CALLED_COUNT)
        constraints, fail_apriori, max_product_of_marginals = generate_iid(n_constr, dummy_vars, degn, clause_sample_repeats=1, f=conf.extra_configs['f'])
        i_effective = -1
        return constraints, max_product_of_marginals, i_effective, fail_apriori

    else:
        assert(False)

def needLong(n_constr):
    """
    Are long constraints needed?
    :param n_constr:
    :return:
    """
    assert n_constr > 0
    n = len(conf.var_set)
    degn = conf.degn
    constr_length = math.ceil((n * degn) / n_constr)
    if constr_length >= n / 2:
        return True
    else:
        return False


def need_dummy_var(l, n, i):
    """
    Do we need an extra "dummy" var in the constraint creation
    :param l:         Left (Variable) degree
    :param n:         Num of variables
    :param i:         Num of constraints
    :return: True/False
    """
    return l * n / i == rfloor(l, n, i) and rfloor(l, n, i) % 2 == 0


def generateLong(n_constr):
    """
    Generate random long XOR constraints of average length n/2
    The independent set if any are taken from the
    configuration.
    :param n_constr: The number of constraints to generate
    :return: A list of constraints repr. as lists of strings ending with nl
    """
    assert n_constr > 0
    n = len(conf.var_set)
    var_list = conf.var_set
    i = n_constr
    log.info('Generating Long: n_constr={}  n_var={}'.format(n_constr, n))
    clauses_l = []
    for j in range(i):
        clause_varlist = [v for v in var_list if random.randint(0, 1) > 0]
        is_negated = (random.randint(0, 1) > 0)
        if is_negated:
            out = ['x -']
        else:
            out = ['x ']
        for v in clause_varlist:
            out.append('{} '.format(v))
        out.append('0\n')
        cl = ''.join(out)
        clauses_l.append(cl)
    return clauses_l

def generateShort(n_constr, f):
    """
    Generate random long XOR constraints of average length n/2
    The independent set if any are taken from the
    configuration.
    :param n_constr: The number of constraints to generate
    :return: A list of constraints repr. as lists of strings ending with nl
    """
    assert n_constr > 0
    n = len(conf.var_set)
    var_list = conf.var_set
    i = n_constr
    log.info('Generating Long: n_constr={}  n_var={}'.format(n_constr, n))
    clauses_l = []
    for j in range(i):
        clause_varlist = [v for v in var_list if np.random.rand() < f]
        is_negated = (random.randint(0, 1) > 0)
        if is_negated:
            out = ['x -']
        else:
            out = ['x ']
        for v in clause_varlist:
            out.append('{} '.format(v))
        out.append('0\n')
        cl = ''.join(out)
        clauses_l.append(cl)
    return clauses_l

def generateLDPC(n_constr, extra_vars, degn):
    """
    Generate the constraints corresponding to an LDPC code.
    The variable set is taken from the configuration. If extra_vars > 0, a
    number of "dummy" vars participate in the generation of the constraints.

    :param n_constr:    The number of constraints to generate
    :param extra_vars:  Number of extra variables
    :param degn:        Average variable (left) degree

    :return: A list of constraints repr. as lists of strings ending with nl
    """
    assert n_constr > 0
    assert extra_vars >= 0
    n0 = len(conf.var_set)
    if extra_vars == 0:
        var_list = conf.var_set
    else:
        var_list = conf.var_set.copy()
        var_list.extend(range(n0 + 1, n0 + extra_vars + 1))
    n = len(var_list)
    i = int(n_constr)
    x = n * degn / i  # Average constraint degree
    log.info('Generating LDPC: n_constr={}  n_var={} var_degree={}  '
             'constr_degree={:.1f}'.format(n_constr, n, degn, x))
    i_degseq = [math.floor(x)] * i
    n_degseq = [math.floor(degn)] * n
    if degn > math.floor(degn):
        x_ = math.ceil((degn - math.floor(degn)) * n)
        indices_ = random.sample(range(n), x_)
        for j_ in indices_:
            n_degseq[j_] += 1
    surplus = sum(n_degseq) - sum(i_degseq)
    assert 0 <= surplus <= n, 'Improper surplus  = ' + str(surplus)
    for c_ in random.sample(range(i), surplus):
        i_degseq[c_] += 1
    # print("i_degseq:")
    # print(i_degseq)
    # print("n_degseq:")
    # print(n_degseq)

    g = bipartite.configuration_model(i_degseq, n_degseq,
                                      create_using=nx.Graph())
    clauses_l = []
    for j in range(i):
        clause_varlist = [var_list[ind_ - i] for ind_ in g.neighbors(j)]
        is_negated = (random.randint(0, 1) > 0)
        if is_negated:
            out = ['x -']
        else:
            out = ['x ']
        for v in clause_varlist:
            out.append('{} '.format(v))
        out.append('0\n')
        cl = ''.join(out)
        clauses_l.append(cl)
    return clauses_l


def get_clause_marginal(satisfying_solutions, clause):
    '''
    Inputs:
    - satisfying_solutions: (np.array with shape (#solutions x #variables)) satisfying solutions to the sat problem
    - clause: (np.array) a parity clause

    Outputs:
    - marginal: (float) the fraction of satsifying solutions that have an even number of 1s intersecting the clause
    '''
    # assert(clause.shape[0] == satisfying_solutions.shape[1])
    global SATISFYING_SOLUTIONS_AS_LIST
    if SATISFYING_SOLUTIONS_AS_LIST:
        if len(satisfying_solutions) == 0:
            return 0
        marginal = 0
        for satisfying_solution in satisfying_solutions:
            if np.dot(satisfying_solution, clause) % 2 == 1:
                marginal += 1
        # marginal = np.sum(np.mod(np.dot(satisfying_solutions[:100],clause), 2))
        marginal /= len(satisfying_solutions)
        assert(marginal >=0 and marginal <= 1.0)
        symmetric_marginal = min(marginal, 1-marginal)
        # return symmetric_marginal, marginal


###########TESTING############

        global SATISFYING_SOLUTIONS_ARRAY
        # marginal1 = np.sum(np.mod(np.dot(SATISFYING_SOLUTIONS_ARRAY[:SOLUTION_COUNT],clause), 2))
        marginal1 = np.sum(np.mod(np.dot(SATISFYING_SOLUTIONS_ARRAY[:SOLUTION_COUNT],clause), 2))
        marginal1 /= SOLUTION_COUNT
        assert(marginal1 >=0 and marginal1 <= 1.0)
        symmetric_marginal1 = min(marginal1, 1-marginal1)
        assert(marginal == marginal1), (SOLUTION_COUNT, marginal, marginal1, satisfying_solutions[0], np.sum(SATISFYING_SOLUTIONS_ARRAY[:10, :], axis=1), clause)
        assert(symmetric_marginal == symmetric_marginal1), (marginal, marginal1, satisfying_solutions[0], clause)
###########TESTING############
        return symmetric_marginal, marginal
    else:
        global SOLUTION_COUNT
        # satisfying_solutions_array = np.array(satisfying_solutions)
        # marginal1 = np.sum(np.mod(np.dot(satisfying_solutions_array[:SOLUTION_COUNT],clause), 2))
        marginal1 = np.sum(np.mod(np.dot(satisfying_solutions[:SOLUTION_COUNT],clause), 2))
        marginal1 /= SOLUTION_COUNT
        assert(marginal1 >=0 and marginal1 <= 1.0)
        symmetric_marginal = min(marginal1, 1-marginal1)
        # assert(marginal == marginal1), (marginal, marginal1, satisfying_solutions[0], clause)

        return symmetric_marginal, marginal1


def get_clause_double_marginal(satisfying_solutions, clause, new_var_idx):
    '''
    Inputs:
    - satisfying_solutions: (np.array with shape (#solutions x #variables)) satisfying solutions to the sat problem
    - clause: (np.array) a parity clause

    Outputs:
    - marginal: (float) the fraction of satsifying solutions that have an even number of 1s intersecting the clause
    '''
    # assert(clause.shape[0] == satisfying_solutions.shape[1])
    # print("get_clause_double_marginal called")
    global SATISFYING_SOLUTIONS_AS_LIST
    # if SATISFYING_SOLUTIONS_AS_LIST:
    if True:
        if len(satisfying_solutions) == 0:
            return 0
        marginal1 = 0
        marginal0 = 0

        for satisfying_solution in satisfying_solutions:
            if np.dot(satisfying_solution, clause) % 2 == 1:
                marginal1 += satisfying_solution[new_var_idx]
            else:
                marginal0 += satisfying_solution[new_var_idx]
        # marginal1 = np.sum(np.mod(np.dot(satisfying_solutions[:100],clause), 2))
        marginal1 /= len(satisfying_solutions)
        marginal0 /= len(satisfying_solutions)
        assert(marginal0 >=0 and marginal0 <= 1.0)
        assert(marginal1 >=0 and marginal1 <= 1.0)
        symmetric_marginal1 = min(marginal1, 1-marginal1)
        symmetric_marginal0 = min(marginal0, 1-marginal0)
        double_marginal = symmetric_marginal1*symmetric_marginal0
        return double_marginal, None

    else:
        global SOLUTION_COUNT
        # satisfying_solutions_array = np.array(satisfying_solutions)
        # marginal1 = np.sum(np.mod(np.dot(satisfying_solutions_array[:SOLUTION_COUNT],clause), 2))
        marginal1 = np.sum(np.mod(np.dot(satisfying_solutions[:SOLUTION_COUNT],clause), 2))
        marginal1 /= SOLUTION_COUNT
        assert(marginal1 >=0 and marginal1 <= 1.0)
        symmetric_marginal = min(marginal1, 1-marginal1)
        # assert(marginal == marginal1), (marginal, marginal1, satisfying_solutions[0], clause)

        return symmetric_marginal, marginal1


def count_hash_collisions(satisfying_solutions, constraint_matrix):
    sum_of_unique_hash_bins = 0.0
    sum_of_max_bin_counts = 0.0
    num_constraints_to_group = 3
    product_of_min_fractions = 1.0
    product_of_variances = 1.0
    for idx in range(np.int(np.floor(constraint_matrix.shape[0]/num_constraints_to_group))): 
        if idx + num_constraints_to_group >= constraint_matrix.shape[0]:
            break
        hash_bins = defaultdict(int)
        for satisfying_solution in satisfying_solutions:
            hash_bins[tuple(np.mod(np.dot(constraint_matrix[num_constraints_to_group*idx:num_constraints_to_group*(idx+1)], satisfying_solution), 2))] += 1
        sum_of_unique_hash_bins += len(hash_bins)
        sum_of_max_bin_counts += max(hash_bins.values())
        hash_bins_fractions = np.zeros(2**num_constraints_to_group)
        hash_bin_idx = 0
        for bin_label, bin_count in hash_bins.items():
            hash_bins_fractions[hash_bin_idx] = bin_count/(len(satisfying_solutions))
            hash_bin_idx += 1
        product_of_min_fractions *= min(hash_bins_fractions)
        product_of_variances *= np.var(hash_bins_fractions)
        # print(product_of_variances, end=' ')
        # print("product of variance:", product_of_variances, "current variance:", np.var(hash_bins_fractions), hash_bins_fractions)
        # print("hash_bins_fractions:", hash_bins_fractions)
    # print('*'*80)
    # print('*'*80)
    # print("unique bins hashed to =", sum_of_unique_hash_bins, 'of max:', satisfying_solutions.shape[0])
    # print("unique bins hashed to =", sum_of_unique_hash_bins, 'of max:', len(satisfying_solutions))
    # print(hash_bins)
    # return product_of_min_fractions
    return product_of_variances

def generateLDPC_bipartiteGraph_parameterizeF_pickEvenSampleSplit(n_constr, f, repeats=1):
    """
    Generate the constraints corresponding to an LDPC code.
    The variable set is taken from the configuration. If extra_vars > 0, a
    number of "dummy" vars participate in the generation of the constraints.

    :param n_constr:    The number of constraints to generate
    :param f:  density of ones

    :return: A list of constraints repr. as lists of strings ending with nl
    """
    max_product_of_marginals = -np.inf
    clauses_with_max_product_of_marginals = None
    i_effective_for_max_product_of_marginals = None
    global SATISFYING_SOLUTIONS

    if VERBOSE:
        print('f=', f, 'var degree=', n_constr*f, 'clause degree=', len(conf.var_set)*f)

    for repeat in range(repeats):
        assert n_constr > 0
        n0 = len(conf.var_set)
        var_list = conf.var_set
        n = len(var_list)
        i = int(n_constr)
        degn = i*f # average variable degree (ones per column)
        x = n * f  # Average constraint degree
        # log.info('Generating LDPC: n_constr={}  n_var={} var_degree={}  '
        #          'constr_degree={:.1f}'.format(n_constr, n, degn, x))
        i_degseq = [math.floor(x)] * i
        n_degseq = [math.floor(degn)] * n
        if degn > math.floor(degn):
            x_ = math.ceil((degn - math.floor(degn)) * n)
            indices_ = random.sample(range(n), x_)
            for j_ in indices_:
                n_degseq[j_] += 1
        surplus = sum(n_degseq) - sum(i_degseq)
        assert 0 <= surplus <= n, 'Improper surplus  = ' + str(surplus)
        for c_ in random.sample(range(i), surplus):
            i_degseq[c_] += 1
        # print("i_degseq:")
        # print(i_degseq)
        # print("n_degseq:")
        # print(n_degseq)
        
        g = bipartite.configuration_model(i_degseq, n_degseq,
                                          create_using=nx.Graph())
        


        clauses_l = []
        symmetric_marginals = []
        marginals = []
        constraint_list = []
        cur_i_effective = 1.0
        global USE_SUM_OF_SATISFYING_SOLUTIONS
        for j in range(i):
            clause_varlist = [var_list[ind_ - i] for ind_ in g.neighbors(j)]
            clause_independentSet_indices = [ind_ - i for ind_ in g.neighbors(j)]
            PRINT_MARGINAL = True
            if PRINT_MARGINAL:
                clause_array = np.zeros(n)
                for v in clause_independentSet_indices:
                    assert(v >= 0 and v <= n-1), (v, n)
                    clause_array[v] = 1.0
                constraint_list.append(clause_array)
                if not USE_SUM_OF_SATISFYING_SOLUTIONS:
                    cur_symmetric_marginal, cur_marginal = get_clause_marginal(SATISFYING_SOLUTIONS, clause_array)
                    # print("clause marginal:", cur_symmetric_marginal)
                    symmetric_marginals.append(cur_symmetric_marginal)
                    marginals.append(cur_marginal)
            # is_negated = (random.randint(0, 1) > 0)
            cur_marginal = .5
            if np.random.rand() < cur_marginal:
                is_negated = False
                cur_i_effective /= cur_marginal
            else:
                is_negated = True
                cur_i_effective /= (1-cur_marginal)
            if is_negated:
                out = ['x -']
            else:
                out = ['x ']
            for v in clause_varlist:
                out.append('{} '.format(v))
            out.append('0\n')
            cl = ''.join(out)
            clauses_l.append(cl)

        constraint_matrix = np.array(constraint_list)
        # product_of_min_fractions = count_hash_collisions(SATISFYING_SOLUTIONS, constraint_matrix)
        if not USE_SUM_OF_SATISFYING_SOLUTIONS:
            product_of_variances = count_hash_collisions(SATISFYING_SOLUTIONS, constraint_matrix)

        # symmetric_marginals.sort()
            symmetric_marginals, constraint_list = [list(x) for x in zip(*sorted(zip(symmetric_marginals, constraint_list), key=lambda pair: pair[0]))]

        # print("all symmetric_marginals", symmetric_marginals)
        # print("smallest 10 symmetric_marginals", symmetric_marginals[:10])
        # print("largest 10 symmetric_marginals", symmetric_marginals[-10:])

        # if len(constraint_list) > 10:
        #     top_10_clauses_combined = np.zeros(n)
        #     bottom_10_clauses_combined = np.zeros(n)
        #     for var_idx in range(n):
        #         top_set_to_1 = False
        #         bottom_set_to_1 = False
        #         for clause_idx in range(10):
        #             if constraint_list[clause_idx][var_idx] == 1.0:
        #                 bottom_set_to_1 = True
        #             if constraint_list[-clause_idx][var_idx] == 1.0:
        #                 top_set_to_1 = True
        #         if top_set_to_1:
        #             top_10_clauses_combined[var_idx] = 1.0
        #         if bottom_set_to_1:
        #             bottom_10_clauses_combined[var_idx] = 1.0

        #     print("marginal for combined top 10 marginal clauses:", get_clause_marginal(SATISFYING_SOLUTIONS, top_10_clauses_combined))
        #     print("marginal for combined bottom 10 marginal clauses:", get_clause_marginal(SATISFYING_SOLUTIONS, bottom_10_clauses_combined))
        #     print()

        if np.product(symmetric_marginals) > max_product_of_marginals:
        # if repeat == 0:
            # print("symmetric_marginals:", symmetric_marginals)

            max_product_of_marginals = np.product(symmetric_marginals)
            clauses_with_max_product_of_marginals = clauses_l
            i_effective_for_max_product_of_marginals = cur_i_effective
            # break
  

    print("i_effective_for_max_product_of_marginals =", np.log(i_effective_for_max_product_of_marginals)/np.log(2))
    if VERBOSE:
        pass
    print("max_product_of_marginals:", max_product_of_marginals, end=' ')
    return clauses_with_max_product_of_marginals, max_product_of_marginals, np.log(i_effective_for_max_product_of_marginals)/np.log(2)

# @do_cprofile
def create_biregular_adding_variables_orderByMarginals_randomInChunks(n, m, total_ones):
    '''
    Inputs:
    - n: (int) the number of variables
    - m: (int) the number of constraints
    - total_ones: (int) the number of ones in the constraint matrix
    Create a biregular constraint matrix.  
    1. find approximate marginals for each variable
    2. add the varables from best marginals to worst (closest to .5 being the best), looping around until the correct density is reached
    '''
    t0 = time.time()
    print("hi", "n:", n, "m:", m, "total_ones:", total_ones)
    global SATISFYING_SOLUTIONS
    global SATISFYING_SOLUTIONS_AS_LIST
    global SOLUTION_COUNT
    global SUM_OF_SATISFYING_SOLUTIONS
    global USE_SUM_OF_SATISFYING_SOLUTIONS

    if USE_SUM_OF_SATISFYING_SOLUTIONS:
        approximate_marginals = SUM_OF_SATISFYING_SOLUTIONS/SOLUTION_COUNT
        symmetric_marginals = np.zeros(n)
        for idx, marginal in enumerate(approximate_marginals):
            symmetric_marginals[idx] = min(approximate_marginals[idx], 1-approximate_marginals[idx])

    else:
        if SATISFYING_SOLUTIONS_AS_LIST:
            assert(len(SATISFYING_SOLUTIONS) > 2) #otherwise sample without marginals
            # compute variable marginals
            approximate_marginals = np.zeros(n)
            for satisfying_solution in SATISFYING_SOLUTIONS:
                approximate_marginals += satisfying_solution
            approximate_marginals /= len(SATISFYING_SOLUTIONS)
            symmetric_marginals = np.zeros(n)
            for idx, marginal in enumerate(approximate_marginals):
                symmetric_marginals[idx] = min(approximate_marginals[idx], 1-approximate_marginals[idx])

        else:
            assert(SATISFYING_SOLUTIONS.shape[0] > 2) #otherwise sample without marginals        
            # compute variable marginals
            approximate_marginals = np.sum(SATISFYING_SOLUTIONS[:SOLUTION_COUNT], axis = 0)
            assert(len(approximate_marginals) == n)
            approximate_marginals /= SOLUTION_COUNT
            symmetric_marginals = np.zeros(n)
            for idx, marginal in enumerate(approximate_marginals):
                symmetric_marginals[idx] = min(approximate_marginals[idx], 1-approximate_marginals[idx])


    #variable_order[0] is the variable with closest marginal to .5 (given by an integer, [0,n-1])
    #variable_order[1] is the variable with 2nd closest marginal to .5 (given by an integer, [0,n-1])
    #...
    #variable_order[n-1] is the variable with marginal furthest from .5 (given by an integer, [0,n-1])
    variable_order = np.argsort(-symmetric_marginals)

    MAKE_VARIANCE_HIGH_FOR_UPPER_BOUND = False
    if MAKE_VARIANCE_HIGH_FOR_UPPER_BOUND:
        ub_variable_order = [0 for i in range(n)]

        block_count = math.ceil(n/m)
        if block_count*m == n:
            small_bin_size = 0
        else:
            small_bin_size = n - block_count*m

        block_idx = 0
        bucket_idx = 0
        for var in variable_order:
            new_idx = block_idx*m + bucket_idx
            ub_variable_order[new_idx] = var
            block_idx += 1
            if block_idx == block_count or block_idx*m + bucket_idx >= n:
                block_idx = 0
                bucket_idx += 1
        zero_count = 0
        for var in ub_variable_order:
            if var == 0:
                zero_count += 1
        assert(zero_count == 1)
        variable_order = ub_variable_order

    t1 = time.time()

    constraints = np.zeros((m,n))
    constraints_as_list = [[] for i in range(m)]
    timeA = 0
    timeB = 0
    timeB2 = 0
    remaining_ones = total_ones # the number of ones left to allocate
    variable_index = 0 # allocate variables variable_order[variable_index:variable_index+ones_to_allocate]
    #at each iteration allocate m ones, one to each constraint.
    #except for the last iteration, where we allocate however many ones remain (<= m)
    while remaining_ones > 0: 
        # print("remaining_ones:", remaining_ones)
        if remaining_ones < m:
            ones_to_allocate = remaining_ones
            remaining_ones = 0
        else:
            ones_to_allocate = m
            remaining_ones -= m
        if variable_index >= n:
            variable_index = variable_index % n
        if variable_index+ones_to_allocate < n:
            variables = variable_order[variable_index:variable_index+ones_to_allocate]
        else:
            ones_at_end = n - variable_index
            ones_at_beginning = ones_to_allocate - ones_at_end
            variables = list(variable_order[variable_index:])
            variables_at_begining = list(variable_order[:ones_at_beginning])
            variables.extend(variables_at_begining)
        assert(len(variables) == ones_to_allocate and len(variable_order) == n)

        variable_index += ones_to_allocate

        valid_shuffling = False #we can't add a variable to a constraint that it is already in
        tries_to_get_valid_assignment=1
        while not valid_shuffling:
            permuted_constraint_indices = [i for i in range(m)]
            np.random.shuffle(permuted_constraint_indices)
            valid_shuffling = True
            t1_a = time.time()
            for variable_idx, constraint_idx in enumerate(permuted_constraint_indices[:ones_to_allocate]):
                assert(constraint_idx <= m)
                # if constraints[constraint_idx, variables[variable_idx]] == 1:
                if variables[variable_idx] in constraints_as_list[constraint_idx]:
                    valid_shuffling = False
                    tries_to_get_valid_assignment += 1
                    # assert(variables[variable_idx] in constraints_as_list[constraint_idx])
                    break
                # else:
                #     assert(variables[variable_idx] not in constraints_as_list[constraint_idx])
            t1_b = time.time()
            timeA += t1_b-t1_a
        if tries_to_get_valid_assignment > 1:
            print("we need", tries_to_get_valid_assignment, "tries to get a valid assignment")

        t1_c = time.time()
        for (variable_idx, constraint_idx) in enumerate(permuted_constraint_indices[:ones_to_allocate]):
            # assert(constraints[constraint_idx, variables[variable_idx]] == 0), (association_list, constraints[constraint_idx, variable_idx], constraints, constraint_idx, variable_idx, cost_matrix)
            # constraints[constraint_idx, variables[variable_idx]] = 1
            t1_c1 = time.time()
            assert(variables[variable_idx] not in constraints_as_list[constraint_idx])
            t1_c2 = time.time()            
            constraints_as_list[constraint_idx].append(variables[variable_idx])
        t1_d = time.time()
        timeB += t1_d - t1_c
        timeB2 += t1_c2 - t1_c1

    t2 = time.time()

    # return constraints
    clauses = []
    list_of_marginals = []
    product_of_marginals = 1.0
    # for constr_idx in range(m):
    #     cur_clause = []
    #     for var_idx in range(n):
    #         if constraints[constr_idx, var_idx] == 1:
    #             cur_clause.append(var_idx)
    #     clauses.append(cur_clause)
    #     if not USE_SUM_OF_SATISFYING_SOLUTIONS:
    #         cur_symmetric_marginal, cur_marginal = get_clause_marginal(SATISFYING_SOLUTIONS, constraints[constr_idx])
    #         list_of_marginals.append(cur_symmetric_marginal)
    #         product_of_marginals *= cur_symmetric_marginal

    # list_of_marginals.sort()

    t3 = time.time()

    # assert(len(clauses) == len(constraints_as_list))
    # for clause_idx in range(len(clauses)):
    #     assert(set(clauses[clause_idx]) == set(constraints_as_list[clause_idx]))
    print()
    print()
    print()
    print()
    print('!'*80)
    print("t3-t0", t3-t0)
    print("t1-t0", t1-t0)
    print("t2-t1", t2-t1)
    print("t3-t2", t3-t2)

    print("timeA:", timeA)
    print("timeB:", timeB)
    print("timeB2:", timeB2)
    # global CONSTRAINT_GENERATION_TIME
    # CONSTRAINT_GENERATION_TIME += t3 - t0 


    # return clauses, product_of_marginals, list_of_marginals
    return constraints_as_list, product_of_marginals, list_of_marginals


# @do_cprofile
def create_biregular_adding_variables_orderByMarginals_randomInChunks_faster(n, m, total_ones):
    '''
    Inputs:
    - n: (int) the number of variables
    - m: (int) the number of constraints
    - total_ones: (int) the number of ones in the constraint matrix
    Create a biregular constraint matrix.  
    1. find approximate marginals for each variable
    2. add the varables from best marginals to worst (closest to .5 being the best), looping around until the correct density is reached
    '''
    t0 = time.time()
    print("hi", "n:", n, "m:", m, "total_ones:", total_ones)
    global SATISFYING_SOLUTIONS
    global SATISFYING_SOLUTIONS_AS_LIST
    global SOLUTION_COUNT
    global SUM_OF_SATISFYING_SOLUTIONS
    global USE_SUM_OF_SATISFYING_SOLUTIONS

    if USE_SUM_OF_SATISFYING_SOLUTIONS:
        approximate_marginals = SUM_OF_SATISFYING_SOLUTIONS/SOLUTION_COUNT
        symmetric_marginals = np.zeros(n)
        for idx, marginal in enumerate(approximate_marginals):
            symmetric_marginals[idx] = min(approximate_marginals[idx], 1-approximate_marginals[idx])

    else:
        if SATISFYING_SOLUTIONS_AS_LIST:
            assert(len(SATISFYING_SOLUTIONS) > 2) #otherwise sample without marginals
            # compute variable marginals
            approximate_marginals = np.zeros(n)
            for satisfying_solution in SATISFYING_SOLUTIONS:
                approximate_marginals += satisfying_solution
            approximate_marginals /= len(SATISFYING_SOLUTIONS)
            symmetric_marginals = np.zeros(n)
            for idx, marginal in enumerate(approximate_marginals):
                symmetric_marginals[idx] = min(approximate_marginals[idx], 1-approximate_marginals[idx])

        else:
            assert(SATISFYING_SOLUTIONS.shape[0] > 2) #otherwise sample without marginals        
            # compute variable marginals
            approximate_marginals = np.sum(SATISFYING_SOLUTIONS[:SOLUTION_COUNT], axis = 0)
            assert(len(approximate_marginals) == n)
            approximate_marginals /= SOLUTION_COUNT
            symmetric_marginals = np.zeros(n)
            for idx, marginal in enumerate(approximate_marginals):
                symmetric_marginals[idx] = min(approximate_marginals[idx], 1-approximate_marginals[idx])


    #variable_order[0] is the variable with closest marginal to .5 (given by an integer, [0,n-1])
    #variable_order[1] is the variable with 2nd closest marginal to .5 (given by an integer, [0,n-1])
    #...
    #variable_order[n-1] is the variable with marginal furthest from .5 (given by an integer, [0,n-1])
    variable_order = np.argsort(-symmetric_marginals)

    MAKE_VARIANCE_HIGH_FOR_UPPER_BOUND = False
    if MAKE_VARIANCE_HIGH_FOR_UPPER_BOUND:
        ub_variable_order = [0 for i in range(n)]

        block_count = math.ceil(n/m)
        if block_count*m == n:
            small_bin_size = 0
        else:
            small_bin_size = n - block_count*m

        block_idx = 0
        bucket_idx = 0
        for var in variable_order:
            new_idx = block_idx*m + bucket_idx
            ub_variable_order[new_idx] = var
            block_idx += 1
            if block_idx == block_count or block_idx*m + bucket_idx >= n:
                block_idx = 0
                bucket_idx += 1
        zero_count = 0
        for var in ub_variable_order:
            if var == 0:
                zero_count += 1
        assert(zero_count == 1)
        variable_order = ub_variable_order

    t1 = time.time()

    dense_constraint_array = -1*np.zeros((m,int(np.ceil(total_ones/m))), dtype=np.int)
    timeA = 0
    timeB = 0
    timeB2 = 0
    remaining_ones = total_ones # the number of ones left to allocate
    variable_index = 0 # allocate variables variable_order[variable_index:variable_index+ones_to_allocate]
    #at each iteration allocate m ones, one to each constraint.
    #except for the last iteration, where we allocate however many ones remain (<= m)
    constraint_variable_idx = 0
    while remaining_ones > 0: 
        # print("remaining_ones:", remaining_ones)
        if remaining_ones < m:
            ones_to_allocate = remaining_ones
            remaining_ones = 0
        else:
            ones_to_allocate = m
            remaining_ones -= m
        if variable_index >= n:
            variable_index = variable_index % n
        if variable_index+ones_to_allocate < n:
            variables = variable_order[variable_index:variable_index+ones_to_allocate]
        else:
            ones_at_end = n - variable_index
            ones_at_beginning = ones_to_allocate - ones_at_end
            variables = list(variable_order[variable_index:])
            variables_at_begining = list(variable_order[:ones_at_beginning])
            variables.extend(variables_at_begining)
        assert(len(variables) == ones_to_allocate and len(variable_order) == n)

        variable_index += ones_to_allocate

        valid_shuffling = False #we can't add a variable to a constraint that it is already in
        tries_to_get_valid_assignment=1
        while not valid_shuffling:
            permuted_constraint_indices = [i for i in range(m)]
            np.random.shuffle(permuted_constraint_indices)
            valid_shuffling = True
            t1_a = time.time()
            # for variable_idx, constraint_idx in enumerate(permuted_constraint_indices[:ones_to_allocate]):
            #     assert(constraint_idx <= m)
            #     # if constraints[constraint_idx, variables[variable_idx]] == 1:
            #     if variables[variable_idx] in constraints_as_list[constraint_idx]:
            #         valid_shuffling = False
            #         tries_to_get_valid_assignment += 1
            #         # assert(variables[variable_idx] in constraints_as_list[constraint_idx])
            #         break
            #     # else:
            #     #     assert(variables[variable_idx] not in constraints_as_list[constraint_idx])
            t1_b = time.time()
            timeA += t1_b-t1_a
        if tries_to_get_valid_assignment > 1:
            print("we need", tries_to_get_valid_assignment, "tries to get a valid assignment")

        t1_c = time.time()
        for (variable_idx, constraint_idx) in enumerate(permuted_constraint_indices[:ones_to_allocate]):
            # # assert(constraints[constraint_idx, variables[variable_idx]] == 0), (association_list, constraints[constraint_idx, variable_idx], constraints, constraint_idx, variable_idx, cost_matrix)
            # # constraints[constraint_idx, variables[variable_idx]] = 1
            # t1_c1 = time.time()
            # assert(variables[variable_idx] not in constraints_as_list[constraint_idx])
            # t1_c2 = time.time()            
            # constraints_as_list[constraint_idx].append(variables[variable_idx])
            dense_constraint_array[constraint_idx, constraint_variable_idx] = variables[variable_idx]
        constraint_variable_idx += 1

        t1_d = time.time()
        timeB += t1_d - t1_c
        # timeB2 += t1_c2 - t1_c1

    t2 = time.time()

    # return constraints
    clauses = []
    list_of_marginals = []
    product_of_marginals = 1.0
    # for constr_idx in range(m):
    #     cur_clause = []
    #     for var_idx in range(n):
    #         if constraints[constr_idx, var_idx] == 1:
    #             cur_clause.append(var_idx)
    #     clauses.append(cur_clause)
    #     if not USE_SUM_OF_SATISFYING_SOLUTIONS:
    #         cur_symmetric_marginal, cur_marginal = get_clause_marginal(SATISFYING_SOLUTIONS, constraints[constr_idx])
    #         list_of_marginals.append(cur_symmetric_marginal)
    #         product_of_marginals *= cur_symmetric_marginal

    # list_of_marginals.sort()

    # constraints_as_list = []
    # for constraint_idx in range(m):
    #     constraints_as_list.append([i for i in dense_constraint_array[constraint_idx]])
    #     if constraints_as_list[constraint_idx][-1] == -1:
    #         del(constraints_as_list[constraint_idx][-1])

    t3 = time.time()

    # assert(len(clauses) == len(constraints_as_list))
    # for clause_idx in range(len(clauses)):
    #     assert(set(clauses[clause_idx]) == set(constraints_as_list[clause_idx]))
    print()
    print()
    print()
    print()
    print('!'*80)
    print("t3-t0", t3-t0)
    print("t1-t0", t1-t0)
    print("t2-t1", t2-t1)
    print("t3-t2", t3-t2)

    print("timeA:", timeA)
    print("timeB:", timeB)
    print("timeB2:", timeB2)
    # global CONSTRAINT_GENERATION_TIME
    # CONSTRAINT_GENERATION_TIME += t3 - t0 


    # return clauses, product_of_marginals, list_of_marginals
    # return constraints_as_list, product_of_marginals, list_of_marginals
    return dense_constraint_array, product_of_marginals, list_of_marginals


# def lapjv():
#     return 5

def non_square_lapjv1(cost_matrix):
    assert(cost_matrix.shape[1] >= cost_matrix.shape[0])
    if cost_matrix.shape[1] == cost_matrix.shape[0]:
        cost, col_ind, junk = lap.lapjv(cost_matrix)
        association_list = zip(range(len(col_ind)), col_ind)
        association_list = [assoc for assoc in association_list]
        return association_list
    else:
        row_count = cost_matrix.shape[0]
        square_cost_matrix = np.zeros((cost_matrix.shape[1],cost_matrix.shape[1]))
        square_cost_matrix[:row_count, :] = cost_matrix
        cost, col_ind, junk = lap.lapjv(square_cost_matrix)                
        association_list = zip(range(row_count), col_ind[:row_count])
        association_list = [assoc for assoc in association_list]
        return association_list

def non_square_lapjv(cost_matrix):
    assert(cost_matrix.shape[1] >= cost_matrix.shape[0])
    if cost_matrix.shape[1] == cost_matrix.shape[0]:
        col_ind, junk, junk1 = lapjv.lapjv(cost_matrix)
        association_list = zip(range(len(col_ind)), col_ind)
        association_list = [assoc for assoc in association_list]
        return association_list
    else:
        row_count = cost_matrix.shape[0]
        square_cost_matrix = np.zeros((cost_matrix.shape[1],cost_matrix.shape[1]))
        square_cost_matrix[:row_count, :] = cost_matrix
        col_ind, junk, junk1 = lapjv.lapjv(square_cost_matrix)                
        association_list = zip(range(row_count), col_ind[:row_count])
        association_list = [assoc for assoc in association_list]
        return association_list

def create_biregular_adding_variables_orderByMarginals_assignmentProblem(n, m, total_ones):
    '''
    Inputs:
    - n: (int) the number of variables
    - m: (int) the number of constraints
    - total_ones: (int) the number of ones in the constraint matrix
    Create a biregular constraint matrix.  
    1. find approximate marginals for each variable
    2. add the varables from best marginals to worst (closest to .5 being the best), looping around until the correct density is reached
    '''
    print("hi", "n:", n, "m:", m, "total_ones:", total_ones)
    global SATISFYING_SOLUTIONS
    global SATISFYING_SOLUTIONS_AS_LIST
    global SOLUTION_COUNT
    if SATISFYING_SOLUTIONS_AS_LIST:
        assert(len(SATISFYING_SOLUTIONS) > 2) #otherwise sample without marginals
        # compute variable marginals
        approximate_marginals = np.zeros(n)
        for satisfying_solution in SATISFYING_SOLUTIONS:
            approximate_marginals += satisfying_solution
        approximate_marginals /= len(SATISFYING_SOLUTIONS)
        symmetric_marginals = np.zeros(n)
        for idx, marginal in enumerate(approximate_marginals):
            symmetric_marginals[idx] = min(approximate_marginals[idx], 1-approximate_marginals[idx])

    else:
        assert(SATISFYING_SOLUTIONS.shape[0] > 2) #otherwise sample without marginals        
        # compute variable marginals
        approximate_marginals = np.sum(SATISFYING_SOLUTIONS[:SOLUTION_COUNT], axis = 0)
        assert(len(approximate_marginals) == n), (len(approximate_marginals), n, SATISFYING_SOLUTIONS.shape, approximate_marginals.shape, SOLUTION_COUNT)
        approximate_marginals /= SOLUTION_COUNT
        symmetric_marginals = np.zeros(n)
        for idx, marginal in enumerate(approximate_marginals):
            symmetric_marginals[idx] = min(approximate_marginals[idx], 1-approximate_marginals[idx])

    #variable_order[0] is the variable with closest marginal to .5 (given by an integer, [0,n-1])
    #variable_order[1] is the variable with 2nd closest marginal to .5 (given by an integer, [0,n-1])
    #...
    #variable_order[n-1] is the variable with marginal furthest from .5 (given by an integer, [0,n-1])
    variable_order = np.argsort(-symmetric_marginals)

    constraints = np.zeros((m,n))
    remaining_ones = total_ones # the number of ones left to allocate
    variable_index = 0 # allocate variables variable_order[variable_index:variable_index+ones_to_allocate]
    #at each iteration allocate m ones, one to each constraint.
    #except for the last iteration, where we allocate however many ones remain (<= m)
    iteration = 0
    while remaining_ones > 0: 
        iteration += 1
        # print("remaining_ones:", remaining_ones)
        if remaining_ones < m:
            ones_to_allocate = remaining_ones
            remaining_ones = 0
        else:
            ones_to_allocate = m
            remaining_ones -= m
        # variables = variable_order[variable_index:variable_index+ones_to_allocate]

        if variable_index >= n:
            variable_index = variable_index % n
        if variable_index+ones_to_allocate < n:
            variables = variable_order[variable_index:variable_index+ones_to_allocate]
        else:
            ones_at_end = n - variable_index
            ones_at_beginning = ones_to_allocate - ones_at_end
            variables = list(variable_order[variable_index:])
            variables_at_begining = list(variable_order[:ones_at_beginning])
            variables.extend(variables_at_begining)
        assert(len(variables) == ones_to_allocate and len(variable_order) == n)




        variable_index += ones_to_allocate

        assignment_solver = 'lapjv1'        
        cost_matrix = construct_cost_matrix(constraints, variables, SATISFYING_SOLUTIONS, iteration, assignment_solver)
        #solve the assignment problem
        # lin_assign = linear_assignment.LinearAssignment(cost_matrix)
        # solution = lin_assign.solution
        # association_list = zip([i for i in range(len(solution))], solution)
        assert((cost_matrix >= 0).all())
        # print(cost_matrix)
        TEST = False
        if assignment_solver == 'scipy':
            row_ind, col_ind = linear_sum_assignment(cost_matrix)
            association_list = zip(row_ind, col_ind)        
            association_list = [assoc for assoc in association_list]
        elif assignment_solver == 'lapjv':
            association_list = non_square_lapjv(cost_matrix)
            # if TEST:
            #     row_ind, col_ind = linear_sum_assignment(cost_matrix)
            #     association_list_check = zip(row_ind, col_ind)        

            #     a_list = [assoc for assoc in association_list]
            #     a_list1 = [assoc for assoc in association_list_check]
            #     cost = 0
            #     cost1 = 0
            #     for (row, col) in a_list:
            #         cost += cost_matrix[row, col]
            #     for (row, col) in a_list1:
            #         cost1 += cost_matrix[row, col]

            #     for idx, assoc in enumerate(a_list):
            #         assert(assoc == a_list1[idx] or cost == cost1), (a_list, a_list1, cost_matrix, cost, cost1)
            #     print('lapjv checked out :)')   

        elif assignment_solver == 'lapjv1':
            association_list = non_square_lapjv1(cost_matrix)
            # if TEST:
            #     association_list_check = non_square_lapjv(cost_matrix)

            #     a_list = [assoc for assoc in association_list]
            #     a_list1 = [assoc for assoc in association_list_check]
            #     cost = 0
            #     cost1 = 0
            #     for (row, col) in a_list:
            #         cost += cost_matrix[row, col]
            #     for (row, col) in a_list1:
            #         cost1 += cost_matrix[row, col]

            #     assert(np.allclose(cost, cost1)), (cost, cost1)
            #     # print("cost:", cost, "cost_matrix:", cost_matrix)

            #     for idx, assoc in enumerate(a_list):
            #         assert(assoc == a_list1[idx] or np.allclose(cost, cost1)), (a_list, a_list1, cost_matrix, cost, cost1)
            #     print('lapjv1 checked out :)')




        # print("association_list:", association_list)
        # print("cost_matrix:", cost_matrix)
        a_list = [assoc for assoc in association_list]
        # print("a_list:", a_list)
        for (var_idx, constr_idx) in association_list:
            # print(var_idx, constr_idx)
            assert(constraints[constr_idx, variables[var_idx]] == 0), (association_list, constraints[constr_idx, var_idx], constraints, constr_idx, var_idx, cost_matrix)
            constraints[constr_idx, variables[var_idx]] = 1
        if assignment_solver != 'scipy':
            assert((np.sum(constraints, axis=1) == iteration).all() or remaining_ones == 0), np.sum(constraints, axis=1)
        # print("np.sum(constraints, axis=1):", np.sum(constraints, axis=1))
    # return constraints
    if assignment_solver != 'scipy':
        assert(np.sum(constraints) == total_ones)

    a = np.sum(SATISFYING_SOLUTIONS, axis=0)[-55:]/SOLUTION_COUNT
    print("marginals =", np.array2string(a, max_line_width=np.inf))
    print("constraint matrix:")
    b = constraints[:, -55:]
    print(np.array2string(b, max_line_width=np.inf))

    clauses = []
    list_of_marginals = []
    product_of_marginals = 1.0
    for constr_idx in range(m):
        cur_clause = []
        for var_idx in range(n):
            if constraints[constr_idx, var_idx] == 1:
                cur_clause.append(var_idx)
        clauses.append(cur_clause)
        cur_symmetric_marginal, cur_marginal = get_clause_marginal(SATISFYING_SOLUTIONS, constraints[constr_idx])
        list_of_marginals.append(cur_symmetric_marginal)
        product_of_marginals *= cur_symmetric_marginal

    list_of_marginals.sort()
    return clauses, product_of_marginals, list_of_marginals

def construct_cost_matrix(constraints, variables, SATISFYING_SOLUTIONS, iteration, assignment_solver):
    '''
    Inputs:
    - constraints: (np.array of size mxn) the constraint matrix
    - variables: (np.array of length ones_to_allocate) each element is an integer in [0, n-1],
                 representing a variable that needs to be assigned to a constraint
    - SATISFYING_SOLUTIONS: (list of np.arrays of length n), 0/1 valued representing solutions

    Outputs:
    - cost_matrix: (np.array of shape ones_to_allocate x m) cost_matrix[i][j] is the -log(marginal)
                   of the jth constraint if the ith variable is added to it.  Note we do not allow
                   adding the same variable twice.  An assignment in the cost matrix is then
                   the negative logarithm of the product of all marginals.  The minimum cost assignent
                   gives the maximum product of marginals.
    '''
    m = constraints.shape[0]
    n = constraints.shape[1]
    # print("n:", n, "m:", m)
    ones_to_allocate = len(variables)
    # print("HIHIHIIHH from construct_cost_matrix", "iteration =", iteration, "ones_to_allocate =", ones_to_allocate)
    cost_matrix = np.zeros((ones_to_allocate, m))
    for col in range(m):
        clause = constraints[col]
        for row in range(ones_to_allocate):
            variable_to_set = variables[row]
            if clause[variable_to_set] == 1:
                cost_matrix[row][col] = np.inf
            else:
                assert(clause[variable_to_set] == 0)
                # double_marginal = False
                global DOUBLE_MARGINAL
                if DOUBLE_MARGINAL:
                    symmetric_marginal, marginal = get_clause_double_marginal(SATISFYING_SOLUTIONS, clause, variable_to_set)

                else:
                    clause[variable_to_set] = 1
                    symmetric_marginal, marginal = get_clause_marginal(SATISFYING_SOLUTIONS, clause)
                    if iteration > 1 and assignment_solver != 'scipy':
                        assert(np.sum(clause) == iteration), (clause, constraints, m, n, ones_to_allocate, iteration)
                        # print("clause:", clause)
                    assert(symmetric_marginal >= 0 and symmetric_marginal <= .5)
                    clause[variable_to_set] = 0
                # print(cost_matrix)
                # print(row, col)
                if symmetric_marginal == 0:
                    cost_matrix[row][col] = 999999999999999999999999
                else:
                    cost_matrix[row][col] = -np.log(symmetric_marginal)
                # print(cost_matrix)
    assert((cost_matrix != 0).all())
    return cost_matrix

def sample_biregular_guaranteeEvenMarginals(constraint_degseq, variable_degseq, marginal_cutoff, prv_constraint_count=-1, prv_constraint_counts=[]):
    '''

    Inputs:
    - prv_constraint_count: (int) the number of constraints in the previous recursive call.  if we get stuck with the same
        number of constraints, don't require marginals close to 1.
    '''
    # print("sample_biregular_guaranteeEvenMarginals called with", len(constraint_degseq), "constraints")
    # print("sum(constraint_degseq) =", sum(constraint_degseq), "constraint_degseq:", constraint_degseq)
    # print("sum(variable_degseq) =", sum(variable_degseq), "variable_degseq:", variable_degseq)
    # print()
    g = bipartite.configuration_model(constraint_degseq, variable_degseq,
                                      create_using=nx.Graph())

    num_constraints = len(constraint_degseq)
    n = len(variable_degseq)

    clauses = []
    number_of_ones = 0
    for j in range(num_constraints):
        clause_independentSet_indices = [ind_ - num_constraints for ind_ in g.neighbors(j)]
        clauses.append(clause_independentSet_indices)
        number_of_ones += len(clauses[-1])

    product_of_marginals = 1.0
    list_of_marginals = []

    global conf
    HACKY = True
    if HACKY:
        # for idx, clause in enumerate(clauses):
        #     assert(constraint_degseq[idx] == len(clause)), (constraint_degseq[idx], len(clause), number_of_ones, sum(constraint_degseq), sum(variable_degseq))

        clauses_with_even_marginals = []
        for clause_idx, clause in reversed(list(enumerate(clauses))): #reverse list so deleteing doesn't mess up indices
            clause_array = np.zeros(n)
            for v in clause:
                assert(v >= 0 and v <= n-1), (v, n)
                clause_array[v] = 1.0
            cur_symmetric_marginal, cur_marginal = get_clause_marginal(SATISFYING_SOLUTIONS, clause_array)    
            # if cur_symmetric_marginal > marginal_cutoff or prv_constraint_count == num_constraints:
            assert(conf.extra_configs['biregular_marginal_problem_constraint'] in ['keep_biregular', 'iid'])
            #'keep_biregular': let the problem clauses remain with bad marginals and keep biregular
            #'iid': sample the problem clausees iid with good marginals giving up biregular

            if ((conf.extra_configs['biregular_marginal_problem_constraint'] == 'keep_biregular') and\
                (cur_symmetric_marginal > marginal_cutoff or (len(prv_constraint_counts) >= 10 and prv_constraint_counts[-10] == num_constraints)))\
               or\
               ((conf.extra_configs['biregular_marginal_problem_constraint'] == 'iid') and (cur_symmetric_marginal > marginal_cutoff)):
                var_count_to_delete = constraint_degseq[clause_idx]
                del constraint_degseq[clause_idx]
                # assert(constraint_degseq[clause_idx] == len(clause)), (constraint_degseq[clause_idx], len(clause))
                for v in clause:
                    # assert(variable_degseq[v] >= 1)
                    if variable_degseq[v] >= 1:
                        variable_degseq[v] -= 1
                        var_count_to_delete -= 1

                    if var_count_to_delete == 0:
                        break
                # print("var_count_to_delete:", var_count_to_delete)
                while var_count_to_delete > 0:
                    # print("hi")
                    var_probs = np.array(variable_degseq)
                    var_probs = var_probs/np.sum(var_probs)
                    var_to_delete = np.random.choice(len(variable_degseq), p=var_probs)

                    # print("var_to_delete:", var_to_delete)
                    assert(variable_degseq[var_to_delete] >= 1), (variable_degseq, var_probs)
                    variable_degseq[var_to_delete] -= 1
                    var_count_to_delete -= 1
                assert(sum(constraint_degseq) == sum(variable_degseq)), (sum(constraint_degseq), sum(variable_degseq), constraint_degseq, variable_degseq)
                clauses_with_even_marginals.append(clause)
                product_of_marginals *= cur_symmetric_marginal
                list_of_marginals.append(cur_symmetric_marginal)
    else:        
        for idx, clause in enumerate(clauses):
            assert(constraint_degseq[idx] == len(clause)), (constraint_degseq[idx], len(clause), number_of_ones, sum(constraint_degseq), sum(variable_degseq))

        clauses_with_even_marginals = []
        for clause_idx, clause in reversed(list(enumerate(clauses))): #reverse list so deleteing doesn't mess up indices
            clause_array = np.zeros(n)
            for v in clause:
                assert(v >= 0 and v <= n-1), (v, n)
                clause_array[v] = 1.0
            cur_symmetric_marginal, cur_marginal = get_clause_marginal(SATISFYING_SOLUTIONS, clause_array)    
            if cur_symmetric_marginal > marginal_cutoff or prv_constraint_count == num_constraints:
                del constraint_degseq[clause_idx]
                assert(constraint_degseq[clause_idx] == len(clause)), (constraint_degseq[clause_idx], len(clause))
                for v in clause:
                    assert(variable_degseq[v] >= 1)
                    variable_degseq[v] -= 1
                assert(sum(constraint_degseq) == sum(variable_degseq)), (sum(constraint_degseq), sum(variable_degseq), constraint_degseq, variable_degseq)
                clauses_with_even_marginals.append(clause)

    if len(constraint_degseq) >= 1:
        if (len(prv_constraint_counts) >= 10 and prv_constraint_counts[-10] == num_constraints):
            # print("we have to sample", len(constraint_degseq), "constraints iid")
            remaining_clauses, remaining_clauses_product_of_marginals, list_of_remaining_marginals = sample_iid_guaranteeEvenMarginals(\
                constraint_degseq, variable_degseq, marginal_cutoff)
        else:
            prv_constraint_counts.append(num_constraints)
            remaining_clauses, remaining_clauses_product_of_marginals, list_of_remaining_marginals = sample_biregular_guaranteeEvenMarginals(\
                constraint_degseq, variable_degseq, marginal_cutoff, prv_constraint_count=num_constraints, prv_constraint_counts=prv_constraint_counts)

        clauses_with_even_marginals.extend(remaining_clauses)
        product_of_marginals *= remaining_clauses_product_of_marginals
        list_of_marginals.extend(list_of_remaining_marginals)


    list_of_marginals.sort()
    if(prv_constraint_count == num_constraints):
        pass
        # print("fudged these constraints, list of marginals:", list_of_marginals)

    return clauses_with_even_marginals, product_of_marginals, list_of_marginals



def sample_iid_guaranteeEvenMarginals(constraint_degseq, variable_degseq, marginal_cutoff):
    '''

    '''
    if marginal_cutoff < .00001:
        marginal_cutoff = 0
    num_constraints = len(constraint_degseq)
    n = len(variable_degseq)
    assert(sum(constraint_degseq) == sum(variable_degseq))
    f = sum(constraint_degseq)/(n*num_constraints)
    def sample_single_clause_iid(f, n):
        clause = []
        for v in range(n):
            if np.random.rand() < f:
                clause.append(v)
        return clause

    clauses_with_even_marginals = []
    product_of_marginals = 1.0
    list_of_marginals = []

    tries = 0
    TRIES_BEFORE_LOWERING_MARGINAL_CUTOFF = 20
    while num_constraints > 0:
        if tries == TRIES_BEFORE_LOWERING_MARGINAL_CUTOFF:
            marginal_cutoff *= .95
            tries = 0
            print("couldn't find an even constraint in", TRIES_BEFORE_LOWERING_MARGINAL_CUTOFF, "tries, lowered marginal cutoff to", marginal_cutoff)
        iid_clause = sample_single_clause_iid(f, n)
        clause_array = np.zeros(n)
        for v in iid_clause:
            assert(v >= 0 and v <= n-1), (v, n)
            clause_array[v] = 1.0
        cur_symmetric_marginal, cur_marginal = get_clause_marginal(SATISFYING_SOLUTIONS, clause_array)    
        if cur_symmetric_marginal >= marginal_cutoff:
            clauses_with_even_marginals.append(iid_clause)
            num_constraints -= 1
            product_of_marginals *= cur_symmetric_marginal
            list_of_marginals.append(cur_symmetric_marginal)
            tries = 0
        else:
            tries += 1
    return clauses_with_even_marginals, product_of_marginals, list_of_marginals


def generateLDPC_bipartiteGraph_parameterizeF_pickEvenSampleSplit_perConstraint(n_constr, f, method, marginal_cutoff=.49):
    """
    Sample a complete biregular matrix, check the marginal for each constraint.  Keep constraints with
    marginals sufficiently close to .5, but resample a new matrix with the same number of ones in each row
    and specified numbers of ones in each column (so when combined with the kept rows, each column has the
    same number of ones).  Recursively keep rows/resample the rest until all rows have marginals close to .5

    Generate the constraints corresponding to an LDPC code.
    The variable set is taken from the configuration. If extra_vars > 0, a
    number of "dummy" vars participate in the generation of the constraints.

    :param n_constr:    The number of constraints to generate
    :param f:  density of ones

    :return: A list of constraints repr. as lists of strings ending with nl

    Inputs:
    - marginal_cutoff: (float in [0,.5)) marginals for each constraint must be in (marginal_cutoff, 1 - marginal_cutoff)
    """
    max_product_of_marginals = -np.inf
    clauses_with_max_product_of_marginals = None
    i_effective_for_max_product_of_marginals = None
    global SATISFYING_SOLUTIONS

    if VERBOSE:
        print('f=', f, 'var degree=', n_constr*f, 'clause degree=', len(conf.var_set)*f)

    assert n_constr > 0
    n0 = len(conf.var_set)
    var_list = conf.var_set
    n = len(var_list)
    i = int(n_constr)
    degn = i*f # average variable degree (ones per column)
    x = n * f  # Average constraint degree
    # log.info('Generating LDPC: n_constr={}  n_var={} var_degree={}  '
    #          'constr_degree={:.1f}'.format(n_constr, n, degn, x))
    constraint_degseq = [math.floor(x)] * i
    variable_degseq = [math.floor(degn)] * n
    if degn > math.floor(degn):
        x_ = math.ceil((degn - math.floor(degn)) * n)
        indices_ = random.sample(range(n), x_)
        for j_ in indices_:
            variable_degseq[j_] += 1
    surplus = sum(variable_degseq) - sum(constraint_degseq)
    assert 0 <= surplus <= n, 'Improper surplus  = ' + str(surplus)
    for c_ in random.sample(range(i), surplus):
        constraint_degseq[c_] += 1
    # print("constraint_degseq:")
    # print(constraint_degseq)
    # print("variable_degseq:")
    # print(variable_degseq)
    
    constraint_deg_dictionary = defaultdict(int)
    for deg in constraint_degseq:
        constraint_deg_dictionary[deg] += 1


    # print("sum(constraint_degseq):", sum(constraint_degseq), "constraint_degseq:", constraint_degseq)
    # print("sum(variable_degseq):", sum(variable_degseq), "variable_degseq:", variable_degseq)

    if method == 'bi_regular_marginals_per_constraint':
        clauses_with_even_marginals, product_of_marginals, list_of_marginals = sample_biregular_guaranteeEvenMarginals(constraint_degseq, variable_degseq, marginal_cutoff)
    elif method in ['bi_regular_order_vars_by_marginals', 'bi_regular_order_vars_by_double_marginals']:
        assert(sum(constraint_degseq) == sum(variable_degseq))
        clauses_with_even_marginals, product_of_marginals, list_of_marginals = create_biregular_adding_variables_orderByMarginals_assignmentProblem(n, m=n_constr, total_ones=sum(constraint_degseq))
    elif method == 'bi_regular_order_vars_by_marginals_randomChunks':
        assert(sum(constraint_degseq) == sum(variable_degseq))
        # clauses_with_even_marginals, product_of_marginals, list_of_marginals = create_biregular_adding_variables_orderByMarginals_randomInChunks(n, m=n_constr, total_ones=sum(constraint_degseq))
        clauses_with_even_marginals, product_of_marginals, list_of_marginals = create_biregular_adding_variables_orderByMarginals_randomInChunks_faster(n, m=n_constr, total_ones=sum(constraint_degseq))

    else:
        assert(False)
    # print("product_of_marginals:", product_of_marginals)
    # if list_of_marginals[0] < .01: #try again
        # print('trying generateLDPC_bipartiteGraph_parameterizeF_pickEvenSampleSplit_perConstraint again')    
        # print("list_of_marginals:", list_of_marginals)
        # print("list_of_marginals[:10]", list_of_marginals[:10])
        # pass
        # return generateLDPC_bipartiteGraph_parameterizeF_pickEvenSampleSplit_perConstraint(n_constr, f, method, marginal_cutoff)

    # print("product_of_marginals:", product_of_marginals, "marginals:", list_of_marginals)

    ###### Testing sample_biregular_guaranteeEvenMarginals ######
    actual_constraint_degseq = [0] * i
    actual_variable_degseq = [0] * n
    for idx, clause in enumerate(clauses_with_even_marginals):
        actual_constraint_degseq[idx] = len(clause)
        for v in clause:
            actual_variable_degseq[v] += 1

    # print("sum(actual_constraint_degseq):", sum(actual_constraint_degseq), "actual_constraint_degseq:", actual_constraint_degseq)
    # print("sum(actual_variable_degseq):", sum(actual_variable_degseq), "actual_variable_degseq:", actual_variable_degseq)

    # for clause in clauses_with_even_marginals:
    #     constraint_deg_dictionary[len(clause)] -= 1
    #     for v in clause:
    #         variable_degseq[v] -= 1
    # for deg in variable_degseq:
    #     assert(deg == 0)
    # for deg, count in constraint_deg_dictionary.iteritems():
    #     assert(count == 0)
    ###### Done Testing sample_biregular_guaranteeEvenMarginals ######

    formated_clauses = []
    cur_i_effective = 1.0
    for clause in clauses_with_even_marginals:
        clause_varlist = [var_list[v] for v in clause] #take independent set indices and translate to variable names
        # is_negated = (random.randint(0, 1) > 0)
        cur_marginal = .5
        if np.random.rand() < cur_marginal:
            is_negated = False
            cur_i_effective /= cur_marginal
        else:
            is_negated = True
            cur_i_effective /= (1-cur_marginal)
        if is_negated:
            out = ['x -']
        else:
            out = ['x ']
        for v in clause_varlist:
            out.append('{} '.format(v))
        out.append('0\n')
        cl = ''.join(out)
        formated_clauses.append(cl)


    return formated_clauses, product_of_marginals, -1


def generateLDPC_bipartiteGraph_parameterizeF(n_constr, f):
    """
    Generate the constraints corresponding to an LDPC code.
    The variable set is taken from the configuration. If extra_vars > 0, a
    number of "dummy" vars participate in the generation of the constraints.

    :param n_constr:    The number of constraints to generate
    :param f:  density of ones

    :return: A list of constraints repr. as lists of strings ending with nl
    """
    assert n_constr > 0
    n0 = len(conf.var_set)
    var_list = conf.var_set
    n = len(var_list)
    i = int(n_constr)
    degn = i*f # average variable degree (ones per column)
    x = n * f  # Average constraint degree
    log.info('Generating LDPC: n_constr={}  n_var={} var_degree={}  '
             'constr_degree={:.1f}'.format(n_constr, n, degn, x))
    i_degseq = [math.floor(x)] * i
    n_degseq = [math.floor(degn)] * n
    if degn > math.floor(degn):
        x_ = math.ceil((degn - math.floor(degn)) * n)
        indices_ = random.sample(range(n), x_)
        for j_ in indices_:
            n_degseq[j_] += 1
    surplus = sum(n_degseq) - sum(i_degseq)
    assert 0 <= surplus <= n, 'Improper surplus  = ' + str(surplus)
    for c_ in random.sample(range(i), surplus):
        i_degseq[c_] += 1
    # print("i_degseq:")
    # print(i_degseq)
    # print("n_degseq:")
    # print(n_degseq)
    
    g = bipartite.configuration_model(i_degseq, n_degseq,
                                      create_using=nx.Graph())
    


    clauses_l = []
    for j in range(i):
        clause_varlist = [var_list[ind_ - i] for ind_ in g.neighbors(j)]
        is_negated = (random.randint(0, 1) > 0)
        if is_negated:
            out = ['x -']
        else:
            out = ['x ']
        for v in clause_varlist:
            out.append('{} '.format(v))
        out.append('0\n')
        cl = ''.join(out)
        clauses_l.append(cl)

    return clauses_l


# approximate marginals sat-grid-pbl-0015.cnf
# iid, repeats=2000, f=.5, m=20
approx_marginals_sat_grid_pbl_0015 = \
 [0.     ,0.     ,0.2515 ,0.187  ,0.576  ,0.3575 ,0.355  ,0.3445 ,0.3175 ,0.311
 ,0.4795 ,0.4625 ,0.4185 ,0.4245 ,0.3375 ,0.335  ,0.3665 ,0.411  ,0.342  ,0.549
 ,0.8085 ,0.8835 ,0.0715 ,0.063  ,0.0955 ,0.1065 ,0.0875 ,0.0885 ,0.893  ,0.8965
 ,0.894  ,0.9025 ,0.105  ,0.098  ,0.0965 ,0.096  ,0.0225 ,0.019  ,0.975  ,0.134
 ,0.9    ,0.895  ,0.906  ,0.904  ,0.1075 ,0.108  ,0.1055 ,0.0975 ,0.006  ,0.004
 ,0.9725 ,0.129  ,0.978  ,0.114  ,0.8875 ,0.901  ,0.902  ,0.906  ,0.98   ,0.1155
 ,0.025  ,0.021  ,0.0035 ,0.0015 ,0.975  ,0.123  ,0.975  ,0.1275 ,0.9775 ,0.143
 ,0.898  ,0.8955 ,0.9005 ,0.9045 ,0.9775 ,0.128  ,0.9745 ,0.128  ,0.     ,0.
 ,0.9745 ,0.1305 ,0.9755 ,0.119  ,0.9765 ,0.1185 ,0.9755 ,0.126  ,0.893  ,0.894
 ,0.907  ,0.9055 ,0.9755 ,0.1265 ,0.9755 ,0.1245 ,0.     ,1.     ,0.     ,1.
 ,0.9685 ,0.1325 ,0.9775 ,0.1275 ,0.9755 ,0.1275 ,0.9805 ,0.1175 ,0.9    ,0.8955
 ,0.8965 ,0.884  ,0.975  ,0.125  ,0.976  ,0.124  ,0.9745 ,0.1225 ,0.974  ,0.1275
 ,0.9805 ,0.115  ,0.971  ,0.1255 ,0.9785 ,0.1095 ,0.976  ,0.12   ,0.9805 ,0.1155
 ,0.9045 ,0.9045 ,0.8995 ,0.902  ,0.9775 ,0.113  ,0.975  ,0.135  ,0.9805 ,0.121
 ,0.9665 ,0.1235 ,0.978  ,0.1225 ,0.9725 ,0.1295 ,0.9785 ,0.119  ,0.9725 ,0.122
 ,0.979  ,0.1215 ,0.978  ,0.129  ,0.902  ,0.904  ,0.8975 ,0.891  ,0.9795 ,0.128
 ,0.973  ,0.1225 ,0.973  ,0.13   ,0.979  ,0.111  ,0.9745 ,0.125  ,0.9755 ,0.109
 ,0.978  ,0.117  ,0.973  ,0.1285 ,0.9765 ,0.122  ,0.9785 ,0.1205 ,0.9775 ,0.1225
 ,0.9115 ,0.8965 ,0.8945 ,0.9    ,0.976  ,0.1275 ,0.978  ,0.114  ,0.979  ,0.1215
 ,0.9775 ,0.124  ,0.9785 ,0.1175 ,0.978  ,0.129  ,0.972  ,0.126  ,0.975  ,0.1215
 ,0.978  ,0.128  ,0.972  ,0.135  ,0.977  ,0.1295 ,0.9755 ,0.126  ,0.8915 ,0.9005
 ,0.6305 ,0.5795 ,0.973  ,0.07   ,0.9725 ,0.061  ,0.9715 ,0.068  ,0.9625 ,0.076
 ,0.9655 ,0.0655 ,0.97   ,0.0655 ,0.9705 ,0.066  ,0.962  ,0.075  ,0.965  ,0.069
 ,0.9655 ,0.0785 ,0.9645 ,0.08   ,0.967  ,0.0735 ,0.9735 ,0.062  ,0.604  ,0.6145]

# approximate marginals sat-grid-pbl-0015.cnf
# bi_regular, repeats=2000, f=.1, m=20
approx_marginals_sat_grid_pbl_0015_2 = \
[0.000e+00 ,0.000e+00 ,2.450e-02 ,2.555e-01 ,4.550e-01 ,3.900e-01 ,3.915e-01
 ,3.815e-01 ,2.920e-01 ,2.970e-01 ,5.085e-01 ,4.990e-01 ,4.740e-01 ,4.780e-01
 ,3.405e-01 ,3.285e-01 ,3.765e-01 ,3.600e-01 ,5.390e-01 ,5.570e-01 ,4.980e-01
 ,5.160e-01 ,3.510e-01 ,3.335e-01 ,3.030e-01 ,3.000e-01 ,3.350e-01 ,3.110e-01
 ,7.495e-01 ,7.755e-01 ,7.035e-01 ,7.380e-01 ,1.640e-01 ,1.355e-01 ,8.950e-02
 ,6.750e-02 ,7.150e-02 ,4.650e-02 ,9.570e-01 ,5.000e-02 ,9.745e-01 ,9.780e-01
 ,9.785e-01 ,9.750e-01 ,2.050e-02 ,2.100e-02 ,2.250e-02 ,1.600e-02 ,1.050e-02
 ,8.500e-03 ,9.865e-01 ,1.950e-02 ,9.905e-01 ,1.650e-02 ,9.935e-01 ,9.910e-01
 ,9.915e-01 ,9.925e-01 ,9.975e-01 ,8.000e-03 ,2.000e-03 ,4.000e-03 ,5.000e-04
 ,2.000e-03 ,9.910e-01 ,1.500e-02 ,9.910e-01 ,1.300e-02 ,9.975e-01 ,1.000e-02
 ,9.935e-01 ,9.940e-01 ,9.955e-01 ,9.935e-01 ,9.965e-01 ,9.500e-03 ,9.955e-01
 ,1.450e-02 ,0.000e+00 ,0.000e+00 ,9.975e-01 ,8.000e-03 ,9.935e-01 ,9.500e-03
 ,9.960e-01 ,8.500e-03 ,9.975e-01 ,8.500e-03 ,9.935e-01 ,9.945e-01 ,9.920e-01
 ,9.905e-01 ,9.965e-01 ,7.500e-03 ,9.955e-01 ,1.150e-02 ,0.000e+00 ,1.000e+00
 ,0.000e+00 ,1.000e+00 ,9.955e-01 ,9.500e-03 ,9.935e-01 ,1.150e-02 ,9.925e-01
 ,1.300e-02 ,9.945e-01 ,1.100e-02 ,9.950e-01 ,9.960e-01 ,9.945e-01 ,9.925e-01
 ,9.980e-01 ,8.500e-03 ,9.975e-01 ,8.000e-03 ,9.970e-01 ,8.000e-03 ,9.945e-01
 ,1.100e-02 ,9.935e-01 ,1.400e-02 ,9.930e-01 ,1.350e-02 ,9.950e-01 ,1.150e-02
 ,9.925e-01 ,1.450e-02 ,9.975e-01 ,7.500e-03 ,9.930e-01 ,9.935e-01 ,9.945e-01
 ,9.955e-01 ,9.950e-01 ,1.150e-02 ,9.960e-01 ,9.000e-03 ,9.920e-01 ,1.250e-02
 ,9.930e-01 ,9.000e-03 ,9.955e-01 ,9.500e-03 ,9.940e-01 ,1.350e-02 ,9.955e-01
 ,9.000e-03 ,9.935e-01 ,1.000e-02 ,9.935e-01 ,1.250e-02 ,9.960e-01 ,8.000e-03
 ,9.940e-01 ,9.940e-01 ,9.935e-01 ,9.940e-01 ,9.965e-01 ,8.500e-03 ,9.935e-01
 ,1.250e-02 ,9.955e-01 ,8.500e-03 ,9.925e-01 ,1.250e-02 ,9.930e-01 ,1.250e-02
 ,9.925e-01 ,1.200e-02 ,9.945e-01 ,1.200e-02 ,9.960e-01 ,1.100e-02 ,9.945e-01
 ,9.500e-03 ,9.950e-01 ,8.000e-03 ,9.940e-01 ,1.200e-02 ,9.945e-01 ,9.910e-01
 ,9.920e-01 ,9.915e-01 ,9.970e-01 ,1.200e-02 ,9.930e-01 ,1.400e-02 ,9.940e-01
 ,1.400e-02 ,9.935e-01 ,1.200e-02 ,9.925e-01 ,1.350e-02 ,9.925e-01 ,1.300e-02
 ,9.960e-01 ,8.500e-03 ,9.915e-01 ,1.450e-02 ,9.940e-01 ,1.400e-02 ,9.920e-01
 ,1.300e-02 ,9.960e-01 ,1.100e-02 ,9.960e-01 ,9.500e-03 ,9.945e-01 ,9.915e-01
 ,9.605e-01 ,9.620e-01 ,9.955e-01 ,9.000e-03 ,9.935e-01 ,1.000e-02 ,9.945e-01
 ,8.500e-03 ,9.905e-01 ,1.100e-02 ,9.920e-01 ,1.150e-02 ,9.920e-01 ,1.250e-02
 ,9.930e-01 ,1.050e-02 ,9.910e-01 ,1.150e-02 ,9.925e-01 ,8.500e-03 ,9.935e-01
 ,1.100e-02 ,9.920e-01 ,1.050e-02 ,9.920e-01 ,1.300e-02 ,9.950e-01 ,7.500e-03
 ,9.620e-01 ,9.575e-01]

def generateLDPC_bipartiteGraph_parameterizeF_proportionalMarginals(n_constr, f, approx_marginals):
    """
    Generate the constraints corresponding to an LDPC code.
    The variable set is taken from the configuration. If extra_vars > 0, a
    number of "dummy" vars participate in the generation of the constraints.

    :param n_constr:    The number of constraints to generate
    :param f:  density of ones

    :return: A list of constraints repr. as lists of strings ending with nl
    """
    # approx_marginals = [1/len(approx_marginals) for i in approx_marginals]
    def get_ones_per_variable(marginals, f, n, m):
        unnormalized_probs = []
        for marginal in marginals:
            unnormalized_probs.append(.5 - abs(.5 - marginal))
            # unnormalized_probs.append(abs(.5 - marginal))
        unnormalized_probs = np.array(unnormalized_probs)
        normalized_probs = unnormalized_probs/sum(unnormalized_probs)
        ones_per_variable = [int(math.floor(prob*f*n*m)) for prob in normalized_probs]
        assert(sum(ones_per_variable) < f*n*m), (sum(ones_per_variable), f*n)
        return ones_per_variable, normalized_probs 

    def get_ones_per_variable1(marginals, f, n, m):
        unnormalized_probs = []
        for marginal in marginals:
            if marginal > .05 and marginal < .95:
                unnormalized_probs.append(1.0)
            else:
                unnormalized_probs.append(0.0)
        unnormalized_probs = np.array(unnormalized_probs)
        normalized_probs = unnormalized_probs/sum(unnormalized_probs)
        ones_per_variable = [int(math.floor(prob*f*n*m)) for prob in normalized_probs]
        assert(sum(ones_per_variable) < f*n*m), (sum(ones_per_variable), f*n)
        return ones_per_variable, normalized_probs 


    assert n_constr > 0
    n0 = len(conf.var_set)
    var_list = conf.var_set
    n = len(var_list)
    assert(n == len(approx_marginals))
    i = int(n_constr)
    degn = i*f # average variable degree (ones per column)
    x = n * f  # Average constraint degree
    log.info('Generating LDPC: n_constr={}  n_var={} var_degree={}  '
             'constr_degree={:.1f}'.format(n_constr, n, degn, x))
    i_degseq = [math.floor(x)] * i
    # n_degseq = [math.floor(degn)] * n
    n_degseq, normalized_probs = get_ones_per_variable(approx_marginals, f, n, i)
    print("n_degseq1:", n_degseq)
    if np.sum(n_degseq) < f*n*i:
        # x_ = math.ceil((degn - math.floor(degn)) * n)
        # indices_ = random.sample(normalized_probs, x_)
        # indices_ = random.sample(range(n), x_)
        x_ = math.ceil(f*n*i - np.sum(n_degseq))
        indices_ = np.random.choice(range(n), size=x_, replace=False, p=normalized_probs)
        for j_ in indices_:
            n_degseq[j_] += 1
    surplus = sum(n_degseq) - sum(i_degseq)
    # print("f*n*i:", f*n*i)
    # print("sum(n_degseq):", sum(n_degseq))
    # print("sum(i_degseq):", sum(i_degseq))
    # print("surplus:", surplus)
    # print("i:", i)
    assert 0 <= surplus <= n, 'Improper surplus  = ' + str(surplus)
    for c_ in random.sample(range(i), surplus):
        i_degseq[c_] += 1
    print("i_degseq:")
    print(i_degseq)
    print("n_degseq:")
    print(n_degseq)
    
    g = bipartite.configuration_model(i_degseq, n_degseq,
                                      create_using=nx.Graph())
    clauses_l = []
    for j in range(i):
        clause_varlist = [var_list[ind_ - i] for ind_ in g.neighbors(j)]
        is_negated = (random.randint(0, 1) > 0)
        if is_negated:
            out = ['x -']
        else:
            out = ['x ']
        for v in clause_varlist:
            out.append('{} '.format(v))
        out.append('0\n')
        cl = ''.join(out)
        clauses_l.append(cl)
    return clauses_l


def add_regular_constraints_subsetOfVariables(n, m, f, f_block, k=None, variable_subset=[], clause_sample_repeats=1):
    """ 
    Add m parity constraints, according to the new combined ensemble without decreasing f,
    adding block 1's with probability f_block, 
    and permuting columns

    Inputs:
    - clause_sample_repeats: (int) sample each clause this many times, then pick the one with 
        the most even marginals (or might experiment with other criteria)

    """
    if variable_subset == []:
        variable_subset=[i+1 for i in range(n)]

    m_effective = m

    if k==None or k*m_effective > len(variable_subset): #use k = n/m_effective
        k_low = int(math.floor(float(len(variable_subset)) / m_effective))
        k_high = int(math.ceil(float(len(variable_subset)) / m_effective))
    else:
        k_low = k
        k_high = k
        
    number_k_high_blocks = len(variable_subset)%m_effective
    k_range = [0]
    for i in range(number_k_high_blocks):
        k_range.append(k_range[i] + k_high)
    for i in range(number_k_high_blocks, m_effective):
        k_range.append(k_range[i] + k_low)            
    if k==None or k*m_effective > len(variable_subset): #use k = n/m_effective
        assert(k_range[-1] == len(variable_subset))
#        print k_range
    block_diag_matrix = np.zeros((m, len(variable_subset)))
    #construct block diagonal 1's matrix
    for i in range(0, m):
        for atom in range(1, len(variable_subset)+1):

            if (atom > k_range[i] and atom <= k_range[i+1]):
                block_diag_matrix[i, atom-1] = 1
#        print block_diag_matrix
    #permute the columns of the block diagonal matrix
    permuted_block_diag_matrix = np.swapaxes(np.random.permutation(np.swapaxes(block_diag_matrix,0,1)),0,1)

#        print permuted_block_diag_matrix
    f_updated = f 

    hash_functions = []
    fail_apriori = False

    curIndex = len(variable_subset) + 1

    total_vars_in_parity_constraints = 0

    global SATISFYING_SOLUTIONS

    product_of_marginals = 1.0
    # print('422find me marginals:', end=' ')
    for i in range(0, m):
        new_functions_to_pick_from = []

        for sample_idx in range(clause_sample_repeats):
            new_function = []

            for (var_idx, var_name) in enumerate(variable_subset):
                #check block diagonal construction
                if (var_idx+1 > k_range[i] and var_idx+1 <= k_range[i+1]):
                    assert (block_diag_matrix[i, var_idx] == 1)
                else:
                    assert (block_diag_matrix[i, var_idx] == 0)

                #is this element part of a permuted block?
                if (permuted_block_diag_matrix[i, var_idx] == 1):
                    #if so, add variable with probabilty f_block
                    if random.random() < f_block:
                        new_function.append(var_idx)
                        total_vars_in_parity_constraints += 1
                #if this element isn't part of a permuted block, add variable with probability f_updated
                elif random.random() < f_updated:
                    new_function.append(var_idx)
                    total_vars_in_parity_constraints += 1

            # if len(new_function) == 0:
                # if random.randint(0, 1) == 0:
                #     continue
                # else:
                #     fail_apriori = True
                #     return hash_functions, fail_apriori, 0
            # if (len(new_function) > 0)  and (random.randint(0, 1) == 0):
            #     new_function[0] = -new_function[0]

            new_functions_to_pick_from.append(new_function)


        closest_marginal_to_point5 = -np.inf
        constraint_with_best_marginal = None
        for new_function in new_functions_to_pick_from:
            if len(new_function) == 0:
                cur_symmetric_marginal = 0
            else:
                clause_array = np.zeros(n)
                for v in new_function:
                    # assert(abs(v) >= 1 and abs(v) <= n), (v, n)
                    assert(v >= 0 and v < len(variable_subset)), (v, n)
                    clause_array[v] = 1.0
                cur_symmetric_marginal, cur_marginal = get_clause_marginal(SATISFYING_SOLUTIONS, clause_array)
            if cur_symmetric_marginal > closest_marginal_to_point5:
                closest_marginal_to_point5 = cur_symmetric_marginal
                constraint_with_best_marginal = [variable_subset[v] for v in new_function]
                if (len(constraint_with_best_marginal) > 0)  and (random.randint(0, 1) == 0):
                    constraint_with_best_marginal[0] = -constraint_with_best_marginal[0]

        if len(constraint_with_best_marginal) == 0:
            if random.randint(0, 1) == 0:
                product_of_marginals *= 0                
                continue
            else:
                fail_apriori = True
                return hash_functions, fail_apriori, 0

        # print(closest_marginal_to_point5, end=' ')
        assert(constraint_with_best_marginal is not None)
        product_of_marginals *= closest_marginal_to_point5
        hash_functions.append(constraint_with_best_marginal)

    # print()
    return hash_functions, fail_apriori, product_of_marginals

def add_bi_regular_constraints_parameterizeByF_corrected(n, f, m):
    """ 
    https://arxiv.org/pdf/1707.09467.pdf

    Add m parity constraints over n variables (parity constraint matrix has shape m x n), where:
        - we have a total of floor(f*m*n) non-zero elements in the matrix
        - every column (variable) has floor(f*m) or ceiling(f*m) non-zero elements
        - every row (parity constraint) has floor(f*n) or ceiling(f*n) non-zero elements
    """
    print('add_bi_regular_constraints_parameterizeByF() called, f=', f, 'm=', m, 'n=', n)

    def permute_rows_and_cols(array, iters):
        #permute the rows and columns of array iters times
        for i in range(iters):
            np.random.shuffle(array)
            array = np.transpose(array)
            np.random.shuffle(array)
            array = np.transpose(array)
        return array        

    def generate_bi_regular_constraints2(total_ones, n, m):
        '''
        Inputs:
        - total_ones: the number of 1's in the entire matrix
        - n: number of variables
        - m: number of constraints

        Outputs:
        - constraints: (list of list of ints) contstraints[i][j] is the jth variable (1 to n) in the ith constraint
        '''
        parity_matrix = np.zeros((m, n))

        
        l_low = int(np.floor(f*m)) #every column has at least this many ones
        l_high = l_low + 1 #the rest of the columns have this many ones
        l_high_count = total_ones - l_low*n # the number of columns with l_high ones
        assert(l_high_count >= 0 and l_high_count <= n)

        if l_high_count > 0:
            cur_l = l_high
        else:
            cur_l = l_low
        col_idx = 0
        l_remaining = cur_l
        for one_idx in range(total_ones):
            row_idx = one_idx % m
            parity_matrix[row_idx][col_idx] = 1.0
            l_remaining -= 1
            if l_remaining == 0:
                col_idx += 1
                if col_idx == l_high_count:
                    cur_l = l_low
                l_remaining = cur_l



        k_low = int(np.floor(total_ones/m)) #number of variables in shorter parity constraints
        k_high = int(k_low + 1) #number of variables in longer parity constraints    
        assert(np.sum(parity_matrix) == total_ones), (np.sum(parity_matrix), total_ones, parity_matrix)                   
        assert(np.logical_or(np.sum(parity_matrix, axis=0) == l_low, np.sum(parity_matrix, axis=0) == l_high).all()), (n, m, np.sum(parity_matrix), total_ones, l_low, l_high, k_low, k_high, np.sum(parity_matrix, axis=0), parity_matrix)
        assert(np.logical_or(np.sum(parity_matrix, axis=1) == k_low, np.sum(parity_matrix, axis=1) == k_high).all()), (n, m, np.sum(parity_matrix), total_ones, l_low, l_high, k_low, k_high, np.sum(parity_matrix, axis=1), parity_matrix)


        parity_matrix = permute_rows_and_cols(parity_matrix, 1)

        # go through each one in the matrix number_of_swaps times and swap column index with a randomly matched one
        number_of_swaps = 10
        for i in range(number_of_swaps):
            list_of_locations_of_ones = []
            for row in range(m):
                for col in range(n):
                    if parity_matrix[row][col] == 1:
                        list_of_locations_of_ones.append((row, col))

            np.random.shuffle(list_of_locations_of_ones)
            for one_index in range(0, len(list_of_locations_of_ones), 2):
                if (one_index < len(list_of_locations_of_ones) - 1):
                    (prev_row_a, prev_col_a) = list_of_locations_of_ones[one_index]
                    (prev_row_b, prev_col_b) = list_of_locations_of_ones[one_index + 1]
                    # list_of_locations_of_ones[one_index] = (prev_row_a, prev_col_b)
                    # list_of_locations_of_ones[one_index + 1] = (prev_row_b, prev_col_a)
                    if parity_matrix[(prev_row_b, prev_col_a)] == 0 and parity_matrix[(prev_row_a, prev_col_b)] == 0:
                        parity_matrix[(prev_row_a, prev_col_a)], parity_matrix[(prev_row_b, prev_col_a)] = parity_matrix[(prev_row_b, prev_col_a)], parity_matrix[(prev_row_a, prev_col_a)]
                        parity_matrix[(prev_row_a, prev_col_b)], parity_matrix[(prev_row_b, prev_col_b)] = parity_matrix[(prev_row_b, prev_col_b)], parity_matrix[(prev_row_a, prev_col_b)]
        
        # parity_matrix = np.zeros((m, n))
        # for (row, col) in list_of_locations_of_ones:
        #     parity_matrix[row][col] = 1.0


        assert(np.logical_or(np.sum(parity_matrix, axis=0) == l_low, np.sum(parity_matrix, axis=0) == l_high).all()), (parity_matrix, np.sum(parity_matrix, axis=0))
        assert(np.logical_or(np.sum(parity_matrix, axis=1) == k_low, np.sum(parity_matrix, axis=1) == k_high).all()), (parity_matrix, np.sum(parity_matrix, axis=1))
        assert(np.sum(parity_matrix) == total_ones), (parity_matrix, np.sum(parity_matrix), total_ones)

        constraints = []
        for row in range(m):
            cur_constraint = []
            for col in range(n):
                if parity_matrix[row][col] == 1.0:
                    cur_constraint.append(col + 1)
            constraints.append(cur_constraint)

        return constraints


    total_ones = int(np.round(f*n*m))
    constraints = generate_bi_regular_constraints2(total_ones, n, m)
    total_vars_in_parity_constraints = 0
    fail_apriori = False
    hash_functions = []
    for i in range(0, m):
        new_function = constraints[i]
        total_vars_in_parity_constraints += len(new_function)
        if len(new_function) == 0:
            if random.randint(0, 1) == 0:
                continue
            else:
                fail_apriori = True
                return hash_functions, fail_apriori
        if random.randint(0, 1) == 0:
            new_function[0] = -new_function[0]
        hash_functions.append(new_function)

    print('empirical density = ', total_vars_in_parity_constraints/(n*m))

    return hash_functions, fail_apriori

def generateLDPC_ours(n_constr, extra_vars, degn, clause_sample_repeats):
    """
    Generate the constraints corresponding to an LDPC code.
    The variable set is taken from the configuration. If extra_vars > 0, a
    number of "dummy" vars participate in the generation of the constraints.

    :param n_constr:    The number of constraints to generate
    :param extra_vars:  Number of extra variables
    :param degn:        Average variable (left) degree

    :return: A list of constraints repr. as lists of strings ending with nl
    """

    print(degn/n_constr, .5)
    f = min(degn/n_constr, .5)
    f = .015

    variable_MIS = conf.var_set
    N = len(variable_MIS)
    m = n_constr

    if m == 1:
        cur_k = 0
    else:
        cur_k = np.floor(N/m)                            
    # cur_k = fw_spec['k']

    k_density = cur_k/N
    # print 'N=', N, 'm=', m, "fw_spec['k']=", fw_spec['k'], 'k_density=', k_density

    #compute the density of ones from k:    
    print("k_density:", k_density, 'n', N, 'm', m, "f", f)
    if k_density > f: #we need to decrease k
        cur_k = np.floor(f*N)
        print("changed cur_k=", cur_k)

    if k_density == f:
        f_prime = 0
    else:
        f_prime = (f*N - cur_k)/(N - 2*cur_k)
    print('f_prime=', f_prime)
    assert(abs((1 - f_prime)*cur_k + f_prime*(N - cur_k) - N*f) < .0001), (f_prime, cur_k, N, (1 - f_prime)*cur_k + f_prime*(N - cur_k) - N*f)

    f_block_val = 1.0-f_prime

    # sat.add_regular_constraints_subsetOfVariables(m=m, f=f_prime, f_block=f_block_val, permute=fw_spec['permute'], k=cur_k,\
    #     change_var_names=fw_spec['change_var_names'], variable_subset=variable_MIS)

    hash_functions, fail_apriori, product_of_marginals = add_regular_constraints_subsetOfVariables(n=N, m=m, f=f_prime, f_block=f_block_val, k=cur_k, variable_subset=variable_MIS, clause_sample_repeats=clause_sample_repeats)
   
    #testing implemenatation of biregular constraints
    # hash_functions, fail_apriori = add_bi_regular_constraints_parameterizeByF_corrected(n=N, f=f, m=m)
    
    print("n:", N, "m:", m, "k", cur_k, "f:", f, "f_prime:", f_prime, "hash_functions:")
    # print(hash_functions)
    # assert(fail_apriori == False)

    clauses_l = []
    for hash_function in hash_functions:
        is_negated = (random.randint(0, 1) > 0)
        # if is_negated:
        #     out = ['x -']
        # else:
        #     out = ['x ']
        out = ['x ']
        for v in hash_function:
            out.append('{} '.format(v))
        out.append('0\n')
        cl = ''.join(out)
        clauses_l.append(cl)
    return clauses_l, fail_apriori, product_of_marginals

def generate_iid(n_constr, extra_vars, degn, clause_sample_repeats=1, f=None):
    """
    Generate the constraints corresponding to an LDPC code.
    The variable set is taken from the configuration. If extra_vars > 0, a
    number of "dummy" vars participate in the generation of the constraints.

    :param n_constr:    The number of constraints to generate
    :param extra_vars:  Number of extra variables
    :param degn:        Average variable (left) degree

    :return: A list of constraints repr. as lists of strings ending with nl
    """

    if f is None:
        f = min(degn/n_constr, .5)

    variable_MIS = conf.var_set
    N = len(variable_MIS)
    m = n_constr

    cur_k = 0                        
    # cur_k = fw_spec['k']


    f_prime = f
    # print('f_prime=', f_prime)
    assert(abs((1 - f_prime)*cur_k + f_prime*(N - cur_k) - N*f) < .0001), (f_prime, cur_k, N)

    f_block_val = 1.0-f_prime

    # sat.add_regular_constraints_subsetOfVariables(m=m, f=f_prime, f_block=f_block_val, permute=fw_spec['permute'], k=cur_k,\
    #     change_var_names=fw_spec['change_var_names'], variable_subset=variable_MIS)

    hash_functions, fail_apriori, product_of_marginals = add_regular_constraints_subsetOfVariables(n=N, m=m, f=f_prime, f_block=f_block_val, k=cur_k, variable_subset=variable_MIS, clause_sample_repeats=clause_sample_repeats)
    print("product_of_marginals:", product_of_marginals)
    # print("n:", N, "m:", m, "k", cur_k, "f:", f, "f_prime:", f_prime, "hash_functions:")
    # print(hash_functions)
    # assert(fail_apriori == False)

    # print("fail_apriori:", fail_apriori)
    # print("hash_functions:", hash_functions)
    if fail_apriori:
        return [], fail_apriori, 0

    clauses_l = []
    for hash_function in hash_functions:
        is_negated = (random.randint(0, 1) > 0)
        # if is_negated:
        #     out = ['x -']
        # else:
        #     out = ['x ']
        out = ['x ']
        for v in hash_function:
            out.append('{} '.format(v))
        out.append('0\n')
        cl = ''.join(out)
        clauses_l.append(cl)
    return clauses_l, fail_apriori, product_of_marginals

def setup_augmented_formula(n_constr, constraints):
    """
    Create a DIMACS file with the original formula augmented with LDPC XOR
    constraints.
    :param n_constr: The number of XORs to be actually added to the original
                     formula. It may be less or equal to the number of
                     constraints supplied with 'constraints' parameter.
    :param constraints: A list of constraints (repr. as strings ending with
                        newline)
    :return: The full pathname of the augmented formula.
    """
    global aug_counter, formula_lines
    if n_constr == 0:
        return conf.formula_filename
    # check_time()
    # check_time_jdk()
    aug_counter += 1
    augmented_filename = 'augform_{}.cnf'.format(aug_counter)
    outputfile = os.path.join(conf.working_dir, augmented_filename)
    print("augmented outputfile:", outputfile)
    # print("3 conf.working_dir:", conf.working_dir)    
    with open(outputfile, 'w') as ofile:
        for iline in formula_lines:
            if str(iline.strip()[0:5]).lower() == 'p cnf':
                fields = iline.strip().split(' ')
                num_variables = int(fields[2])
                num_clauses = int(fields[3]) + n_constr
                ofile.write('p cnf {} {}\n'.format(num_variables, num_clauses))
                continue
            ofile.write(iline)
        print("n_constr:", n_constr)
        for l in constraints[:n_constr]:
            ofile.write(l)
    return outputfile

def setup_formula_xorOnly(constraints):
    """
    setup a formula of only xor constraints
    """
    global aug_counter, formula_lines
    # check_time()
    # check_time_jdk()
    aug_counter += 1
    augmented_filename = 'augform_{}.cnf'.format(aug_counter)
    outputfile = os.path.join(conf.working_dir, augmented_filename)
    print("xor only outputfile:", outputfile)
    # print("3 conf.working_dir:", conf.working_dir)    
    with open(outputfile, 'w') as ofile:
        for iline in formula_lines:
            if str(iline.strip()[0:5]).lower() == 'p cnf':
                fields = iline.strip().split(' ')
                num_variables = int(fields[2])
                num_clauses = len(constraints)
                ofile.write('p cnf {} {}\n'.format(num_variables, num_clauses))
                continue
        for l in constraints:
            ofile.write(l)
    return outputfile

def get_non_zero_constraints(constraints, max_time=100):
    '''
    return only those constraint matrices that do not fail a priori 
    '''
    good_constraints = []
    for cur_constraints in constraints:

        aug_form = setup_formula_xorOnly(cur_constraints)
        time_out, b, n, ctime = sat_count(aug_form, max_time, 1)
        if n > 0:
            good_constraints.append(cur_constraints)
    return good_constraints


def nCr(n, r):
    #https://stackoverflow.com/questions/4941753/is-there-a-math-ncr-function-in-python
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer / denom

def setup_duplicated_augmented_formula(constraints):
    """
    Create a DIMACS file with the original formula augmented with LDPC XOR
    constraints.
    :param n_constr: The number of XORs to be actually added to the original
                     formula. It may be less or equal to the number of
                     constraints supplied with 'constraints' parameter.
    :param constraints: constraints[i] is a list of constraints for the ith duplication
                       (repr. as strings ending with newline)
    :return: The full pathname of the augmented formula.
    """
    global aug_counter, formula_lines
    n_constr = 0
    for cur_constraints in constraints:
        n_constr += len(cur_constraints)

    aug_counter += 1
    augmented_filename = 'augform_{}.cnf'.format(aug_counter)
    outputfile = os.path.join(conf.working_dir, augmented_filename)
    print("duplicated outputfile:", outputfile)
    duplication_factor = len(constraints)

    num_variables = None
    with open(outputfile, 'w') as ofile:
        for iline in formula_lines:
            if str(iline.strip()[0:5]).lower() == 'p cnf':
                fields = iline.strip().split(' ')
                num_variables = int(fields[2])
                num_variables_with_duplication = num_variables*duplication_factor + duplication_factor
                # num_clauses = int(fields[3]) + n_constr
                original_num_clauses = int(fields[3])
                num_clauses = original_num_clauses*duplication_factor + 1 + nCr(duplication_factor, 2) + n_constr + num_variables*duplication_factor
                # num_clauses = original_num_clauses*duplication_factor + n_constr
                ofile.write('p cnf {} {}\n'.format(num_variables_with_duplication, num_clauses))
                continue
            if iline.split(' ')[:1] == ['c', 'ind']:
                assert(num_variables is not None)                
                for dup_idx in range(duplication_factor):
                    converted_independent_set_clause = convert_independentSet_with_duplicates(clause=iline, n_vars=num_variables, dup_idx=dup_idx)
                    ofile.write(converted_independent_set_clause)
            elif iline.split(' ')[0] == 'c':
                continue
            else:
                assert(num_variables is not None)
                if iline == '\n':
                    continue                
                # double check that variables are named 1 through num_variables with no gaps
                # print("iline:", iline)
                for var in iline.split(' ')[:-1]:
                    # print(var)
                    var_as_int = abs(int(var))
                    assert(var_as_int >= 0 and var_as_int <= num_variables), (var, var_as_int)
                # done double check that variables are named 1 through num_variables with no gaps

                for dup_idx in range(duplication_factor):
                    duplicate_converted_clause = convert_standardClause_with_duplicates(clause=iline, n_vars=num_variables, dup_idx=dup_idx, duplication_factor=duplication_factor)
                    ofile.write(duplicate_converted_clause)
        assert(num_variables is not None) 

        #add the (n choose 2) + n + 1 clauses for picking a single problem
        #we have 2 new variables for each duplicate problem so we can control parity constraints
        problem_vars_a = [var for var in range(num_variables*duplication_factor+1, num_variables*duplication_factor + duplication_factor + 1)]
        problem_vars_b = [var for var in range(num_variables*duplication_factor + duplication_factor + 1, num_variables*duplication_factor + 2*duplication_factor + 1)]
        #at least one of the duplicate problems must be true
        ofile.write('c at least one extra variable must be true\n')   
        clause_at_least_one = []
        for v in problem_vars_a:
            clause_at_least_one.append('{} '.format(v))
        clause_at_least_one.append('0\n')
        clause_at_least_one = ''.join(clause_at_least_one)  
        ofile.write(clause_at_least_one)   
        #at most one of the duplicate problems must be true, use variables b
        ofile.write('c at most one extra variable must be true\n')           
        # comb = combinations(problem_vars_b, 2)
        comb = combinations(problem_vars_a, 2)
        for v1, v2 in list(comb): 
            cur_clause = []        
            cur_clause.append('{} '.format(-v1))
            cur_clause.append('{} '.format(-v2))
            cur_clause.append('0\n')
            cur_clause = ''.join(cur_clause)  
            ofile.write(cur_clause)   
        # #link variables a and b
        # for var_a_idx, var_a in enumerate(problem_vars_a):
        #     var_b = problem_vars_b[var_a_idx]
        #     cur_clause = []        
        #     cur_clause.append('{} '.format(-var_a))
        #     cur_clause.append('{} '.format(var_b))
        #     cur_clause.append('0\n')
        #     cur_clause = ''.join(cur_clause)  
        #     ofile.write(cur_clause)   



        # force (hopefully) a different variable to begin each clause per constraint matrix
        # negated_vars_for_each_constraint_matrix = []
        # reordered_constraints = []
        # for cur_constraints in constraints:
        #     reordered_cur_constraints, negated_vars = re_order_vars_within_constraints(constraints=cur_constraints, n_vars=num_variables, num_constraints=len(cur_constraints))
        #     negated_vars_for_each_constraint_matrix.append(negated_vars)
        #     reordered_constraints.append(reordered_cur_constraints)

        #add the num_variables*duplication_factor clauses for making sure that the non-selected problems have only one solution
        for problem_idx, problem_var in enumerate(problem_vars_a):
            for v in range(1, num_variables+1):
                dup_v = problem_idx*num_variables + v
                cur_clause = []        
                cur_clause.append('{} '.format(problem_var))
                # if v in negated_vars_for_each_constraint_matrix[problem_idx]:
                cur_clause.append('{} '.format(dup_v)) 
                cur_clause.append('0\n')
                cur_clause = ''.join(cur_clause)  
                ofile.write(cur_clause)                   


        #add xor constraints        
        for dup_idx, cur_constraints in enumerate(constraints):
            for l in cur_constraints:
                if l == '\n':
                    continue
                duplicate_converted_constraint = convert_xor_clause_with_duplicates(clause=l, n_vars=num_variables, dup_idx=dup_idx, problem_var_a=problem_vars_a[dup_idx], problem_var_b=problem_vars_b[dup_idx])
                ofile.write(duplicate_converted_constraint)

            
    return outputfile

def re_order_vars_within_constraints(constraints, n_vars, num_constraints):
    '''
    Reorder the constraints so that each variable has a different constraint first.

    Fail if this isn't possible
    '''
    cost_matrix = np.zeros((num_constraints, n_vars))
    assert(num_constraints == len(constraints))
    for constraint_idx, constraint in enumerate(constraints):
        original_clause_as_list = constraint.split(' ')
        assert(original_clause_as_list[0] == 'x')
        assert(original_clause_as_list[-1] == '0\n'), original_clause_as_list
        for char_v in original_clause_as_list[1:-1]:        
            v = abs(int(char_v))
            cost_matrix[constraint_idx][v-1] = -1
            association_list = non_square_lapjv1(cost_matrix)

    cost = 0
    for (row, col) in association_list:
        cost += cost_matrix[row][col]
    assert(cost == -num_constraints), ("whoops, looks like we can't have a different variable beginning each clause")

    negated_vars = []
    output_constraints = []
    for constraint_idx, constraint in enumerate(constraints):
        original_clause_as_list = constraint.split(' ')
        assert(original_clause_as_list[0] == 'x')
        assert(original_clause_as_list[-1] == '0\n'), original_clause_as_list
        assert(association_list[constraint_idx][0] == constraint_idx)
        first_var = association_list[constraint_idx][1] + 1
        assert(str(first_var) in original_clause_as_list)

        if np.sign(int(original_clause_as_list[1])) < 0:
            is_negated = True
        else:
            is_negated = False
        if is_negated:
            output_clause_as_list = ['x -']
            negated_vars.append(first_var)
        else:
            output_clause_as_list = ['x ']

        output_clause_as_list.append('{} '.format(first_var))    
        for char_v in original_clause_as_list[1:-1]:
            v = abs(int(char_v))
            if v == first_var:
                continue
            output_clause_as_list.append('{} '.format(v))
        output_clause_as_list.append('0\n')
        assert(len(output_clause_as_list) == len(original_clause_as_list))
        output_clause = ''.join(output_clause_as_list)
        output_constraints.append(output_clause)
        print("output_constraints:")
        print(output_constraints)

    return output_constraints, negated_vars


def convert_xor_clause_with_duplicates(clause, n_vars, dup_idx, problem_var_a, problem_var_b):
    '''
    Inputs:
    - clause: (str) a string representation of a clause
    - n_var: (int) the number of variables in the original formula
    - dup_idx: (int) the duplicate index of this set of variables.  e.g.
        if this is the ith duplicate, variables [0, n_var-1] will be translated
        to [dup_idx*n_var, (dup_idx+1)*n_var - 1] 

    Outputs:
    - output_clause: (str) a string representation of the clause with variable names increased appropriately
    '''
    original_clause_as_list = clause.split(' ')
    assert(original_clause_as_list[0] == 'x')
    assert(original_clause_as_list[-1] == '0\n'), original_clause_as_list
    clause_length = len(original_clause_as_list) - 2
    # print(original_clause_as_list[1])
    if np.sign(int(original_clause_as_list[1])) < 0:
        is_negated = True
        clause_length_for_determining_problem_var_sign = clause_length - 1
    else:
        is_negated = False
        clause_length_for_determining_problem_var_sign = clause_length

    if is_negated:
        output_clause_as_list = ['x -']
    else:
        output_clause_as_list = ['x ']
    for char_v in original_clause_as_list[1:-1]:
        v = abs(int(char_v))
        v = v + dup_idx*n_vars
        output_clause_as_list.append('{} '.format(v))


    if clause_length_for_determining_problem_var_sign % 2 == 0:
        output_clause_as_list.append('{} '.format(-problem_var_a))
    else:
        pass
        # output_clause_as_list.append('{} '.format(problem_var))

    output_clause_as_list.append('0\n')
    output_clause = ''.join(output_clause_as_list)
    return output_clause

def convert_independentSet_with_duplicates(clause, n_vars, dup_idx):
    '''
    Inputs:
    - clause: (str) a string representation of a clause
    - n_var: (int) the number of variables in the original formula
    - dup_idx: (int) the duplicate index of this set of variables.  e.g.
        if this is the ith duplicate, variables [0, n_var-1] will be translated
        to [dup_idx*n_var, (dup_idx+1)*n_var - 1] 

    Outputs:
    - output_clause: (str) a string representation of the clause with variable names increased appropriately
    '''
    original_clause_as_list = clause.split(' ')
    assert(original_clause_as_list[0:1] == ['c', 'ind'])
    assert(original_clause_as_list[-1] == '0\n'), original_clause_as_list
    # print(original_clause_as_list[1])

    output_clause_as_list = ['c', 'ind ']

    for char_v in original_clause_as_list[1:-1]:
        assert(np.sign(int(char_v)) > 0)
        v = int(char_v)
        v = v + dup_idx*n_vars
        output_clause_as_list.append('{} '.format(v))
    output_clause_as_list.append('0\n')
    output_clause = ''.join(output_clause_as_list)
    return output_clause

def convert_standardClause_with_duplicates(clause, n_vars, dup_idx, duplication_factor):
    '''
    Inputs:
    - clause: (str) a string representation of a clause
    - n_var: (int) the number of variables in the original formula
    - dup_idx: (int) the duplicate index of this set of variables.  e.g.
        if this is the ith duplicate, variables [0, n_var-1] will be translated
        to [dup_idx*n_var, (dup_idx+1)*n_var - 1] 

    Outputs:
    - output_clause: (str) a string representation of the clause with variable names increased appropriately
    '''
    original_clause_as_list = clause.split(' ')
    assert(original_clause_as_list[-1] == '0\n'), original_clause_as_list
    # print(original_clause_as_list[1])

    output_clause_as_list = []

    for char_v in original_clause_as_list[:-1]:
        var_sign = np.sign(int(char_v))
        v = abs(int(char_v))
        v = v + dup_idx*n_vars
        if var_sign > 0:
            output_clause_as_list.append('{} '.format(v))
        else:
            output_clause_as_list.append('{} '.format(-v))

    extra_var_for_this_duplicate = -(n_vars*duplication_factor + dup_idx + 1)
    output_clause_as_list.append('{} '.format(extra_var_for_this_duplicate))
    output_clause_as_list.append('0\n')
    output_clause = ''.join(output_clause_as_list)
    return output_clause

def parse_cryptominisat_output(output):
    # print('-'*80)
    # print('hi, parse_cryptominisat_output')
    t_begin = time.time()    
    if output is None:
        return 'ERROR', math.nan, math.nan  # Shouldn't reach anywhere,normally
    version = conf.cmsat_major_version
    ctime = math.nan
    nsol = 0
    res_type = "ERROR"
    solutions = []
    global SATISFYING_SOLUTIONS
    global SATISFYING_SOLUTIONS_ARRAY
    global SOLUTION_COUNT
    global SATISFYING_SOLUTIONS_AS_LIST
    global USE_SUM_OF_SATISFYING_SOLUTIONS
    global SUM_OF_SATISFYING_SOLUTIONS
    var_list = conf.var_set

    variable_index_to_val_dict = {}
    for idx, val in enumerate(var_list):
        variable_index_to_val_dict[int(val)] = idx

    addition_time = 0
    solution_construction_time = 0
    append_time = 0
    array_creation_time = 0
    iteration_time = 0
    array_setting_time = 0
    checking_if_ind_var_time = 0
    t_06 = -np.inf
    if version == 5:
        # print('output from cryptominisat5:')
        # print(output)
        # print()            
        current_solution = []
        # print('&'*80)
        # print('&'*80)
        # print('&'*80)
        # print('&'*80)
        # print('cryptominisat output')
        for line in output.split('\n'):
            # print(line)
            if not line.startswith('v ') and len(current_solution) > 0:
                t_00 = 0#time.time()
                solutions.append(current_solution)
                t_01 = 0#time.time()
                append_time += t_01 - t_00
                assert(current_solution[-1] == '0')
                t_02 = 0#time.time()
                solution_as_array = np.zeros(len(var_list))
                t_04 = 0#time.time()                
                for idx, val in enumerate(current_solution):
                    assert(idx + 1 == abs(int(val)) or idx + 1 == len(current_solution)), (idx, val, len(current_solution))
                    t_05 = 0#time.time()
                    iteration_time += t_05 - t_04
                    # if int(val) > 0 and int(val) in var_list:
                    if int(val) > 0 and int(val) in variable_index_to_val_dict:
                        t_06 = 0#time.time()
                        # solution_as_array[idx] = 1.0
                        # solution_as_array[var_list.index(int(val))] = 1.0
                        # assert(var_list.index(int(val)) == variable_index_to_val_dict[int(val)])
                        solution_as_array[variable_index_to_val_dict[int(val)]] = 1.0
                        t_07 = 0#time.time()
                        array_setting_time += t_07 - t_06

                    t_04 = 0#time.time()
                    if t_06 > t_05:
                        checking_if_ind_var_time +=t_06 - t_05
                    else:
                        checking_if_ind_var_time +=t_04 - t_05

                t_03 = 0#time.time()
                array_creation_time += t_03 - t_02   
                # if nsol == 1:
                if nsol >= 1:
                    t0 = 0#time.time()
                    if USE_SUM_OF_SATISFYING_SOLUTIONS:
                        if SUM_OF_SATISFYING_SOLUTIONS is None:
                            SUM_OF_SATISFYING_SOLUTIONS = np.zeros(len(solution_as_array))
                            assert(SOLUTION_COUNT == 0)
                        else:
                            SUM_OF_SATISFYING_SOLUTIONS += solution_as_array
                        SOLUTION_COUNT += 1
                        t1 = 0#time.time()
                        addition_time += t1 - t0
                    elif SATISFYING_SOLUTIONS_AS_LIST:
                        SATISFYING_SOLUTIONS.append(solution_as_array)
                        if SATISFYING_SOLUTIONS_ARRAY is None:
                            SATISFYING_SOLUTIONS_ARRAY = np.zeros((100000, len(solution_as_array)))
                            SATISFYING_SOLUTIONS_ARRAY[SOLUTION_COUNT, :] = solution_as_array                            
                            assert(SOLUTION_COUNT == 0)
                        else:
                            SATISFYING_SOLUTIONS_ARRAY[SOLUTION_COUNT, :] = solution_as_array
                        SOLUTION_COUNT += 1
                    else:
                        if SATISFYING_SOLUTIONS is None:
                            SATISFYING_SOLUTIONS = np.zeros((100000, len(solution_as_array)))
                            SATISFYING_SOLUTIONS[SOLUTION_COUNT, :] = solution_as_array                            
                            assert(SOLUTION_COUNT == 0)
                        else:
                            SATISFYING_SOLUTIONS[SOLUTION_COUNT, :] = solution_as_array
                        SOLUTION_COUNT += 1
                # print("len(SATISFYING_SOLUTIONS):", len(SATISFYING_SOLUTIONS))
                current_solution = []
            t2 = 0#time.time()
            line = line.strip()
            if line.startswith('c Number of solutions found until now:'):
                nsol = int(line.split('now:')[1].strip())
            elif line.startswith('c Total time'):
                ctime = float(line.split(':')[1].strip())
            elif line.startswith('s SAT'):
                res_type = 'SAT'
                nsol += 1
            elif line.startswith('s UNSAT'):
                res_type = 'UNSAT'
            elif line.startswith('s INDET'):
                res_type = 'INDET'
            elif line.startswith('v '):
                current_solution.extend(line.split()[1:])
            t3 = 0#time.time()
            solution_construction_time += t3 - t2
        # print("all solutions:")
        # print(solutions)
        # print(SATISFYING_SOLUTIONS)
    elif version == 2:
        # print('hi!')
        # print(output)
        res_type = 'INDET'
        for line in output.split('\n'):
            line = line.strip()
            if line.startswith('c Number of solutions found until now:'):
                nsol = int(line.split('now:')[1].strip())
            elif line.startswith('c CPU time'):
                ctime = float(line.split(':')[1].strip().split('s')[0].strip())
            elif line.startswith('c SAT'):
                res_type = 'SAT'
                nsol += 1
            elif line.startswith('c UNSAT'):
                res_type = 'UNSAT'
            elif (line.startswith('cryptominisat:') or line ==
                    'Memory manager cannot handle the load. Sorry. Exiting.'):
                log.warning('*** cryptominisat failed!')
                res_type = 'ERROR'
        if nsol > 0:
            res_type = 'SAT'  # Because a 'c UNSAT' is given by CMSAT if maxsol
#                               is not reached.
    else:
        raise Exception('cmsat parsing for this version not yet implemented')
    # print(res_type, nsol, ctime)

    t_end = time.time()    
    # global C_MINISAT_OUPUT_PARSE_TIME
    # print('!'*80)
    # print('cur parse time:', t_end - t_begin)
    # print('solution_construction_time:', solution_construction_time)
    # print('addition_time:', addition_time)
    # print('append_time:', append_time)
    # print('array_creation_time:', array_creation_time)
    # print('iteration_time:', iteration_time)
    # print('array_setting_time:', array_setting_time)
    # print('checking_if_ind_var_time:', checking_if_ind_var_time)
    # C_MINISAT_OUPUT_PARSE_TIME += t_end - t_begin

    return res_type, nsol, ctime


def parse_sharpsat_output(output):
    """
    Parse the output of sharpsat and remove the data.out file it
    generates in the current dir.
    If it finds the number of solutions it is reported in the first
    element of the returning pair, otherwise it is math.nan.
    The second element of the returning pair is 0, to conform with the
    protocol for the last element of the returning tuple from output parsers
    :param output: The sharpsat output. If None then a timeout occured
    :return: A pair of values
    """
    if output is None:
        return math.nan, 0
    dataout_file = os.path.join(os.getcwd(), 'data.out')
    if os.path.isfile(dataout_file):
        os.remove(dataout_file)
    ctime = 0
    nsol = math.nan
    all_lines = output.split('\n')
    for i, line in enumerate(all_lines):
        line = line.strip()
        if line.startswith('# solutions'):
            nsol = int(all_lines[i+1].strip())
    return nsol, ctime


def execute_cmd(sh_cmd, timeout, output_parser, plain_timeout=False, count_time=True):
    """ Do the following:
    * Execute the command
    * Save and parse the output
    :param sh_cmd: The cryptominisat command to execute
    :param timeout: A marginaly larger value than the one passed to cmsat
                    directly
    :param output_parser: function to parse the output of the command. The
                          function should return a tuple with last element the
                          processor time that consumed.
    :param plain_timeout: [boolean] If it is True then it is normal the timeout
                          set for the subprocess to be triggered i.e. there
                          is no other 'internal' timeout in the command
                          that should be triggered earlier. Thus if
                          'plain_timeout' is true and a timeout occurs it is
                          not reported as unrecoverable error.
    :return: The parser output
    """
    global conf, subprocess_cumulative_time
    # check_time()
    # check_time_jdk()
    cmd_l = shlex.split(sh_cmd)
    return_code = math.nan
    # noinspection PyUnusedLocal
    res_type = 'ERROR'
    # noinspection PyUnusedLocal
    nsol = math.nan
    ctime = math.nan
    output = None
    log.debug('Executing cmd: "{}"'.format(sh_cmd))
    try:
        cmd_res = subprocess.run(
            cmd_l, cwd=conf.working_dir, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, timeout=timeout)
        return_code = cmd_res.returncode
        output = cmd_res.stdout.decode('utf-8')
        # print("sh_cmd")
        # print(sh_cmd)
        # print("output")
        # print(output)
    except subprocess.TimeoutExpired as e:
        if not plain_timeout:
            # print('** ', str(e))
            # assert False, 'Wrapper timeout triggered before command finishes.'
            result, n_sol, ctime = 'INDET', None, None #timed out
            return (result, n_sol, ctime)
    finally:
        log.debug('Finished cmd: "{}" with return code:{} time:{}'.format(
                   sh_cmd, return_code, ctime))
        log.debug('TotalTime(s):{}'.format(ctime))
        parsed_results = output_parser(output)
        ctime = parsed_results[-1]
    assert (not math.isnan(ctime))
    subprocess_cumulative_time += ctime

    if count_time:
        global SAT_SOLVER_TIME
        SAT_SOLVER_TIME += ctime
    # print("%"*80)
    # print("ctime low precision:", ctime)

    return parsed_results


def sat_count(formula, time_limit, sol_limit):
    """
    Count the number of solutions of xor+cnf formula.

    :param formula: The full pathname of the formula file
    :param time_limit: Maximum run time (in sec) allowed
    :param sol_limit: Maximum number of solutions to be returned
    :return: a tuple (t, b, n, ctime) such that:
                t: [boolean] has the time_limit been reached?
                b: [boolean] has the max number of solutions be reached?
                n: - min{num_of_solutions, sol_limit}, if 't' is False, or
                   - num of solutions found until time out, if 't' is True
                ctime: The processor time consumed in the subprocess
    """
    # check_time()
    # check_time_jdk()
    cmsat5_cmd = ('{cmsat_exe} --printsol 1  --verb 1'
                  ' --maxtime={timeout}'
                  ' --maxsol={maxsol} {augm_form}')
    cmsat2_cmd = ('{cmsat_exe} --nosolprint  --verbosity=1'
                  ' --gaussuntil=400'
                  ' --maxsolutions={maxsol} {augm_form}')
    t0 = time.time()
    t0prime = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime
    # print("info.ru_utime1:", resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime)
    # print()    
    # print("time_limit:", time_limit)
    if conf.cmsat_major_version == 5:
        sh_cmd = cmsat5_cmd.format(timeout=time_limit, augm_form=formula,
                                   cmsat_exe=conf.cmsat_exe, maxsol=sol_limit)
    elif conf.cmsat_major_version == 2:
        sh_cmd = cmsat2_cmd.format(timeout=time_limit, augm_form=formula,
                                   cmsat_exe=conf.cmsat_exe, maxsol=sol_limit)
    else:
        raise Exception('Unknown cmsat version')

    # print(sh_cmd)
    # sleep(3)

    result, n_sol, ctime = execute_cmd(sh_cmd, time_limit + 5,
                                       parse_cryptominisat_output)
    t1 = time.time()
    t1prime = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime    
    #for getting times with more precision.  
    #cryptominisat rounds to hundredth, this doesn't and looks like it matches cryptominisat well
    ctime = t1prime-t0prime 
    global SAT_SOLVER_TIME_BETTER_PRECISION
    SAT_SOLVER_TIME_BETTER_PRECISION += ctime    
    print('-'*10, "n_sol", n_sol, "ctime", ctime, "remaining_time()", remaining_time(), '-'*10)
    # print("result", result, "n_sol", n_sol, "ctime", ctime, 'manual time:', t1-t0, 'manual time2:', t1prime-t0prime, t1prime, t0prime)
    info = resource.getrusage(resource.RUSAGE_CHILDREN)
    # print(info)
    # print("info.ru_utime2:", info.ru_utime)
    # print("info.ru_utime3:", resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime)

    # print()
    if result == 'ERROR':
        print('sat_count: *** CMSAT ERROR ***')
    # print('#$'*40)
    # print("result:", result, 'time now:', remaining_time())

    rt = remaining_time()
    if result == 'INDET' or rt < 0:
        time_out = True
        b = False
    else:
        time_out = False
        b = n_sol >= sol_limit

    return time_out, b, n_sol, ctime


def sharpsat_count(formula, time_limit):
    """
    Count the number of solutions of a cnf formula by invoking
    sharpsat.
    If sharpsat returns number of solutions equal to 1, the formula may
    also be UNSATISFIABLE. (bug in sharpsat)
    :param formula: The full pathname of the formula file
    :param time_limit: Maximum run time (in sec) allowed
    :return: a triplet (t, n, ctime) such that:
                t: [boolean] has the time_limit been reached without finishing?
                n: number of solutions (if t is False)
                ctime: The processor time consumed in the subprocess
    """
    global subprocess_cumulative_time
    # check_time()
    sh_cmd = '{sharpsat_exe} {form}'.format(form=formula,
                                            sharpsat_exe=conf.sharpsat_exe)
    n_sol, _ = execute_cmd(sh_cmd, time_limit, parse_sharpsat_output,
                           plain_timeout=True)
    # A pessimistic take is that all time is consumed.
    # However if sharpsat terminated earlier then it has successfully
    # counted solutions so we don't care about time any more.
    ctime = time_limit
    subprocess_cumulative_time += ctime
    return n_sol is math.nan, n_sol, ctime


def algorithm1(n_constr, n_iter, sampleRepeats=1, method='original', VERBOSE=False):
    """
    Determine if 2^(n_constr) is a lower bound with error probability
    exp(-n_iter/8)
    :param n_constr: Number of XOR's to add to the formula
    :param n_iter: Number of iterations ('t' in the paper)
    :return: True or False

    Outputs:
    - bound_holds: (bool) True -> the bound holds
    - time_out: (bool) True -> experienced timeout (bound_holds must be False)
    """
    print()
    # print("sampleRepeats:", sampleRepeats)
    global SOLUTION_COUNT
    global SAT_SOLVER_TIME_BETTER_PRECISION_PARALLEL
    print("algorithm 1 called with", n_constr, "constraints (there are", SOLUTION_COUNT, "satisfying_solutions found currently) and n_iter =", n_iter)
    marginal_products = []
    solution_counts = []

    global SAT_SOLVER_TIME
    begin_time = SAT_SOLVER_TIME_BETTER_PRECISION
    assert(n_constr > 0)
    n_var = len(conf.var_set)
    if n_constr >= n_var:
        print("n_constr >= n_var, n_constr=", n_constr, "n_var=", n_var)
        bound_holds = False
        time_out = False
        return (bound_holds, time_out)

    assert(conf.alg1_maxsol == 4)
    threshold = conf.alg1_threshold_coeff * n_iter
    assert(threshold == n_iter*conf.alg1_maxsol/2)
    z = 0  # The capital Z in the paper
    bound_holds = False
    # print("n_iter:", n_iter, "threshold", threshold)
    list_of_i_effectives = []
    all_duplicate_constraints = []
    total_individual_solution_count = 0
    for i in range(1, n_iter+1):
        cur_solution_counts = []
        cur_runtimes = []
        for repeat in range(conf.extra_configs['sum_of_T_solutions']):
            t1 = time.process_time()
            constraints, product_of_marginals, i_effective, fail_apriori = generateConstraints(n_constr, conf.degn, sampleRepeats=sampleRepeats, method=method)
            all_duplicate_constraints.append(constraints)
            list_of_i_effectives.append(i_effective)
            marginal_products.append(product_of_marginals)


            if fail_apriori:
                # t: [boolean] has the time_limit been reached?
                # b: [boolean] has the max number of solutions be reached?
                # n: - min{num_of_solutions, sol_limit}, if 't' is False, or
                #    - num of solutions found until time out, if 't' is True
                # ctime: The processor time consumed in the subprocess
                t = False
                b = False
                n = 0
                ctime = 0.0
                t2 = time.process_time()
                solution_counts.append(0)
                cur_solution_counts.append(0)

            else:
                global SAT_SOLVER_TIME_BETTER_PRECISION
                rt = conf.total_max_time - SAT_SOLVER_TIME_BETTER_PRECISION
                # print('about to call setup_augmented_formula, remaining time:', remaining_time())


                aug_form = setup_augmented_formula(n_constr, constraints)
                t2 = time.process_time()
                # max_time = min(conf.alg1_loop_timeout, remaining_time() + 1)
                # max_time = min(conf.alg1_loop_timeout, remaining_time_jdk())
                max_time = min(conf.alg1_loop_timeout, remaining_time())
                time_out, b, n, ctime = sat_count(aug_form, max_time, conf.alg1_maxsol*conf.extra_configs['sum_of_T_solutions'])
                if not conf.keep_cnf:
                    os.remove(aug_form)

                if time_out:
                    assert(bound_holds == False)
                    cur_solution_counts.append(0)
                    cur_runtimes.append(max_time)
                    # return (bound_holds, time_out)

                # t, b, n, ctime = sat_count(aug_form, max_time, 16) #jdk temp edit
                if VERBOSE:
                    print("time:", ctime, "solutions found:", n, end=' ')
                solution_counts.append(n)
                cur_solution_counts.append(n)
                cur_runtimes.append(ctime)
                # if np.sum(cur_solution_counts)/conf.extra_configs['sum_of_T_solutions'] >= 4:
                #     break
                # total_individual_solution_count += n

        SAT_SOLVER_TIME_BETTER_PRECISION_PARALLEL += np.max(cur_runtimes)
        assert(len(cur_solution_counts) == conf.extra_configs['sum_of_T_solutions'])
        n = min(conf.alg1_maxsol, np.sum(cur_solution_counts)/conf.extra_configs['sum_of_T_solutions'])
        print("cur_solution_counts:", cur_solution_counts)

        # print('!'*80)
        # print(t, b, n, ctime)
             
        # log.info('Algorithm1:{itr}:{n_constr}:{n_sol}:{setup:.2f}:{run:.2f}:'
        #          '{timeout}:{hit_bound}'
        #          ''.format(itr=i, n_constr=n_constr, n_sol=n,
        #                    setup=t2 - t1, run=ctime, timeout=t, hit_bound=b))
        z += n
        if VERBOSE:
            print("Z=", z)
        if z >= threshold:
            bound_holds = True
            break #jdk temp edit

    # positive_constraints = get_non_zero_constraints(all_duplicate_constraints)
    # aug_form = setup_duplicated_augmented_formula(positive_constraints)
    if VERBOSE:
        if n_iter == 24:
            pass
            # time_out, b, n, ctime = sat_count(aug_form, max_time, 48)
            # print('!'*80)
            # print("with duplication, time:", ctime, "solutions found:", n)
            # print("total_individual_solution_count:", total_individual_solution_count)
            # print()


        else:
            pass
            # time_out, b, n, ctime = sat_count(aug_form, max_time, conf.alg1_maxsol*n_iter)
            # print('!'*80)
            # print("with duplication, time:", ctime, "solutions found:", n)
            # print("total_individual_solution_count:", total_individual_solution_count)
            # print()        
            # assert(n >= total_individual_solution_count), (n, total_individual_solution_count)



    end_time = SAT_SOLVER_TIME_BETTER_PRECISION
    if VERBOSE:
        print("mean i_effective:", np.mean(list_of_i_effectives))
        print("algorithm 1 spent", end_time - begin_time, "seconds calling the sat solver")
        print("correlation coefficient between product of marginals and satisfying solutions")
        print(np.corrcoef(marginal_products, solution_counts))
        print("sum of solutions found in first 24 problems:", sum(solution_counts[:24]))
        print("sum of solutions found in all", len(solution_counts), "problems:", sum(solution_counts))
        print("number of SAT found in first 24 problems:", 24 - solution_counts[:24].count(0))

        if len(solution_counts) == 240:
            sat_count_group_by_10 = 0
            sum_of_solutions_pick_groups_of_ten = 0
            for grp_idx in range(0, 240, 10):
                cur_max_idx = marginal_products[grp_idx:grp_idx+10].index(max(marginal_products[grp_idx:grp_idx+10]))
                sum_of_solutions_pick_groups_of_ten += solution_counts[grp_idx:grp_idx+10][cur_max_idx]
                if solution_counts[grp_idx:grp_idx+10][cur_max_idx] > 0:
                    sat_count_group_by_10 += 1
            print("sum of solutions from max product marginals grouped by 10:", sum_of_solutions_pick_groups_of_ten)
            print("number of SAT from max product marginals grouped by 10:", sat_count_group_by_10)
        else:
            print("len(solution_counts):", len(solution_counts))
        # marginal_products, solution_counts = (list(x) for x in zip(*sorted(zip(marginal_products, solution_counts))))
        print("marginal_products:", marginal_products)
        print("solution_counts:", solution_counts)
        marginal_products, solution_counts = [list(x) for x in zip(*sorted(zip(marginal_products, solution_counts), key=lambda pair: pair[0]))]
        print("sum of solutions found in 24 problems with largest marginal_products:", sum(solution_counts[-24:]))
        print("number of SAT from 24 problems with largest marginal_products:", 24-solution_counts[-24:].count(0))
        print()

    time_out = False
    return (bound_holds, time_out)


def algorithm2(lb, delta: float=math.nan, num_iter: int=math.nan):
    """
    Calculate an upper bound which is correct with probability 1-delta
    given a correct lower bound.
    :param lb: Floor of the binary log of the lower bound ('l' in the paper)
    :param delta: error probability
    :param num_iter: bypass the rigorous computation and give num of iterations
    :return: (a,b) such that a*2^b is a rigorous upper bound
    """
    assert lb >= 1
    force_long = False
    # Only one of delta or num_iter can be specified
    assert((not math.isnan(delta) and math.isnan(num_iter)) or
           (math.isnan(delta) and not math.isnan(num_iter)))
    if not math.isnan(num_iter):
        n_iter = num_iter
    else:
        bst = 1
        if needLong(lb):
            force_long = True
        else:
            bst = get_boost(lb)  # B in the paper
            if bst == 1:
                force_long = True
        n_iter = int(math.ceil(8*(bst+1)*math.log(1/float(delta))))
    z = 0
    for j in range(1, n_iter+1):
        t1 = time.process_time()
        if not need_dummy_var(conf.degn, len(conf.var_set), lb):
            constraints = generateConstraints(lb, conf.degn, force_long)
        else:
            constraints = generateConstraints(lb, conf.degn, force_long,
                                              dummy_vars=1)
        aug_form = setup_augmented_formula(lb, constraints)
        t2 = time.process_time()
        max_time = min(conf.alg2_loop_timeout, remaining_time())
        t, b, n, ctime = sat_count(aug_form, max_time, conf.alg2_maxsol)
        log.info('Algorithm1:{itr}:{n_constr}:{n_sol}:{setup:.2f}:{run:.2f}:'
                 '{timeout}:{hit_bound}'
                 ''.format(itr=j, n_constr=lb, n_sol=n,
                           setup=t2 - t1, run=ctime, timeout=t, hit_bound=b))
        if not conf.keep_cnf:
            os.remove(aug_form)
        if t or b:
            raise SolverBoundsExceededException(
                'Algorithm2: timeout={}, cutoff={}'.format(t, b))
        z += n
    b = z / float(n_iter)
    if not need_dummy_var(conf.degn, len(conf.var_set), lb):
        e = lb + 1
    else:
        e = lb
    return b, e


def F2_algorithm(lb_in, ub_in, delta_in: float, theta: float):
    """
    Calculates an estimate for the number of solutions with a specified
    confidence and accuracy
    :param lb_in: binary logarithm of the lower bound
    :param ub_in: binary logarithm of the upper bound
    :param delta_in: confidence
    :param theta: probability of error
    :return: a, b such that the estimate is a * 2^b
    """
    check_time()
    log.info('Starting F2_algorithm( lb={:.2f}, ub={:.2f}, delta={}, theta={}'
             ''.format(lb_in, ub_in, delta_in, theta))
    if lb_in < 2 - math.log2(delta_in):  # Small LB phase
        t, b, n, ctime = sat_count(conf.formula_filename, conf.total_max_time,
                                   4 / delta_in)
        if t:
            raise SolverBoundsExceededException(
                'F2_algorithm(Small LB phase): Uncontainable Solver Timeout'
                ' occurred')
    lb = int(math.floor(math.log2(delta_in) + lb_in - 2))
    ub = int(math.ceil(ub_in))
    assert ub - 1 > lb  # The range(lb, ub-1) has at least one element
    bst = get_boost(list(range(lb, ub-1)))

    delta = min(delta_in, 1/3)
    xi = 8 / delta
    b = math.ceil(xi + 2*(xi + xi**2*(bst-1)))
    n_iter = int(math.ceil((2*b**2/9)*math.log(5.0/theta)))

    z_list = [0 for _ in range(lb, ub+1)]

    f2_alg_loop_timeout = remaining_time() / min(5, n_iter*(ub-lb+1))
    force_long = needLong(lb)
    for j in range(1, n_iter+1):
        if not need_dummy_var(conf.degn, len(conf.var_set), lb):
            constraints = generateConstraints(lb, conf.degn, force_long)
        else:
            constraints = generateConstraints(lb, conf.degn, force_long,
                                              dummy_vars=1)
        for i in range(lb, ub+1):
            aug_form = setup_augmented_formula(i, constraints)
            max_time = min(f2_alg_loop_timeout, remaining_time())
            tp, _, y, ctime = sat_count(aug_form, max_time, b)
            if not conf.keep_cnf:
                os.remove(aug_form)
            if tp:
                raise SolverBoundsExceededException(
                    'F2_algorithm(Main phase): Uncontainable Solver Timeout '
                    'occurred')
            z_list[i-lb] += y

    thres = n_iter * (1 - delta) * (4 / delta)
    kl = [i for i in range(lb, ub+1) if z_list[i - lb] > thres]
    if not kl:
        return -1, -1
    else:
        if not need_dummy_var(conf.degn, len(conf.var_set), lb):
            j = kl[-1]
        else:
            j = kl[-1] - 1
        return z_list[j]/n_iter, j


def print_sharpsat_warning():
    msg = """
*** Important Notice: The result is produced by
SharpSAT (https://github.com/marcthurley/sharpSAT/) 
which ignores header lines starting with 'c ind' in the CNF input file.
If the variables specified there represent an INDEPENDENT SUPPORT then you
don't need to do anything. Otherwise ignore the result and run F2 again, 
adding the '--skip-sharpsat' option to the command line.\
"""
    print(msg)


# Boost Calculation related code below
def rfloor(l, n, i):
    return int(sp.floor(l * n / i))


def i0(l, n, i):
    return int(-l * n + i * sp.floor(l * n / i) + i)


def i1(l, n, i):
    return int(l * n - i * sp.floor(l * n / i))


def pl(r):
    return [(sp.binomial(r, 2 * j).evalf(), 2 * j) for j in range(0, r + 1)]


def ql(r):
    return [(sp.binomial(r + 1, 2 * j).evalf(), 2 * j)
            for j in range(0, r + 1)]


def h(t):
    return -t * sp.log(t, 2) - (1 - t) * sp.log(1 - t, 2)


def zp(i, n):
    # Solve h(x)-(i-1)/n == 0 with starting point x_0=0.000000001 where
    # h(x) the binary entropy function: h(x):=-x*log2(x)-(1-x)*log2(1-x)
    t = scipy.optimize.fsolve(
        lambda x: (-x) * math.log2(x) -(1-x) * math.log2(1-x) - (i - 1) / n,
        0.000000001
    )[0] * n
    return int(math.ceil(t))


def power_productl(cel1, ex1, cel2, ex2, order):
    assert order > 0
    assert ex1 > 0
    # Sort the list by increasing exponent
    cel1s = sorted(cel1, key=lambda x: x[1])
    cel2s = sorted(cel2, key=lambda x: x[1])
    # Prune the exponents above 'order' in the initial list
    cel1s = [x for x in cel1s if x[1] < order]
    cel2s = [x for x in cel2s if x[1] < order]
    # Coefficients of powers x^0, x^1, ... x^\(order-1\)
    coef_accumulators = [0] * order
    # Initialize accumulators by first the first factor
    for c, e in cel1s:
        if e < order:
            coef_accumulators[e] += c
    ex1 -= 1
    while ex1 > 0:
        coef_accumulators_snapshot = coef_accumulators.copy()
        for c, e in cel1s:
            for i in range(0, order - e):
                if coef_accumulators_snapshot[i] != 0:
                    if e == 0:
                        coef_accumulators[i] *= c
                    else:
                        coef_accumulators[i + e] += coef_accumulators_snapshot[
                                                        i] * c
        ex1 -= 1
    while ex2 > 0:
        coef_accumulators_snapshot = coef_accumulators.copy()
        for c, e in cel2s:
            for i in range(0, order - e):
                if coef_accumulators_snapshot[i] != 0:
                    if e == 0:
                        coef_accumulators[i] *= c
                    else:
                        coef_accumulators[i + e] += coef_accumulators_snapshot[
                                                        i] * c
        ex2 -= 1
    return coef_accumulators


def codewordsl(l, n, i, dlist):
    ordr = max(dlist) * l + 1
    e1 = pl(rfloor(l, n, i))
    e2 = ql(rfloor(l, n, i))
    expansion = power_productl(e1, i0(l, n, i), e2, i1(l, n, i), ordr)
    cwlist = []
    for d in dlist:
        c = expansion[d * l]
        b = sp.binomial(n, d).evalf() / sp.binomial(n * l, d * l).evalf()
        cwlist.append(b * c)
    return cwlist


def boostl(l, n, i):
    zeta = zp(i, n)
    if zeta < 2:
        return 1
    d_l = list(range(1, zeta))
    s1 = sum(codewordsl(l, n, i, d_l))
    s2 = sum([sp.binomial(n, d2).evalf() for d2 in d_l])
    return s1 / s2


def boost(l, n, i):
    # if the constraint length is >= n/2 then we use long constraints
    if math.ceil((n * l) / i) >= n/2:
        return 1

    if l * n / i == rfloor(l, n, i) and rfloor(l, n, i) % 2 == 0:
        return boostl(l, n + 1, i) * 2 ** i
    else:
        return boostl(l, n, i) * 2 ** i

# #####


if __name__ == '__main__':
    # time_out, solution_count, sharp_sat_time = sharp_sat_call_from_python(problem_name='54.sk_12_97.cnf.gz.no_w.cnf', time_limit=2)
    # print(time_out, solution_count, sharp_sat_time)
    # exit(0)

    extra_configs = {
                     #density of ones in constraint matrix for method = 'iid'
                     'f': .5,
                     #when sampling biregular matrices such that each constraint has a good marginal,
                     #how do we deal with the 'problem constraints' at the end?
                     # - 'iid': sample them iid, give up biregular
                     # - 'keep_biregular': leave them as biregular, give up good marginals
                     'biregular_marginal_problem_constraint': 'keep_biregular',
                     #boolean
                     # - True: sample b_i according to the constraint's marginal
                     # - False: sample b_i with probability .5
                     'adaptive_b': True,
                     'sum_of_T_solutions':10,
                    }


    # method:'original', 'bi_regular_marginals_per_constraint', 'bi_regular_marginals_joint_constraint', 
    # 'iid', 'bi_regular_order_vars_by_marginals', 'bi_regular_order_vars_by_double_marginals',
    # 'bi_regular_order_vars_by_marginals_randomChunks'


    # cProfile.run("find_lower_bound_call_from_python(problem_name='75-22-6-q.cnf.gz.no_w.cnf',\
    #     random_seed=0, var_degree=1, method='original', extra_configs=extra_configs)")
    # exit(0)

#s13207a_15_7 is large and original 1 is better than 1.5, runtime 1.5: 3.9239999999999995 bound 1.5: 536.0 runtime 1: 1.8480000000000008 bound 1: 695.0 
#min-12 is large and original 1 is better than 1.5, runtime 1.5: 1.6480000000000006 bound 1.5: 317.0 runtime 1: 1.0079999999999996 bound 1: 355.0
    

    global CONSTRAINT_GENERATION_TIME
    CONSTRAINT_GENERATION_TIME = 0   

    global C_MINISAT_OUPUT_PARSE_TIME
    C_MINISAT_OUPUT_PARSE_TIME = 0  

    t_begin = time.time()
    # before cProfile was just using this!
    lb, sat_solver_time, timeout, parallel_runtime = find_lower_bound_call_from_python(problem_name='90-34-3-q.cnf.gz.no_w.cnf',\
    # lb, sat_solver_time, timeout, parallel_runtime = find_lower_bound_call_from_python(problem_name='17A-6.cnf.gz.no_w.cnf',\


    # lb, sat_solver_time, timeout = find_lower_bound_call_from_python(problem_name='50-12-5-q.cnf.gz.no_w.cnf',\
    # lb, sat_solver_time, timeout = find_lower_bound_call_from_python(problem_name='75-14-2-q.cnf.gz.no_w.cnf',\
    # lb, sat_solver_time, timeout = find_lower_bound_call_from_python(problem_name='75-12-5-q.cnf.gz.no_w.cnf',\
    # lb, sat_solver_time, timeout = find_lower_bound_call_from_python(problem_name='75-12-3-q.cnf.gz.no_w.cnf',\
    # lb, sat_solver_time, timeout = find_lower_bound_call_from_python(problem_name='min-12.cnf.gz.no_w.cnf',\
    # lb, sat_solver_time, timeout = find_lower_bound_call_from_python(problem_name='s1423a_15_7.cnf.gz.no_w.cnf',\
    # lb, sat_solver_time, timeout = find_lower_bound_call_from_python(problem_name='01B-5.cnf.gz.no_w.cnf',\
    # lb, sat_solver_time, timeout = find_lower_bound_call_from_python(problem_name='50-12-10-q.cnf.gz.no_w.cnf',\
    # lb, sat_solver_time, timeout = find_lower_bound_call_from_python(problem_name='sat-grid-pbl-0015.cnf',\
    # lb, sat_solver_time, timeout = find_lower_bound_call_from_python(problem_name='tire-1.cnf',\
    # lb, sat_solver_time, timeout = find_lower_bound_call_from_python(problem_name='hypercube.cnf',\
     # random_seed=1, var_degree=1.3, method='original', extra_configs=extra_configs)
     random_seed=2, var_degree=1, method='bi_regular_order_vars_by_marginals_randomChunks', extra_configs=extra_configs)
     # random_seed=0, var_degree=1.5, method='original', extra_configs=extra_configs)
    t_end = time.time()
    print("lower bound:", lb, "sat_solver_time:", sat_solver_time, "timeout:", timeout, "parallel_runtime:", parallel_runtime)


    # cProfile.run("find_lower_bound_call_from_python(problem_name='90-34-3-q.cnf.gz.no_w.cnf',\
    #     random_seed=0, var_degree=1, method='bi_regular_order_vars_by_marginals_randomChunks', extra_configs=extra_configs)")

    print("lower bound:", lb, "sat_solver_time:", sat_solver_time, "timeout:", timeout, "parallel_runtime:", parallel_runtime)  
    print("CONSTRAINT_GENERATION_TIME:", CONSTRAINT_GENERATION_TIME)  
    print("C_MINISAT_OUPUT_PARSE_TIME:", C_MINISAT_OUPUT_PARSE_TIME)  
    print("total time:", t_end - t_begin)
    exit(0)


    global SAT_SOLVER_TIME
    global SAT_SOLVER_TIME_BETTER_PRECISION
    global CONSTRAINTS_CALLED_COUNT
    global SATISFYING_SOLUTIONS
    global SOLUTION_COUNT
    global SATISFYING_SOLUTIONS_AS_LIST
    SAT_SOLVER_TIME = 0  
    SAT_SOLVER_TIME_BETTER_PRECISION = 0  
    CONSTRAINTS_CALLED_COUNT = 0
    SOLUTION_COUNT = 0
    SATISFYING_SOLUTIONS_AS_LIST = False
    # with open('/atlas/u/jkuck/XORModelCount/SATModelCount/satisfying_solutions_25sat-grid-pbl-0015', 'rb') as file:
    #     SATISFYING_SOLUTIONS = pickle.load(file)

    # SATISFYING_SOLUTIONS = np.array(SATISFYING_SOLUTIONS)
    SATISFYING_SOLUTIONS = []
    main()
