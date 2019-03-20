from pyeda.parsing.dimacs import parse_cnf
from pyeda.boolalg.expr import ast2expr
# from pyeda.expr import DimacsCNF
from pyeda.boolalg.expr import expr2dimacscnf
import pyeda
import numpy as np




def redistribute_SAT_sols_oneDimesion(input_cnf, split_var):
    '''
    Inputs:
    - input_cnf: (pyeda expression)
    - split_var: (int) positive integer representing the variable to split on

    Outputs:
    - output_cnf: (pyeda expression) CNF with the same #SAT count, but solutions
        redistributed around split_var

    '''
    assert(split_var > 0)
    litmap, nvars, clauses = input_cnf.encode_cnf()

    #first half of the new problem with split_var = True
    variable_random_flips = np.random.binomial(n=1, p=.5, size=nvars+1)
    for idx, value in enumerate(variable_random_flips):
        if value == 0:
            variable_random_flips[idx] = -1

    half1_orig_cnf_as_expr = litmap[split_var]
    for clause in clauses:
        cur_vars = []
        clause_always_true = False
        for var in clause:
            # if var == -split_var: #simplify by removing from clause
            #     continue
            # elif var == split_var: #this clause is covered by the clause that is simply split_var
            #     clause_always_true = True
            #     break
            # if -var in clause:
            #     continue
            # elif abs(var) == split_var:
            #     cur_vars.append(litmap[var])
            if abs(var) == split_var:
                cur_vars.append(litmap[var])            
            else:
                cur_vars.append(litmap[variable_random_flips[abs(var)]*var])
                # cur_vars.append(litmap[var])
            # print(var, end=' ')
        if clause_always_true:
        # if clause_always_true or len(cur_vars) == 0:        
            continue
        cur_or = pyeda.boolalg.expr.Or(*cur_vars)
        # print("cur_or:", cur_or)
        # if half1_orig_cnf_as_expr is None:
        #   half1_orig_cnf_as_expr = cur_or
        # else:
        half1_orig_cnf_as_expr = pyeda.boolalg.expr.And(*[half1_orig_cnf_as_expr, cur_or])  

    #second half of the new problem with split_var = False
    variable_random_flips2 = np.random.binomial(n=1, p=.5, size=nvars+1)
    for idx, value in enumerate(variable_random_flips2):
        if value == 0:
            variable_random_flips2[idx] = -1

    half2_orig_cnf_as_expr = litmap[-split_var]
    for clause in clauses:
        cur_vars = []
        clause_always_true = False
        for var in clause:
            # if var == -split_var: #this clause is covered by the clause that is simply -split_var
            #     clause_always_true = True
            #     break               
            # elif var == split_var: #simplify by removing from clause
            #     continue
            # if -var in clause:
            #     continue
            # elif abs(var) == split_var:
            #     cur_vars.append(litmap[var])
            if abs(var) == split_var:
                cur_vars.append(litmap[var])
            else:
                cur_vars.append(litmap[variable_random_flips2[abs(var)]*var])
                # cur_vars.append(litmap[variable_random_flips[var]*var])
                # cur_vars.append(litmap[var])
            # print(var, end=' ')
        # assert(len(cur_vars)>0)
        # if clause_always_true or len(cur_vars) == 0:
        if clause_always_true:
            continue
        cur_or = pyeda.boolalg.expr.Or(*cur_vars)
        # print("cur_or:", cur_or)
        # if half2_orig_cnf_as_expr is None:
        #   half2_orig_cnf_as_expr = cur_or
        # else:
        # print("clause:", clause, "cur_or:", cur_or, "half2_orig_cnf_as_expr:", half2_orig_cnf_as_expr)
        half2_orig_cnf_as_expr = pyeda.boolalg.expr.And(*[half2_orig_cnf_as_expr, cur_or])


    # print("half1_orig_cnf_as_expr:", half1_orig_cnf_as_expr)
    # print("half2_orig_cnf_as_expr:", half2_orig_cnf_as_expr)
    total_expression = pyeda.boolalg.expr.Or(*[half1_orig_cnf_as_expr, half2_orig_cnf_as_expr])
    # print("total_expression:", total_expression)
    output_cnf = total_expression.to_cnf()
    return output_cnf

def redistribute_SAT_sols_allDimesions(cnf_input_file, cnf_output_file):
    with open(cnf_input_file, 'r') as f:
        cnf_as_string = f.read()
        cnf_abstract_syntax_tree = parse_cnf(cnf_as_string)

    cnf_expression = ast2expr(cnf_abstract_syntax_tree)
    print(cnf_expression)
    print()

    litmap, nvars, clauses = cnf_expression.encode_cnf()

    # for var in range(1, nvars + 1):
    for var in range(1, 10):
        output_cnf = redistribute_SAT_sols_oneDimesion(input_cnf=cnf_expression, split_var=1)   
        litmap, nvars, clauses = output_cnf.encode_cnf()
        print('split on var', var, 'now we have', len(clauses), 'clauses')
        # for clause in clauses:
        #     for var in clause:
        #         print(var, end=' ')         
        #     print()
        # print()
        cnf_expression = output_cnf

    cnf_string = expr2dimacscnf(output_cnf)
    # print(cnf_string[1])
    # print(type(cnf_string[1]))
    with open(cnf_output_file, 'w') as ofile:
        ofile.write(str(cnf_string[1]))
        
    return output_cnf



if __name__ == '__main__':
    problem_name = 'c432.isc'
    problem_name = 'test_cnf2.cnf'
    problem_name = 'hypercube.cnf'
    # cnf_input_file = '/atlas/u/jkuck/low_density_parity_checks/SAT_problems_cnf/%s' % 'test_cnf.cnf'
    # cnf_input_file = '/atlas/u/jkuck/low_density_parity_checks/SAT_problems_cnf/%s' % 'hypercube.cnf'
    cnf_input_file = '/atlas/u/jkuck/low_density_parity_checks/SAT_problems_cnf/%s' % problem_name

    cnf_output_file = '/atlas/u/jkuck/low_density_parity_checks/SAT_problems_cnf/redistributed_%s' % problem_name

    output_cnf = redistribute_SAT_sols_allDimesions(cnf_input_file, cnf_output_file)
    print("output_cnf:", output_cnf)
    litmap, nvars, clauses = output_cnf.encode_cnf()
    print("HIIHIHIHIHI")
    print("litmap:", litmap)
    print("nvars:", nvars)
    print("clauses:", clauses)

    for clause in clauses:
        for var in clause:
            print(var, end=' ')         
        print()

    exit(0)

    # /atlas/u/jkuck/sharpSAT/build/Release/sharpSAT /atlas/u/jkuck/low_density_parity_checks/SAT_problems_cnf/hypercube.cnf
    # /atlas/u/jkuck/sharpSAT/build/Release/sharpSAT /atlas/u/jkuck/low_density_parity_checks/SAT_problems_cnf/test_redistribute_cnf.cnf
    # /atlas/software/cluster/cryptominisat/build/cryptominisat5 --printsol 0  --verb 1 --maxtime=8 --maxsol=1 /atlas/u/jkuck/low_density_parity_checks/SAT_problems_cnf/test_redistribute_cnf.cnf

    '''
python3 /atlas/u/jkuck/F2_clean/src/python/f2.py --cmsat-exe /atlas/u/jkuck/software/cryptominisat_BIRD/cryptominisat-5.6.6/build/cryptominisat5 --cmsat-version 5 --skip-sharpsat --sharpsat-exe /atlas/u/jkuck/sharpSAT/build/Release/sharpSAT --mode lb --var-degree 1 --max-time 3600 --error-prob .05  --random-seed 0 /atlas/u/jkuck/low_density_parity_checks/SAT_problems_cnf/test_redistribute_cnf.cnf
python3 /atlas/u/jkuck/F2_clean/src/python/f2.py --cmsat-exe /atlas/u/jkuck/software/cryptominisat_BIRD/cryptominisat-5.6.6/build/cryptominisat5 --cmsat-version 5 --skip-sharpsat --sharpsat-exe /atlas/u/jkuck/sharpSAT/build/Release/sharpSAT --mode lb --var-degree 3 --max-time 3600 --error-prob .05  --random-seed 0 /atlas/u/jkuck/low_density_parity_checks/SAT_problems_cnf/c432.isc

    '''
    cnf_file = '/atlas/u/jkuck/low_density_parity_checks/SAT_problems_cnf/%s' % 'test_cnf.cnf'
              
    with open(cnf_file, 'r') as f:
        cnf_as_string = f.read()
        abstract_syntax_tree = parse_cnf(cnf_as_string)

    expression = ast2expr(abstract_syntax_tree)
    print(expression)
    print()

    litmap, nvars, clauses = expression.encode_cnf()
    print("litmap:", litmap)
    print("nvars:", nvars)
    print("clauses:", clauses)
    print()

    half1_orig_cnf_as_expr = None
    for clause in clauses:
        cur_vars = []
        for var in clause:
            cur_vars.append(litmap[var])
            print(var, end=' ')
        cur_or = pyeda.boolalg.expr.Or(*cur_vars)
        print("cur_or:", cur_or)
        if half1_orig_cnf_as_expr is None:
            half1_orig_cnf_as_expr = cur_or
        else:
            half1_orig_cnf_as_expr = pyeda.boolalg.expr.And(*[half1_orig_cnf_as_expr, cur_or])
        print("half1_orig_cnf_as_expr:", half1_orig_cnf_as_expr)
            
        print()

    print()
    print("half1_orig_cnf_as_expr:", half1_orig_cnf_as_expr)
    exit(-1)    



    for clause in expression.iter_dfs():
        if clause == litmap[1]:
            print('HI!!')
        # if litmap[1] in clause:
        #   print('HI2!!')

        print(clause)
        print(type(clause))
        print(isinstance(clause, pyeda.boolalg.expr.AndOp))

    print()
    exit(0)
    print(expr2dimacscnf(expression))

    print(clauses)  