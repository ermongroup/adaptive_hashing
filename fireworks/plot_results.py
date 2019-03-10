import numpy as np
import os
from collections import defaultdict
import matplotlib
import matplotlib.pyplot as plt
import math

def get_results(results_folder, random_seed):
    '''
    Outputs:
    - runtimes_for_each_method: (dictionary)
        key: (string) method name
        value: (list of floats) runtimes for the method
               the ith entry of each list contains the runtime for the same method
    - lower_bounds_for_each_method: (dictionary)
        key: (string) method name
        value: (list of floats) lower bounds for the method
        the ith entry of each list contains the lower bound for the same method
    '''
    number_of_problems_sharp_sat_solved_fast = 0
    number_of_unfinished_problems = 0
    number_of_finished_problems = 0

    runtimes_for_each_method = defaultdict(list)
    lower_bounds_for_each_method = defaultdict(list)

    for filename in os.listdir(results_folder):
        if filename.endswith(".txt"): 
            problem_name = filename.split('.txt')[0]
            # print(problem_name)
            results_for_each_method = []
            results_for_each_method_trial1 = []
            results_for_each_method_trial0 = []
            with open(results_folder + '/' + filename,'r') as f:
                for line in f:
                    if line.split(' ')[0] == "repeats_of_randomized_hashing_methods:":
                        continue
                    else:
                        # print(line)
                        method = line.split(' ')[0]
                        timeout = line.split(' ')[2]
                        if method == 'sharpSAT':
                            # print()
                            # print('-'*80)
                            # print(problem_name)
                            # print(float(line.split(' ')[4]))
                            # print(line)
                            if float(line.split(' ')[4]) == 0 or line.split(' ')[4] == 'nan':
                                lower_bound = -np.inf
                            else:
                                lower_bound = math.log(np.long(line.split(' ')[4]))/math.log(2)         
                            if lower_bound > 20000:
                                lower_bound = -np.inf               
                        else:
                            cur_random_seed = int(line.split(' ')[8])
                            if line.split(' ')[4] == 'False':
                                lower_bound = 0
                            else:
                                lower_bound = float(line.split(' ')[4])
                        if timeout == 'True':
                            sat_solver_time = 1000.0
                            if method == 'biregular_order_vars_by_marginals_assignmentProblem_variable_degree_1' and float(line.split(' ')[6]) < 100:
                                # print("problem", problem_name, "timed out for method biregular_order_vars_by_marginals_assignmentProblem_variable_degree_1")
                                print(problem_name, float(line.split(' ')[6]))
                            continue
                        else:
                            assert(timeout == 'False')
                            sat_solver_time = float(line.split(' ')[6])
                        if method == 'sharpSAT' or cur_random_seed == random_seed:
                            results_for_each_method.append((method, timeout, lower_bound, sat_solver_time))
                        if method != 'sharpSAT' and cur_random_seed == 1:
                            results_for_each_method_trial1.append((method, timeout, lower_bound, sat_solver_time))
                        elif method != 'sharpSAT' and cur_random_seed == 0:
                            results_for_each_method_trial0.append((method, timeout, lower_bound, sat_solver_time))

                for (method, timeout, lower_bound, sat_solver_time) in results_for_each_method:
                    if method == 'long_iid_.5' and lower_bound == 2:
                        # print("problem_name:", problem_name)
                        for (method1, timeout1, lower_bound1, sat_solver_time1) in results_for_each_method:
                            # print(method1)
                            if method1 == 'sharpSAT':
                                true_count = lower_bound1
                        print()
                        print('-'*80)
                        print('iid .5 returned bound of 2 true count =', true_count)
                        print(results_folder + '/' + filename)
                        # for line1 in f:
                        #     print(line1)
                        print('-'*80)
                        print()

            # if len(results_for_each_method_trial0) == 7 and len(results_for_each_method_trial1) != 7:
            #     print("trial 0 completed, but not trial 1 for problem:", problem_name)


            method1_found = False
            method2_found = False
            method3_found = False
            method4_found = False
            for cur_method_results in results_for_each_method:
                if cur_method_results[0] == 'biregular_variable_degree_1_Tsol_1':
                    method1_found = True
                if cur_method_results[0] == 'biregular_variable_degree_1_Tsol_10':
                    method2_found = True
                if cur_method_results[0] == 'biregular_order_vars_by_marginals_assignmentProblem_variable_degree_1_Tsol_10':
                    method3_found = True
                    if cur_method_results[1] == 'True':
                        print("problem", problem_name, "timed out for adaptive with 10 solutions")
                        print(cur_method_results)
                        print()
                if cur_method_results[0] == 'sharpSAT':
                    method4_found = True



            #sharpSAT solved the problem in < 2 seconds
            if len(results_for_each_method) == 1 and results_for_each_method[0][0] == 'sharpSAT' and sat_solver_time < 2:
                    number_of_problems_sharp_sat_solved_fast += 1
            # only get results where the particular random seed finished                  
            # elif len(results_for_each_method) != 8:
            # elif len(results_for_each_method) != 2:
            elif not method1_found or not method2_found or not method3_found or method4_found:
            # only get results where both random seeds finished
            # elif len(results_for_each_method) != 2 or len(results_for_each_method_trial0) != 2 or len(results_for_each_method_trial1) != 2:
                finished_methods = []
                for cur_method_results in results_for_each_method:
                    finished_methods.append(cur_method_results[0])
                # print(problem_name, "unfinished, finished methods are:", finished_methods)
                number_of_unfinished_problems += 1
                # print(len(results_for_each_method), len(results_for_each_method_trial0), len(results_for_each_method_trial1), problem_name)
            else:
                number_of_finished_problems += 1
                fast_problem_for_1point5 = False
                slow_problem_for_3 = False
                large_problem = False
                time_1p5 = None
                time_1 = None
                bound_1p5 = None
                bound_1 = None
                # for cur_method_results in results_for_each_method:
                #     if cur_method_results[0] == 'biregular_variable_degree_1.5' and cur_method_results[3] < 10.0:
                #         fast_problem_for_1point5 = True
                #         time_1p5 = cur_method_results[3]
                #         bound_1p5 = cur_method_results[2]
                #     if cur_method_results[0] == 'biregular_variable_degree_3' and cur_method_results[3] > 10.0:
                #         slow_problem_for_3 = True
                #     if cur_method_results[0] == 'biregular_variable_degree_1' and cur_method_results[2] > 300:
                #         large_problem = True
                #     if cur_method_results[0] == 'biregular_variable_degree_1':
                #         time_1 = cur_method_results[3]
                #         bound_1 = cur_method_results[2]


                # if not large_problem:
                #     pass
                # #     continue
                # elif large_problem and fast_problem_for_1point5:
                #     pass
                #     # print('!'*20, problem_name, "is large and fast for 1.5, runtime 1.5:", time_1p5, "bound 1.5:", bound_1p5, "runtime 1:", time_1, "bound 1:", bound_1)
                # if large_problem and (bound_1p5 is not None) and (bound_1 is not None) and bound_1p5 < bound_1 and bound_1 < 400:
                #     print('@'*20, problem_name, "is large and original 1 is better than 1.5, runtime 1.5:", time_1p5, "bound 1.5:", bound_1p5, "runtime 1:", time_1, "bound 1:", bound_1)
                #     exit(0)
                # # if not fast_problem_for_1point5 or slow_problem_for_3:
                # #     continue

                # # if results_for_each_method[0][2] > 200:
                # #     print(problem_name, 'is fast for 1.5, slow for 3, and has a lower bound > 200')
                
                biregular_1p5 = None
                biregular_3 = None
                found_methods = []
                for cur_method_results in results_for_each_method:
                    method_name = cur_method_results[0]
                    # if method_name in found_methods:
                    #     continue
                    # else:
                    #     found_methods.append(method_name)
                    lower_bound = cur_method_results[2]
                    sat_solver_time = cur_method_results[3]
                    runtimes_for_each_method[method_name].append(sat_solver_time)
                    lower_bounds_for_each_method[method_name].append(lower_bound)
                    timeout = cur_method_results[1]
                    if method_name == 'sharpSAT' and lower_bound > 1000:
                        print('sharpSAT has bound of', lower_bound, "on problem", problem_name)                    
                    if method_name == 'biregular_order_vars_by_marginals_assignmentProblem_variable_degree_1_Tsol_10' and lower_bound > 2000:
                        print('biregular_order_vars_by_marginals_assignmentProblem_variable_degree_1_Tsol_10 has bound of', lower_bound, "on problem", problem_name)                    
                    # if timeout:
                    #     print("method", method_name, "timed out on problem", problem_name)
                    if method_name == 'biregular_variable_degree_1.5':
                        biregular_1p5 = lower_bound
                    elif method_name == 'biregular_variable_degree_3':
                        biregular_3 = lower_bound
                if biregular_1p5 is not None and biregular_3 is not None and biregular_3 > biregular_1p5 + 15 and biregular_3 < 100:
                    print('!'*30, 'big difference', '!'*30)
                    print(problem_name)
                    print(results_for_each_method)
    
    print("sharpSAT solved", number_of_problems_sharp_sat_solved_fast, "problems in less than 2 seconds")
    print("there are", number_of_unfinished_problems, "unfinished problems")
    print("there are", number_of_finished_problems, "finished problems")
    return runtimes_for_each_method, lower_bounds_for_each_method


def plot_results_for_2_methods(runtimes_for_each_method, lower_bounds_for_each_method, method1, method2, plot_folder):
    if method1 == 'sharpSAT':
        method_label1 = 'Exact Model Count'
    elif method1 == 'biregular_variable_degree_1_Tsol_1':
        method_label1 = 'Regular Lower Bound'
    elif method1 == 'biregular_variable_degree_1_Tsol_10':
        method_label1 = 'Regular Lower Bound, K = 10'
    elif method1 == 'biregular_order_vars_by_marginals_assignmentProblem_variable_degree_1_Tsol_10':
        method_label1 = 'Adaptive Regular Lower Bound, K = 10'

    if method2 == 'sharpSAT':
        method_label2 = 'Exact Model Count'
    elif method2 == 'biregular_variable_degree_1_Tsol_1':
        method_label2 = 'Regular Lower Bound'
    elif method2 == 'biregular_variable_degree_1_Tsol_10':
        method_label2 = 'Regular Lower Bound, K = 10'
    elif method2 == 'biregular_order_vars_by_marginals_assignmentProblem_variable_degree_1_Tsol_10':
        method_label2 = 'Adaptive Regular Lower Bound, K = 10'


    if not os.path.exists(plot_folder):
        os.makedirs(plot_folder)

    method1_period_removed = method1.replace('.','point')
    method2_period_removed = method2.replace('.','point')
    print("plot_results_for_2_methods called for", method1, method2)
    ########## RUNTIME PLOT ##########
    print("len(runtimes_for_each_method[method1]):", len(runtimes_for_each_method[method1]))
    print("len(runtimes_for_each_method[method2]):", len(runtimes_for_each_method[method2]))
    print("runtimes_for_each_method[method1][:10]:", runtimes_for_each_method[method1][:10])
    print("runtimes_for_each_method[method2][:10]:", runtimes_for_each_method[method2][:10])
    # print("runtimes_for_each_method", runtimes_for_each_method)
    plt.loglog(runtimes_for_each_method[method1], runtimes_for_each_method[method2], marker='x', linestyle="None")

    #plot y=x for comparison
    max_runtime = 1000
    plt.loglog([0,max_runtime], [0,max_runtime], linestyle='-')

    # matplotlib.rcParams.update({'font.size': 30})
    plt.ylabel('%s runtime' % method2, fontsize=18)
    plt.xlabel('%s runtime' % method1, fontsize=18)
    plt.title('Runtime Comparison', fontsize=22)
    # plt.show()
    plt.savefig('%s%sVS%s_runtime_comparison' % (plot_folder, method_label1, method_label2), bbox_inches = "tight")
    plt.close() 
    ########## RUNTIME log ratio PLOT ##########    
    ratios = [math.log(runtimes_for_each_method[method2][i]/runtimes_for_each_method[method1][i])/math.log(2) for i in range(len(runtimes_for_each_method[method1]))]

    # print("runtimes_for_each_method", runtimes_for_each_method)
    plt.plot(runtimes_for_each_method[method1], ratios, marker='x', linestyle="None")

    median_ratio = np.median(ratios)
    mean_ratio = np.mean(ratios)

    plt.axhline(y=median_ratio, label='median')
    plt.axhline(y=mean_ratio, label='mean')
    #plot y=x for comparison
    plt.legend()

    # matplotlib.rcParams.update({'font.size': 30})
    plt.ylabel('runtime old', fontsize=18)
    plt.xlabel('log_2(new_runtime/old_runtime)', fontsize=18)
    plt.title('Runtime Comparison', fontsize=22)
    # plt.show()
    plt.savefig('%s%sVS%s_runtime_logratio_comparison' % (plot_folder, method_label1, method_label2), bbox_inches = "tight")
    plt.close()     

    ########## PLOT log(LB1) vs log(LB2/LB1) ########## APPEARS IN PAPER
    print("len(lower_bounds_for_each_method[method1]):", len(lower_bounds_for_each_method[method1]))
    print("len(lower_bounds_for_each_method[method2]):", len(lower_bounds_for_each_method[method2]))
    print("lower_bounds_for_each_method[method1][:10]:", lower_bounds_for_each_method[method1][:10])
    print("lower_bounds_for_each_method[method2][:10]:", lower_bounds_for_each_method[method2][:10])

    ratios = [lower_bounds_for_each_method[method2][i]-lower_bounds_for_each_method[method1][i] for i in range(len(lower_bounds_for_each_method[method1]))]

    median_ratio = np.median(ratios)
    mean_ratio = np.mean(ratios)

    matplotlib.rcParams.update({'font.size': 10})

    plt.xscale('symlog')
    plt.yscale('symlog')
    plt.plot(lower_bounds_for_each_method[method1], ratios, marker='x', linestyle="None")
    plt.axhline(y=mean_ratio, c='r', linestyle="--", label='mean =%s' % mean_ratio)
    plt.axhline(y=median_ratio, c='b', linestyle="-", label='median = %s' % median_ratio)
    plt.axhline(0, color='k')

    # plt.axhline(y=0, c='k', label='0')
    #plot y=x for comparison
    plt.legend()
    # matplotlib.rcParams.update({'font.size': 30})
    # plt.ylabel('log(%slower_bound/%slower_bound)' % (method2, method1), fontsize=18)
    plt.ylabel('Log Bound Ratio', fontsize=18)
    plt.xlabel('%s' % method_label1, fontsize=18)


    if method_label1 == 'Regular Lower Bound' and method_label2 == 'Regular Lower Bound, K = 10':
        plt.title('Baseline vs. Variance Reduced LB', fontsize=22)
    elif method_label1 == 'Regular Lower Bound, K = 10' and method_label2 == 'Adaptive Regular Lower Bound, K = 10':
        plt.title('Variance Reduced vs. Adaptive LB', fontsize=22)
    else:
        print("method_label1:", method_label1)
        print("method_label2:", method_label2)
        plt.title('Lower Bound Comparison', fontsize=22)


    # plt.title('Lower Bound Comparison', fontsize=22)
    # plt.show()
    plt.savefig('%s%sVS%s_lower_bound_CLEAR_RATIO_comparison' % (plot_folder, method1_period_removed, method2_period_removed), bbox_inches = "tight")
    plt.close() 

    ########## LOWERBOUND PLOT ##########
    print("len(lower_bounds_for_each_method[method1]):", len(lower_bounds_for_each_method[method1]))
    print("len(lower_bounds_for_each_method[method2]):", len(lower_bounds_for_each_method[method2]))
    print("lower_bounds_for_each_method[method1][:10]:", lower_bounds_for_each_method[method1][:10])
    print("lower_bounds_for_each_method[method2][:10]:", lower_bounds_for_each_method[method2][:10])

    plt.plot(lower_bounds_for_each_method[method1], lower_bounds_for_each_method[method2], marker='x', linestyle="None")

    #plot y=x for comparison
    max_lower_bound = max(max(lower_bounds_for_each_method[method1]), max(lower_bounds_for_each_method[method2]))
    plt.plot([0,max_lower_bound], [0,max_lower_bound], linestyle='-')

    # matplotlib.rcParams.update({'font.size': 30})
    plt.ylabel('%s lower_bound' % method2, fontsize=18)
    plt.xlabel('%s lower_bound' % method1, fontsize=18)
    # plt.show()
    plt.savefig('%s%sVS%s_lower_bound_comparison' % (plot_folder, method_label1, method_label2), bbox_inches = "tight")
    plt.close() 

    ########## LOWERBOUND loglog PLOT ##########
    print("len(lower_bounds_for_each_method[method1]):", len(lower_bounds_for_each_method[method1]))
    print("len(lower_bounds_for_each_method[method2]):", len(lower_bounds_for_each_method[method2]))
    print("lower_bounds_for_each_method[method1][:10]:", lower_bounds_for_each_method[method1][:10])
    print("lower_bounds_for_each_method[method2][:10]:", lower_bounds_for_each_method[method2][:10])

    plt.plot(lower_bounds_for_each_method[method1], lower_bounds_for_each_method[method2], marker='x', linestyle="None")

    #plot y=x for comparison
    max_lower_bound = max(max(lower_bounds_for_each_method[method1]), max(lower_bounds_for_each_method[method2]))
    plt.loglog([0,max_lower_bound], [0,max_lower_bound], linestyle='-')

    # matplotlib.rcParams.update({'font.size': 30})
    plt.ylabel('%s lower_bound' % method2, fontsize=18)
    plt.xlabel('%s lower_bound' % method1, fontsize=18)
    plt.title('Lower Bound Comparison', fontsize=22)
    # plt.show()
    plt.savefig('%s%sVS%s_lower_bound_loglog_comparison' % (plot_folder, method_label1, method_label2), bbox_inches = "tight")
    plt.close()     


def plot_bounds_multiple_methods_lb(lower_bounds_for_each_method, x_axis_method, y_axis_methods, plot_folder):
    PLOT_LOG_DIFF = True
    if not os.path.exists(plot_folder):
        os.makedirs(plot_folder)

    method1_period_removed = x_axis_method.replace('.','point')

    ########## LOWERBOUND PLOT ##########
    matplotlib.rcParams.update({'font.size': 10})

    for method_y in y_axis_methods:
        if method_y == 'sharpSAT':
            method_label = 'Exact Model Count'
        elif method_y == 'biregular_variable_degree_1_Tsol_10':
            method_label = 'Regular Lower Bound, K = 10'
        elif method_y == 'biregular_order_vars_by_marginals_assignmentProblem_variable_degree_1_Tsol_10':
            method_label = 'Adaptive Regular Lower Bound, K = 10'

        if PLOT_LOG_DIFF:
            log_difference = [lower_bounds_for_each_method[method_y][i] - lower_bounds_for_each_method[x_axis_method][i] for i in range(len(lower_bounds_for_each_method[x_axis_method]))]
            plt.plot(lower_bounds_for_each_method[x_axis_method], log_difference, marker='x', linestyle="None", label=method_label)
        else:
            plt.plot(lower_bounds_for_each_method[x_axis_method], lower_bounds_for_each_method[method_y], marker='x', linestyle="None", label=method_label)

    plt.legend()

    #plot y=x for comparison
    # max_lower_bound = max(max(lower_bounds_for_each_method[x_axis_method]), max(lower_bounds_for_each_method[method_y]))
    # plt.semilogx([0,max_lower_bound], [0,max_lower_bound], linestyle='-')

    # matplotlib.rcParams.update({'font.size': 30})
    if PLOT_LOG_DIFF:
        plt.ylabel('Log Model Count Ratio', fontsize=18)

    else:
        plt.ylabel('Model Count', fontsize=18)
    plt.xlabel('Regular Lower Bound', fontsize=18)
    # plt.title('Exact Model Count Comparison', fontsize=22)
    plt.title('SharpSAT Timeouts', fontsize=22)
    # plt.title('Lower Bound Comparison', fontsize=22)
    # plt.show()
    plt.savefig('%s%s_x_axis_lower_bound_comparison_noSharpSAT_log' % (plot_folder, method1_period_removed), bbox_inches = "tight")
    # plt.savefig('%s%s_x_axis_lower_bound_comparison_log' % (plot_folder, method1_period_removed), bbox_inches = "tight")
    plt.close() 

 


def plot_method_1_better_than_2(lower_bounds_for_each_method0, lower_bounds_for_each_method1, method1, method2, min_bound, max_bound, plot_folder='./plots_fixSharpSat/'):
    method1_period_removed = method1.replace('.','point')
    method2_period_removed = method2.replace('.','point')
    print("plot_results_for_2_methods called for", method1, method2)

    method1_bounds_run0 = lower_bounds_for_each_method0[method1]
    method2_bounds_run0 = lower_bounds_for_each_method0[method2]
    method1_bounds_run1 = lower_bounds_for_each_method1[method1]
    method2_bounds_run1 = lower_bounds_for_each_method1[method2]


    run0_method1_bounds = []
    run0_method2_bounds = []
    run1_method1_bounds = []
    run1_method2_bounds = []

    assert(len(method1_bounds_run0) == len(method2_bounds_run0)), (len(method1_bounds_run0), len(method2_bounds_run0))
    assert(len(method1_bounds_run1) == len(method2_bounds_run1)), (len(method1_bounds_run1), len(method2_bounds_run1))
    assert(len(method1_bounds_run0) == len(method1_bounds_run1)), (len(method1_bounds_run0), len(method1_bounds_run1))

    for instance_idx, bound in enumerate(method1_bounds_run0):
        if bound > min_bound and bound < max_bound and bound > method2_bounds_run0[instance_idx]:
            print('found a bound!!!')
            run0_method1_bounds.append(bound)
            run0_method2_bounds.append(method2_bounds_run0[instance_idx])

            run1_method1_bounds.append(method1_bounds_run1[instance_idx])
            run1_method2_bounds.append(method2_bounds_run1[instance_idx])

    print("found a total of", len(run0_method1_bounds), "bounds")
    ########## LOWERBOUND PLOT ##########
    plt.plot(run0_method1_bounds, run0_method2_bounds, marker='x', linestyle="None")

    #plot y=x for comparison
    max_lower_bound = max((max(run0_method1_bounds), max(run0_method2_bounds)))
    plt.plot([0,max_lower_bound], [0,max_lower_bound], linestyle='-')

    # matplotlib.rcParams.update({'font.size': 30})
    plt.ylabel('%s lower_bound' % method2, fontsize=18)
    plt.xlabel('%s lower_bound' % method1, fontsize=18)
    plt.title('Lower Bound Comparison', fontsize=22)
    # plt.show()
    plt.savefig('%sRUN0%sVS%s_lower_bound_comparison' % (plot_folder, method1_period_removed, method2_period_removed), bbox_inches = "tight")
    plt.close() 

    ########## LOWERBOUND PLOT ##########
    plt.plot(run1_method1_bounds, run1_method2_bounds, marker='x', linestyle="None")

    #plot y=x for comparison
    max_lower_bound = max((max(run1_method1_bounds), max(run1_method2_bounds)))
    plt.plot([0,max_lower_bound], [0,max_lower_bound], linestyle='-')

    # matplotlib.rcParams.update({'font.size': 30})
    plt.ylabel('%s lower_bound' % method2, fontsize=18)
    plt.xlabel('%s lower_bound' % method1, fontsize=18)
    plt.title('Lower Bound Comparison', fontsize=22)
    # plt.show()
    plt.savefig('%sRUN1%sVS%s_lower_bound_comparison' % (plot_folder, method1_period_removed, method2_period_removed), bbox_inches = "tight")
    plt.close()   


def plot_computation_tradeoff(plot_folder):
    baseline_runtimes = [1.976, 10.9, 63.31999999999999, 332.072, 1053.824, 2502.3]
    baseline_bounds = [1861, 2062, 2277, 2320, 2525, 2962]

    matplotlib.rcParams.update({'font.size': 10})
    plt.semilogx(baseline_runtimes, baseline_bounds, marker='x', linestyle="None", label='Baseline Regular Bounds', markersize=10)

    plt.semilogx([5.492000000000026], [2673], marker='*', linestyle="None", label='Adaptive Regular Bound, Naive Parallel', markersize=10)
    plt.semilogx([45.888000000000005], [2673], marker='*', linestyle="None", label='Adaptive Regular Bound, Sequential', markersize=10)
    plt.legend(loc='lower right')

    # plt.ylim(bottom=0)
    # matplotlib.rcParams.update({'font.size': 30})
    plt.ylabel('Log Lower Bound', fontsize=18)
    plt.xlabel('Runtime (seconds)', fontsize=18)
    plt.title('Runtime vs. Bound Tradeoff', fontsize=22)
    # plt.show()
    plt.savefig('%s/runtime_bound_tradeoff90-34-3-q' % plot_folder , bbox_inches = "tight")
    plt.close()


if __name__=="__main__":
    runtimes_for_each_method_0, lower_bounds_for_each_method_0 = get_results('/atlas/u/jkuck/F2/fireworks/cluster_results_orderVarsByMarginals_randomInChunks_Tsol10_parallelRuntime', random_seed=0)
    # runtimes_for_each_method_1, lower_bounds_for_each_method_1 = get_results('/atlas/u/jkuck/F2/fireworks/cluster_results_orderVarsByMarginals_randomInChunks_Tsol10_parallelRuntime', random_seed=1)
    
    plot_folder='./plots_orderVarsByMarginals_randomInChunks_Tsol10_parallelRuntime_paper/'
    # plot_computation_tradeoff(plot_folder)
    # exit(0)
    # plot_method_1_better_than_2(lower_bounds_for_each_method_0, lower_bounds_for_each_method_1,\
    #                             method1='biregular_variable_degree_1',\
    #                             method2='biregular_order_vars_by_marginals_assignmentProblem_variable_degree_1', \
    #                             min_bound=0, max_bound=10000, plot_folder=plot_folder)
    # exit(0)

    # biregular_order_vars_by_marginals_assignmentProblem_variable_degree_1_Tsol_10
    # biregular_variable_degree_1_Tsol_10
    # biregular_variable_degree_1_Tsol_1


    # these next 3 are the ones for the paper !!
    # plot_results_for_2_methods(runtimes_for_each_method_0, lower_bounds_for_each_method_0, method1='biregular_variable_degree_1_Tsol_1', method2='biregular_variable_degree_1_Tsol_10', plot_folder=plot_folder)
    # plot_results_for_2_methods(runtimes_for_each_method_0, lower_bounds_for_each_method_0, method1='biregular_variable_degree_1_Tsol_10', method2='biregular_order_vars_by_marginals_assignmentProblem_variable_degree_1_Tsol_10', plot_folder=plot_folder)
    # plot_results_for_2_methods(runtimes_for_each_method_0, lower_bounds_for_each_method_0, method1='biregular_variable_degree_1_Tsol_1', method2='biregular_order_vars_by_marginals_assignmentProblem_variable_degree_1_Tsol_10', plot_folder=plot_folder)
    

    # plot_results_for_2_methods(runtimes_for_each_method_0, lower_bounds_for_each_method_0, method1='sharpSAT', method2='biregular_variable_degree_1_Tsol_1', plot_folder=plot_folder)
    # plot_results_for_2_methods(runtimes_for_each_method_0, lower_bounds_for_each_method_0, method1='sharpSAT', method2='biregular_variable_degree_1_Tsol_10', plot_folder=plot_folder)
    # plot_results_for_2_methods(runtimes_for_each_method_0, lower_bounds_for_each_method_0, method1='sharpSAT', method2='biregular_order_vars_by_marginals_assignmentProblem_variable_degree_1_Tsol_10', plot_folder=plot_folder)
    
    # y_axis_methods = ['biregular_variable_degree_1_Tsol_1', 'biregular_variable_degree_1_Tsol_10', 'biregular_order_vars_by_marginals_assignmentProblem_variable_degree_1_Tsol_10']
    # y_axis_methods = ['sharpSAT', 'biregular_variable_degree_1_Tsol_10', 'biregular_order_vars_by_marginals_assignmentProblem_variable_degree_1_Tsol_10']
    y_axis_methods = ['biregular_variable_degree_1_Tsol_10', 'biregular_order_vars_by_marginals_assignmentProblem_variable_degree_1_Tsol_10']
    plot_bounds_multiple_methods_lb(lower_bounds_for_each_method_0, x_axis_method='biregular_variable_degree_1_Tsol_1', y_axis_methods=y_axis_methods, plot_folder=plot_folder)
 
    # plot_results_for_2_methods(runtimes_for_each_method_0, lower_bounds_for_each_method_0, method1='biregular_variable_degree_1', method2='biregular_order_vars_by_marginals_assignmentProblem_variable_degree_1', plot_folder=plot_folder)
    
    # plot_results_for_2_methods(runtimes_for_each_method_0, lower_bounds_for_each_method_0, method1='sharpSAT', method2='biregular_variable_degree_3', plot_folder=plot_folder)
    # plot_results_for_2_methods(runtimes_for_each_method_0, lower_bounds_for_each_method_0, method1='sharpSAT', method2='long_iid_.5', plot_folder=plot_folder)
    # plot_results_for_2_methods(runtimes_for_each_method_0, lower_bounds_for_each_method_0, method1='long_iid_.5', method2='biregular_variable_degree_3', plot_folder=plot_folder)

    # plot_results_for_2_methods(runtimes_for_each_method_0, lower_bounds_for_each_method_0, method1='biregular_variable_degree_3', method2='biregular_variable_degree_1.5', plot_folder=plot_folder)
    # plot_results_for_2_methods(runtimes_for_each_method_0, lower_bounds_for_each_method_0, method1='biregular_variable_degree_1.5', method2='bi_regular_marginals_joint_constraint_variable_degree_1.5', plot_folder=plot_folder)
    # plot_results_for_2_methods(runtimes_for_each_method_0, lower_bounds_for_each_method_0, method1='biregular_variable_degree_1.5', method2='bi_regular_marginals_per_constraint_iid_variable_degree_1.5', plot_folder=plot_folder)
    # plot_results_for_2_methods(runtimes_for_each_method_0, lower_bounds_for_each_method_0, method1='biregular_variable_degree_1.5', method2='bi_regular_marginals_per_constraint_keep_biregular_variable_degree_1.5', plot_folder=plot_folder)
    # plot_results_for_2_methods(runtimes_for_each_method_0, lower_bounds_for_each_method_0, method1='biregular_variable_degree_3', method2='biregular_order_vars_by_marginals_variable_degree_1.5', plot_folder=plot_folder)

    # plot_results_for_2_methods(runtimes_for_each_method_0, lower_bounds_for_each_method_0, method1='bi_regular_marginals_joint_constraint_variable_degree_1.5', method2='bi_regular_marginals_per_constraint_iid_variable_degree_1.5', plot_folder=plot_folder)
    # plot_results_for_2_methods(runtimes_for_each_method_0, lower_bounds_for_each_method_0, method1='bi_regular_marginals_joint_constraint_variable_degree_1.5', method2='bi_regular_marginals_per_constraint_keep_biregular_variable_degree_1.5', plot_folder=plot_folder)

    # plot_results_for_2_methods(runtimes_for_each_method_0, lower_bounds_for_each_method_0, method1='biregular_variable_degree_1.5', method2='biregular_order_vars_by_marginals_assignmentProblem_variable_degree_1.5', plot_folder=plot_folder)
    # plot_results_for_2_methods(runtimes_for_each_method_0, lower_bounds_for_each_method_0, method1='biregular_variable_degree_1.5', method2='biregular_order_vars_by_marginals_assignmentProblem_variable_degree_1', plot_folder=plot_folder)
    # plot_results_for_2_methods(runtimes_for_each_method_0, lower_bounds_for_each_method_0, method1='biregular_variable_degree_1.5', method2='biregular_variable_degree_1', plot_folder=plot_folder)
    # plot_results_for_2_methods(runtimes_for_each_method_0, lower_bounds_for_each_method_0, method1='biregular_variable_degree_1', method2='biregular_order_vars_by_marginals_assignmentProblem_variable_degree_1', plot_folder=plot_folder)
    # plot_results_for_2_methods(runtimes_for_each_method_0, lower_bounds_for_each_method_0, method1='biregular_variable_degree_1', method2='biregular_order_vars_by_doubleMarginals_assignmentProblem_variable_degree_1', plot_folder=plot_folder)
    # plot_results_for_2_methods(runtimes_for_each_method_0, lower_bounds_for_each_method_0, method1='biregular_variable_degree_1.5', method2='biregular_order_vars_by_doubleMarginals_assignmentProblem_variable_degree_1', plot_folder=plot_folder)
    # plot_results_for_2_methods(runtimes_for_each_method_0, lower_bounds_for_each_method_0, method1='biregular_order_vars_by_marginals_assignmentProblem_variable_degree_1.5', method2='biregular_order_vars_by_marginals_assignmentProblem_variable_degree_1', plot_folder=plot_folder)

# biregular_order_vars_by_marginals_variable_degree_1
# biregular_variable_degree_1

# biregular_order_vars_by_marginals_variable_degree_1.5



