import pandas as pd
import numpy as np
from estimator import Estimator
import argparse

##########
# Helper function to get the ratio of the sublines given the group name
# @param group Group name
# @return List containing (in order), the fraction of C1, the fraction of C11, and the initial T cell population as dictated by the config file.
##########
def get_init(group):
    if "A1" in group or "B1" in group: return [1, 0]
    if "A2" in group or "B2" in group: return [0.8, 0.2]
    if "A3" in group or "B3" in group: return [0.5, 0.5]
    if "A4" in group or "B4" in group: return [0.2, 0.8]
    if "A5" in group or "B5" in group: return [0, 1]

##########
# Helper function to run the model.
# @param es Estimator object
# @param init List of initial values of C1, C11, and T cells
# @param g1 Growth rate of C1
# @param g11 Growth rate of C11
# @param k Impact of C11 on C1
# @param m Impact of C1 on C11
# @param d Impact of T cells on C1
# @param a Impact of C1 on T cells
# @param f Exhaustion/death rate of T cells
# @param l Michaelis–Menten constant (limiting factor on T cell recruitment rate)
# @return The solution, the time to switch from equation 1 to equation 2, the time to switch from equation 2 to equation 3.
#         The solution is in the form of a 4-element list.
#         sol[0] is a list of the sizes of C1 at each time point.
#         sol[1] is a list of the sizes of C11 at each time point.
#         sol[2] is a list of the sizes of the T cell population at each time point.
#         sol[3] is a list of the time points.
def run_model(es, init, g1, g11, k, m, d, a, f, l, end_time):
    sol = es.run_before_t(init, g1, g11, k, m, d, end_time)
    t1 = sol[3][-1]
    sol = es.run_incr_t(sol, g1, g11, k, m, d, a, f, l, end_time)
    t2 = sol[3][-1]
    sol = es.run_dec_t(sol, g1, g11, k, m, d, a, f, l, end_time) 
    return sol, t1, t2

##########
# Function to generate curves of all parameters included in df and find the minimum 
# and maximum values at each time point across all generated curves
# @param es Estimator object
# @param group Mouse group name (e.g. "A2" or "Group A2 B6 (80% C1; 20% C11)")
# @param time_range The time range to generate the curves. List of two elements. time_range[0] is start time, time_range[1] is end time.
# @param t_init Initial value of T cell population
# @param g1 Growth rate of C1
# @param g11 Growth rate of C11
# @param df The dataframe containing all parameter combinations to evaluate
# @return List of minimum sizes of C1, list of maximum sizes of C1, 
#         list of minimum sizes of C11, list of maximum sizes of C11,
#         list of time points
##########
def find_min_max(es, group, time_range, t_init, g1, g11, df):
    init = get_init(group)
    # Add T cell initial size to the inital value list
    init += [t_init]

    min_sol, _, _ = run_model(es, init, g1, g11, df.loc[0,"k"], df.loc[0,"m"], df.loc[0,"d"], df.loc[0,"a"], df.loc[0,"f"], df.loc[0,"l"], time_range[1])     
    min_c1 = min_sol[0]
    max_c1 = min_c1.copy()
    min_c11 = min_sol[1]
    max_c11 = min_c11.copy()

    for idx, row in df.iterrows():
        sol, _, _ = run_model(es, init, g1, g11, row["k"], row["m"], row["d"], row["a"], row["f"], row["l"], time_range[1])     
        
        min_c1 = np.amin([min_c1, sol[0]], axis=0)
        max_c1 = np.amax([max_c1, sol[0]], axis=0)
        min_c11 = np.amin([min_c11, sol[1]], axis=0)
        max_c11 = np.amax([max_c11, sol[1]], axis=0)

    return min_c1, max_c1, min_c11, max_c11, sol[3]

##########
# Function to find the distribution of extinction times of all curves generated from a set of parameter combinations
# @param es Estimator object
# @param group Mouse group name (e.g. "A2" or "Group A2 B6 (80% C1; 20% C11)")
# @param time_range The time range to generate the curves. List of two elements. time_range[0] is start time, time_range[1] is end time.
# @param t_init Initial value of T cell population
# @param g1 Growth rate of C1
# @param g11 Growth rate of C11
# @param df The dataframe containing all parameter combinations to evaluate
# @return List of minimum sizes of C1, list of maximum sizes of C1, 
#         list of minimum sizes of C11, list of maximum sizes of C11,
#         list of time points
##########
def min_max_dist(es, group, time_range, t_init, g1, g11, df):
    if "A2" in group: c="c11" # C11 goes extinct
    else: c="c1" # C1 goes extinct

    init = get_init(group)
    # Add T cell initial size to the inital value list
    init += [t_init]

    # Run the model
    extinct_times = []

    for idx, row in df.iterrows():
        sol, _, _ = run_model(es, init, g1, g11, row["k"], row["m"], row["d"], row["a"], row["f"], row["l"], time_range[1])     

        if c == "c1":
            extinct_times += [sol[3][np.where(sol[0]==0)[0][0]]]
        else:
            extinct_times += [sol[3][np.where(sol[1]==0)[0][0]]]

    return np.asarray(extinct_times)

def main(r_group, m_group):
    config = "config_tcell.json"
    results_file = "result_avg_b6.csv"
    
    es = Estimator(config)
    
    results = pd.read_csv(results_file)

    g1 = 0.135  # Growth rate of C1 (proliferative subclone)
    g11 = 0.107 # Growth rate of C11 (invasive subclone)
    t_init = 0.001 # Initial T cell population size; 0.001
    time_range = [0,125] # Time range to run the simulations
    
    # Subset the parameter combinations to those in the group we want to evaluate
    group_results = results[results["group"] == r_group].reset_index()
    
    # Find minimum and maximum values across all generated curves in that group
    min_c1, max_c1, min_c11, max_c11, t = find_min_max(es, m_group, time_range, t_init, g1, g11, group_results)
    min_c1.tofile("../min_max_curve_results/{}_{}_min_c1.csv".format(r_group, m_group), sep=",")
    max_c1.tofile("../min_max_curve_results/{}_{}_max_c1.csv".format(r_group, m_group), sep=",")
    min_c11.tofile("../min_max_curve_results/{}_{}_min_c11.csv".format(r_group, m_group), sep=",")
    max_c11.tofile("../min_max_curve_results/{}_{}_max_c11.csv".format(r_group, m_group), sep=",")

    # Find extinction time distribution for that group
    extinct_times = min_max_dist(es, m_group, time_range, t_init, g1, g11, group_results)
    extinct_times.tofile("../min_max_dist/{}_{}.csv".format(r_group, m_group), sep=",")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", type=int, action="store", description="Number of the group (1a = 1, 1b = 2, 2a = 3, 2b = 4")
    parser.add_argument("-m", type=str, action="store", description="Mouse group ID (A1, A2, A3, A4)")
    args = parser.parse_args()
    main(args.r, args.m)