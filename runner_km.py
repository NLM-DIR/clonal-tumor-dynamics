import pandas as pd
import numpy as np
from estimator import Estimator
import argparse
import time
import json

##########
# This script estimates the interaction terms k and m for the experiments with
# mixtures of the proliferative and invasive subclones in nude mice.
# It uses a configuration file that includes the details for the estimations, 
# which includes the file with subclonal growth rates to use.
# # User input defines the range and increment of values to try.
# The result file is saved to the user defined folder with the file name
# 'k_{kmin}_{kmax}_{kinc}_m_{mmin}_{mmax}_{minc}.csv'
##########


def main(config, m_min, m_max, m_inc, k_min, k_max, k_inc, save_path):
    # Create an estimator with the configuration file
    es = Estimator(config)

    m_precision = max(str(m_min)[::-1].find("."), str(m_max)[::-1].find("."), str(m_inc)[::-1].find("."))
    k_precision = max(str(k_min)[::-1].find("."), str(k_max)[::-1].find("."), str(k_inc)[::-1].find("."))

    m_list = [round(m_min+m_inc*i, m_precision) for i in range(round((m_max-m_min)/m_inc+1))]
    k_list = [round(k_min+k_inc*i, k_precision) for i in range(round((k_max-k_min)/k_inc+1))]

    # Print information
    print("m list:", m_list)
    print("k list:", k_list)
    print("Number of grid elements:", len(m_list)*len(k_list), flush=True)

    with open(config, "r") as f:
            config_params = json.load(f)
    
    # Read in growth rates
    gs = pd.read_csv(config_params["g_result_file"])
    # To lists
    g1s = gs[["g1"]].to_numpy()
    g11s = gs[["g11"]].to_numpy()
    g1s = g1s[g1s != 0]
    g11s = g11s[g11s != 0]

    print("Number of g1s:", len(g1s))
    print("Number of g11s:", len(g11s))
    print("Total grid size:", len(g1s)*len(g11s)*len(m_list)*len(k_list), flush=True)
    
    # Get the k and m estimates for all admixtures
    start_time = time.time()
    results = es.nude_admix_all(g1s, g11s, m_list, k_list)
    print("Time:", time.time() - start_time)

    # Drop results that have the wrong subline winning
    results = results[results["error"] != np.inf]
    print(results)
    
    # Save results
    print("Saving to {}m_{}_{}_{}_k_{}_{}_{}.csv".format(save_path, m_min, m_max, m_inc, k_min, k_max, k_inc), flush=True)
    results.to_csv("{}m_{}_{}_{}_k_{}_{}_{}.csv".format(save_path, m_min, m_max, m_inc, k_min, k_max, k_inc), index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", type=str, action="store",
                        help="The configuration file that includes the details for the estimation." \
                        "See the readme for details.")
    parser.add_argument("--mmin", type=float, action="store",
                        help="The minimum value of m to evaluate.")
    parser.add_argument("--mmax", type=float, action="store",
                        help="The maximum value of k to evaluate.")
    parser.add_argument("--minc", type=float, action="store",
                        help="The amount to increment between m values to evaluate.")
    parser.add_argument("--kmin", type=float, action="store",
                        help="The minimum value of k to evaluate.")
    parser.add_argument("--kmax", type=float, action="store",
                        help="The maximum value of k to evaluate.")
    parser.add_argument("--kinc", type=float, action="store",
                        help="The amount to increment between k values to evaluate.")
    parser.add_argument("-s", "--save_path", type=str, action="store",
                        help="The path in which to save the result file.")
    args = parser.parse_args()

    # If we only have one value we want to run, we need to make the bounds to be able to make an array with one value
    if args.mmin == args.mmax or args.minc == 0: 
        args.mmax = args.mmin
        args.minc = 0.1
    if args.kmin == args.kmax or args.kinc == 0: 
        args.kmax = args.kmin
        args.kinc = 0.1


    main(args.config, args.mmin, args.mmax, args.minc, args.kmin, args.kmax, args.kinc, args.save_path)