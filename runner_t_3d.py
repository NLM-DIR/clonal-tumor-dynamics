import pandas as pd
import numpy as np
from estimator import Estimator
import argparse
import time
import json

##########
# This script estimates the results used to make the 3d figures in Figure 6.
# For the optimal values of g1, g11, k, m, and b for each mouse provided in the bests_file, 
# determine the outcome of the simulation for a range of c, a, and f values.
# The result file is saved to the user defined folder with the file name
# "c_{cmin}_{cmax}_{cinc}_a_{amin}_{amax}_{ainc}_f_{fmin}_{fmax}_{finc}.csv"
##########


def main(config, bests_file, c_min, c_max, c_inc, a_min, a_max, a_inc, f_min, f_max, f_inc, save_path):
# def main(config, bests_file):
    # Create an estimator with the configuration file
    es = Estimator(config)

    c_precision = max(str(c_min)[::-1].find("."), str(c_max)[::-1].find("."), str(c_inc)[::-1].find("."))
    a_precision = max(str(a_min)[::-1].find("."), str(a_max)[::-1].find("."), str(a_inc)[::-1].find("."))
    f_precision = max(str(f_min)[::-1].find("."), str(f_max)[::-1].find("."), str(f_inc)[::-1].find("."))
    print("Precision: c=", c_precision, "a=", a_precision, "f=", f_precision)

    c_list = [round(c_min+c_inc*i, c_precision) for i in range(round((c_max-c_min)/c_inc+1))]
    a_list = [round(a_min+a_inc*i, a_precision) for i in range(round((a_max-a_min)/a_inc+1))]
    f_list = [round(f_min+f_inc*i, f_precision) for i in range(round((f_max-f_min)/f_inc+1))]
    
    # Print information
    print("c list:", c_list)
    print("a list:", a_list)
    print("f list:", f_list)

    print("Number of grid elements:", len(c_list)*len(a_list)*len(f_list), flush=True)

    # Read in g, k and m pairs
    gkm_pairs_df = pd.read_csv(bests_file)
    gkm_pairs = gkm_pairs_df[["g1", "g11", "m", "k"]].to_numpy()
    
    print("Number of k and m pairs:", len(gkm_pairs))
    print("Total grid size:", len(gkm_pairs)*len(c_list)*len(a_list)*len(f_list), flush=True)

    # Run simulations and get results
    start_time1 = time.time()
    results = []
    for idx, row in gkm_pairs_df.iterrows():
        start_time = time.time()
        gkm_pairs = row[["g1", "g11", "m", "k"]].to_numpy()

        print(idx)
        curr = es.b6_group(row["group"], [gkm_pairs], c_list, a_list, f_list, [row["b"]], [0.01])
        results += [curr]
    
        print("Time:", time.time() - start_time)
    print("Time:", time.time() - start_time1)
    results = pd.concat(results)
    print(results)

    results.to_csv("{}c_{}_{}_{}_a_{}_{}_{}_f_{}_{}_{}.csv".format(
         save_path, c_min, c_max, c_inc, a_min, a_max, a_inc,
         f_min, f_max, f_inc), index=False)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", type=str, action="store",
                        help="The configuration file that includes the details for the estimation." \
                        "See the readme for details.")
    parser.add_argument("-b", "--bests", type=str, action="store",
                        help="The file containing the optimal parameter values for each mouse.")
    parser.add_argument("--cmin", type=float, action="store",
                        help="The minimum value of c to evaluate.")
    parser.add_argument("--cmax", type=float, action="store",
                        help="The maximum value of c to evaluate.")
    parser.add_argument("--cinc", type=float, action="store",
                        help="The amount to increment between c values to evaluate.")
    parser.add_argument("--amin", type=float, action="store",
                        help="The minimum value of a to evaluate.")
    parser.add_argument("--amax", type=float, action="store",
                        help="The maximum value of a to evaluate.")
    parser.add_argument("--ainc", type=float, action="store",
                        help="The amount to increment between a values to evaluate.")
    parser.add_argument("--fmin", type=float, action="store",
                        help="The minimum value of f to evaluate.")
    parser.add_argument("--fmax", type=float, action="store",
                        help="The maximum value of f to evaluate.")
    parser.add_argument("--finc", type=float, action="store",
                        help="The amount to increment between f values to evaluate.")
    parser.add_argument("--eq", type=str, action="store",
                        help="The equation set to use.")
    parser.add_argument("-s", "--save_path", type=str, action="store",
                        help="The path in which to save the result file.")
    args = parser.parse_args()

    main(args.config, args.bests, args.cmin, args.cmax, args.cinc, args.amin, args.amax, args.ainc, 
         args.fmin, args.fmax, args.finc, args.save_path)