import pandas as pd
import numpy as np
from estimator import Estimator
import argparse
import time
import json

##########
# This script estimates the terms c, a, f, b, l, and h for the experiments with
# mixtures of the proliferative and invasive subclones in B6 mice including the interaction
# between the invasive subclone and T cells.
# It uses a configuration file that includes the details for the estimations,
# including the file with the tuples of g1, g11, k, and m values to use.
# User input defines the range and increment of values to try.
# The result file is saved to the user defined folder with the file name
# 'c_{cmin}_{cmax}_{cinc}_a_{amin}_{amax}_{ainc}_f_{fmin}_{fmax}_{finc}_b_{bmin}_{bmax}_{binc}_l_{lmin}_{lmax}_{linc}_h_{hmin}_{hmax}_{hinc}.csv'
##########


def main(config, c_min, c_max, c_inc, a_min, a_max, a_inc, f_min, f_max, f_inc, b_min, b_max, b_inc, l_min, l_max, l_inc, h_min, h_max, h_inc, save_path):
    # Create an estimator with the configuration file
    es = Estimator(config)

    c_precision = max(str(c_min)[::-1].find("."), str(c_max)[::-1].find("."), str(c_inc)[::-1].find("."))
    a_precision = max(str(a_min)[::-1].find("."), str(a_max)[::-1].find("."), str(a_inc)[::-1].find("."))
    f_precision = max(str(f_min)[::-1].find("."), str(f_max)[::-1].find("."), str(f_inc)[::-1].find("."))
    b_precision = max(str(b_min)[::-1].find("."), str(b_max)[::-1].find("."), str(b_inc)[::-1].find("."))
    l_precision = max(str(l_min)[::-1].find("."), str(l_max)[::-1].find("."), str(l_inc)[::-1].find("."))
    h_precision = max(str(h_min)[::-1].find("."), str(h_max)[::-1].find("."), str(h_inc)[::-1].find("."))
    print("Precision: d=", c_precision, "a=", a_precision, "f=", f_precision, "b=", b_precision, "l=", l_precision,  "h=", h_precision)

    c_list = [round(c_min+c_inc*i, c_precision) for i in range(round((c_max-c_min)/c_inc))]
    a_list = [round(a_min+a_inc*i, a_precision) for i in range(round((a_max-a_min)/a_inc))]
    f_list = [round(f_min+f_inc*i, f_precision) for i in range(round((f_max-f_min)/f_inc))]
    b_list = [round(b_min+b_inc*i, b_precision) for i in range(round((b_max-b_min)/b_inc))]
    l_list = [round(l_min+l_inc*i, l_precision) for i in range(round((l_max-l_min)/l_inc))]
    h_list = [round(h_min+h_inc*i, h_precision) for i in range(round((h_max-h_min)/h_inc))]
    
    if len(c_list)==0: c_list = [c_min]
    if len(a_list)==0: a_list = [a_min]
    if len(f_list)==0: f_list = [f_min]
    if len(b_list)==0: b_list = [b_min]
    if len(l_list)==0: l_list = [l_min]
    if len(h_list)==0: h_list = [h_min]

    # Print information
    print("c list:", c_list)
    print("a list:", a_list)
    print("f list:", f_list)
    print("b list:", b_list)
    print("l list:", l_list)
    print("h list:", h_list)

    print("Number of grid elements:", len(c_list)*len(a_list)*len(f_list)*len(b_list)*len(l_list)*len(h_list), flush=True)

    with open(config, "r") as f:
        config_params = json.load(f)
    
    # Read in g, k and m pairs
    gkm_pairs = pd.read_csv(config_params["km_result_file"])
    gkm_pairs = gkm_pairs[["g1", "g11", "m", "k"]].to_numpy()
    
    print("Number of k and m pairs:", len(gkm_pairs))
    print("Total grid size:", len(gkm_pairs)*len(c_list)*len(a_list)*len(f_list)*len(b_list)*len(l_list), flush=True)
    
    # Run simulations and get results
    start_time = time.time()
    results = es.b6_all_it(gkm_pairs, c_list, a_list, f_list, b_list, l_list, h_list)
    
    print("Time:", time.time() - start_time)
    results = results[results["error"] != np.inf]
    results = results.astype({"id": "uint16", "g1":"float32", "g11":"float32", "m":"float32", "k":"float32", 
                              "c":"float32", "a":"float32", "f":"float32", "b":"float32", "l":"float32", "h":"float32",
                              "error":"float32", "winner":"uint8"})

    # Save results
    results.to_csv("{}c_{}_{}_{}_a_{}_{}_{}_f_{}_{}_{}_b_{}_{}_{}_l_{}_{}_{}_h_{}_{}_{}.csv".format(
         save_path, c_min, c_max, c_inc, a_min, a_max, a_inc,
         f_min, f_max, f_inc, b_min, b_max, b_inc, l_min, l_max, l_inc, h_min, h_max, h_inc), index=False)
    print(results)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", type=str, action="store",
                        help="The configuration file that includes the details for the estimation." \
                        "See the readme for details.")
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
    parser.add_argument("--bmin", type=float, action="store",
                        help="The minimum value of b to evaluate.")
    parser.add_argument("--bmax", type=float, action="store",
                        help="The maximum value of b to evaluate.")
    parser.add_argument("--binc", type=float, action="store",
                        help="The amount to increment between b values to evaluate.")
    parser.add_argument("--lmin", type=float, action="store",
                        help="The minimum value of l to evaluate.")
    parser.add_argument("--lmax", type=float, action="store",
                        help="The maximum value of l to evaluate.")
    parser.add_argument("--linc", type=float, action="store",
                        help="The amount to increment between l values to evaluate.")
    parser.add_argument("--hmin", type=float, action="store",
                        help="The minimum value of h to evaluate.")
    parser.add_argument("--hmax", type=float, action="store",
                        help="The maximum value of h to evaluate.")
    parser.add_argument("--hinc", type=float, action="store",
                        help="The amount to increment between h values to evaluate.")
    parser.add_argument("-s", "--save_path", type=str, action="store",
                        help="The path in which to save the result file.")
    args = parser.parse_args()

    main(args.config, args.cmin, args.cmax, args.cinc, args.amin, args.amax, args.ainc, 
         args.fmin, args.fmax, args.finc, args.bmin, args.bmax, args.binc, args.lmin, args.lmax, args.linc, 
         args.hmin, args.hmax, args.hinc, args.save_path)