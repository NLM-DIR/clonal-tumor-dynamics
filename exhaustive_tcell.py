import pandas as pd
import numpy as np
from estimator import Estimator
import argparse
import time

##########
# This script estimates m, d, a, f, and l for the experiments with
# mixtures of the proliferative and invasive subclones in B6 mice.
# It uses a configuration file that includes the details for the estimations,
# including the growth rates of each subclone.
# User input defines the range and increment of values to try.
# The result file is saved to the user defined folder with the file name
# 'm_{mmin}_{mmax}_{minc}_d_{dmin}_{dmax}_{dinc}_a_{amin}_{amax}_{ainc}_f_{fmin}_{fmax}_{finc}_l_{lmin}_{lmax}_{linc}.csv'
##########

def main(config, m_range, m_incr, d_range, d_incr, a_range, a_incr, f_range, f_incr, l_range, l_incr, save_path):
    # Create an estimator with the configuration file
    es = Estimator(config)

    # Print information
    print("m:", m_range, m_incr)
    print("d:", d_range, d_incr)
    print("a:", a_range, a_incr)
    print("f:", f_range, f_incr)
    print("l:", l_range, l_incr)
    print("Number of grid elements:", len(np.arange(m_range[0], m_range[1], m_incr))*len(np.arange(d_range[0], d_range[1], d_incr))*len(np.arange(a_range[0], a_range[1], a_incr))*len(np.arange(f_range[0], f_range[1], f_incr))*len(np.arange(l_range[0], l_range[1], l_incr)), flush=True)

    # Get the estimates for all parameters, given the growth rates and k and m relationship in the config file
    start_time = time.time()
    results = es.all_mice_b6(m_range, m_incr, d_range, d_incr, a_range, a_incr, f_range, f_incr, l_range, l_incr,
                             g1=es.params["g1"], g11=es.params["g11"], km_eq=[es.params["km_intercept"], es.params["km_slope"]])
    print("Time:", time.time() - start_time)

    # Save results
    results.to_csv("{}m_{}_{}_{}_d_{}_{}_{}_a_{}_{}_{}_f_{}_{}_{}_l_{}_{}_{}.csv".format(
        save_path, m_range[0], m_range[1], m_incr, d_range[0], d_range[1], d_incr, a_range[0], a_range[1], a_incr, 
        f_range[0], f_range[1], f_incr, l_range[0], l_range[1], l_incr), index=False)

    print(results)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", type=str, action="store",
                        help="The configuration file that includes the details for the estimation." \
                        "See the readme for details.")
    parser.add_argument("--mmin", type=float, action="store", default=-0.1,
                        help="The minimum value of m to evaluate.")
    parser.add_argument("--mmax", type=float, action="store", default=-0.1,
                        help="The maximum value of m to evaluate.")
    parser.add_argument("--minc", type=float, action="store", default=0,
                        help="The amount to increment between m values to evaluate.")
    parser.add_argument("--dmin", type=float, action="store", default=-1,
                        help="The minimum value of d to evaluate.")
    parser.add_argument("--dmax", type=float, action="store", default=-1,
                        help="The maximum value of d to evaluate.")
    parser.add_argument("--dinc", type=float, action="store", default=0,
                        help="The amount to increment between d values to evaluate.")
    parser.add_argument("--amin", type=float, action="store", default=1,
                        help="The minimum value of a to evaluate.")
    parser.add_argument("--amax", type=float, action="store", default=1,
                        help="The maximum value of a to evaluate.")
    parser.add_argument("--ainc", type=float, action="store", default=0,
                        help="The amount to increment between a values to evaluate.")
    parser.add_argument("--fmin", type=float, action="store", default=1,
                        help="The minimum value of f to evaluate.")
    parser.add_argument("--fmax", type=float, action="store", default=1,
                        help="The maximum value of k to evaluate.")
    parser.add_argument("--finc", type=float, action="store", default=0,
                        help="The amount to increment between f values to evaluate.")
    parser.add_argument("--lmin", type=float, action="store", default=0.001,
                        help="The minimum value of l to evaluate.")
    parser.add_argument("--lmax", type=float, action="store", default=0.001,
                        help="The maximum value of l to evaluate.")
    parser.add_argument("--linc", type=float, action="store", default=0,
                        help="The amount to increment between l values to evaluate.")
    parser.add_argument("-s", "--save_path", type=str, action="store",
                        help="The path in which to save the result file.")
    args = parser.parse_args()

    # If we only have one value we want to run, we need to make the bounds to be able to make an array with one value
    if args.mmin == args.mmax or args.minc == 0: 
        args.mmax = args.mmin
        args.minc = 0.1
    if args.dmin == args.dmax or args.dinc == 0: 
        args.dmax = args.dmin
        args.dinc = 0.1
    if args.amin == args.amax or args.ainc == 0: 
        args.amax = args.amin
        args.ainc = 0.1
    if args.fmin == args.fmax or args.finc == 0: 
        args.fmax = args.fmin
        args.finc = 0.1
    if args.lmin == args.lmax or args.linc == 0: 
        args.lmax = args.lmin
        args.linc = 0.1

    # Call main with the parameters
    # (Adding inc to the max so the range calls include the max value)
    main(args.config, 
         [args.mmin, args.mmax+args.minc], args.minc,
         [args.dmin, args.dmax+args.minc], args.dinc,
         [args.amin, args.amax+args.ainc], args.ainc, 
         [args.fmin, args.fmax+args.finc], args.finc,
         [args.lmin, args.lmax+args.linc], args.linc,
         args.save_path)