import pandas as pd
import numpy as np
from estimator import Estimator
import argparse
import time

##########
# This script estimates the interaction terms k and m for the experiments with
# mixtures of the proliferative and invasive subclones in nude mice.
# It uses a configuration file that includes the details for the estimations,
# including the growth rates of each subclone.
# User input defines the range and increment of values to try.
# The result file is saved to the user defined folder with the file name
# 'k_{kmin}_{kmax}_{kinc}_m_{mmin}_{mmax}_{minc}.csv'
##########


def main(config, k_range, k_incr, m_range, m_incr, w, save_path):
    # Create an estimator with the configuration file
    es = Estimator(config)

    # Print information
    print("k:", k_range, k_incr)
    print("m:", m_range, m_incr)
    print("Determine subclone winner only?",  w)
    print("Number of grid elements:", len(np.arange(k_range[0], k_range[1], k_incr))*len(np.arange(m_range[0], m_range[1], m_incr)), flush=True)

    # Get the k and m estimates for all admixtures
    start_time = time.time()
    results = es.all_mice_nude_admix(k_range, k_incr, m_range, m_incr, g1=0.135, g11=0.107, win_only=w)
    print("Time:", time.time() - start_time)

    # Save results
    results.to_csv("{}k_{}_{}_{}_m_{}_{}_{}.csv".format(save_path, k_range[0], k_range[1], k_incr, m_range[0], m_range[1], m_incr), index=False)
    print(results)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", type=str, action="store",
                        help="The configuration file that includes the details for the estimation." \
                        "See the readme for details.")
    parser.add_argument("--kmin", type=float, action="store",
                        help="The minimum value of k to evaluate.")
    parser.add_argument("--kmax", type=float, action="store",
                        help="The maximum value of k to evaluate.")
    parser.add_argument("--kinc", type=float, action="store",
                        help="The amount to increment between k values to evaluate.")
    parser.add_argument("--mmin", type=float, action="store",
                        help="The minimum value of m to evaluate.")
    parser.add_argument("--mmax", type=float, action="store",
                        help="The maximum value of k to evaluate.")
    parser.add_argument("--minc", type=float, action="store",
                        help="The amount to increment between m values to evaluate.")
    parser.add_argument("-w", action="store_true",
                        help="Use this flag to only get results for which subline wins (error values reported will be incorrect)")
    parser.add_argument("-s", "--save_path", type=str, action="store",
                        help="The path in which to save the result file.")
    args = parser.parse_args()
    main(args.config, [args.kmin, args.kmax], args.kinc, [args.mmin, args.mmax], args.minc, args.w, args.save_path)