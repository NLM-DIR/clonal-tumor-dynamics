import numpy as np
from estimator import Estimator
import argparse
import time

##########
# This script estimates the growth rates of the experiments with 100% proliferative or 100% invasive subclones.
# It uses a configuration file that includes the details for the estimations.
# User input defines the range and increment of values to try.
##########


def main(config, g_range, g_incr, save_path):
    # Create an estimator with the configuration file
    es = Estimator(config)

    # Print information
    print("g:", g_range, g_incr)
    print("Number of grid elements:", len(np.arange(g_range[0], g_range[1], g_incr)), flush=True)

    # Get C1 (proliferative) growth rate
    print("C1 (proliferative)")
    start_time = time.time()
    results = es.all_mice_nude_pure("Grp. B1 nude (100% C1)", g_range, g_incr)
    print("Time:", time.time() - start_time)
    
    # Save results
    results.to_csv("{}C1_g_{}_{}_{}.csv".format(save_path, g_range[0], g_range[1], g_incr), index=False)
    print(results)

    # Get C11 (invasive) growth rate
    print("C11 (invasive)")
    start_time = time.time()
    results = es.all_mice_nude_pure("Grp. B5 nude (100% C11)", g_range, g_incr)
    print("Time:", time.time() - start_time)

    # Save results
    results.to_csv("{}C11_g_{}_{}_{}.csv".format(save_path, g_range[0], g_range[1], g_incr), index=False)
    print(results)

    # Get the admixture total tumor growth rates
    print("Admixtures")
    start_time = time.time()
    results = es.all_admix_nude_growth(g_range, g_incr)
    print("Time:", time.time() - start_time)

    # Save results
    results.to_csv("{}admix_g_{}_{}_{}.csv".format(save_path, g_range[0], g_range[1], g_incr), index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", type=str, action="store", 
                        help="The configuration file that includes the details for the estimation." \
                        "See the readme for details.")
    parser.add_argument("--gmin", type=float, action="store",
                        help="The minimum growth rate to evaluate.")
    parser.add_argument("--gmax", type=float, action="store",
                        help="The maximum growth rate to evaluate.")
    parser.add_argument("--ginc", type=float, action="store",
                        help="The amount to increment between growth values to evaluate.")
    parser.add_argument("-s", "--save_path", type=str, action="store",
                        help="The file in which to save the results.")
    args = parser.parse_args()
    main(args.config, [args.gmin, args.gmax], args.ginc, args.save_path)