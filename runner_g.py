from estimator import Estimator
import argparse
import time

##########
# This script estimates the growth rates gP and gI for the experiments with
# only the proliferative or invasive subclones in nude mice.
# It uses a configuration file that includes the details for the estimations
# # User input defines the range and increment of values to try.
# The result file is saved to the user defined folder with the file name
# 'g_{gmin}_{gmax}_{ginc}.csv'
##########


def main(config, g_min, g_max, g_inc, save_path):
    # Create an estimator with the configuration file
    es = Estimator(config)

    g_precision = max(str(g_min)[::-1].find("."), str(g_max)[::-1].find("."), str(g_inc)[::-1].find("."))
    print("Precision:", g_precision)

    g_list = [round(g_min+g_inc*i, g_precision) for i in range(round((g_max-g_min)/g_inc+1))]

    # Print information
    print("g list:", g_list)
    print("Number of grid elements:", len(g_list), flush=True)

    # Get the g estimates for C1 and C11
    start_time = time.time()
    results = es.nude_pure_all(g_list)
    print("Time:", time.time() - start_time)

    # Save results
    results.to_csv("{}_g_{}_{}_{}.csv".format(save_path, g_min, g_max, g_inc), index=False)
    print(results)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", type=str, action="store",
                        help="The configuration file that includes the details for the estimation." \
                        "See the readme for details.")
    parser.add_argument("--gmin", type=float, action="store",
                        help="The minimum value of g to evaluate.")
    parser.add_argument("--gmax", type=float, action="store",
                        help="The maximum value of g to evaluate.")
    parser.add_argument("--ginc", type=float, action="store",
                        help="The amount to increment between g values to evaluate.")
    parser.add_argument("-s", "--save_path", type=str, action="store",
                        help="The path in which to save the result file.")
    args = parser.parse_args()

    # If we only have one value we want to run, we need to make the bounds to be able to make an array with one value
    if args.gmin == args.gmax or args.ginc == 0: 
        args.gmax = args.gmin
        args.ginc = 0.1

    main(args.config, args.gmin, args.gmax, args.ginc, args.save_path)