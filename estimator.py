import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error
import json
from itertools import product
from src import pyxfunc

##########
# Class to generate and get errors of simulated curves.
# Takes as input a configuration file. Configuration files used are in config_files/.
##########

class Estimator():
    def __init__(self, config_file):
        with open(config_file, "r") as f:
            self.params = json.load(f)
        self.growth_df = pd.read_csv(self.params["growth_file"])
        if "t_init" in self.params.keys(): self.t_init = self.params["t_init"]
        else: self.t_init = 0
        self.dt = self.params["dt"]
        self.cutoff_0_raw =  self.params["cutoff_0_raw"]
        self.cutoff_0_percent = self.params["cutoff_0_percent"]
        if "max_loser_subline_percent" in self.params.keys(): self.max_loser_subline_percent = self.params["max_loser_subline_percent"]
        if "delay" in self.params.keys(): self.delay = self.params["delay"]
        exclude = [456, 458, 461, 464, 471, 474, 478, 482, 483, 426, 428, 429, 430, 434, 437, 438, 442, 443, 451, 642, 662]
        self.growth_df = self.growth_df[~(self.growth_df["id"].isin(exclude))]

    # Get the initial sizes of the sublines for the given mouse and the appropriate ratio for that group.
    # @param group The name of the group the mouse is in to get the initial subclone ratio
    # @param mid The id of the mouse to get the initial total tumor size
    # @return A three element list. The first element is the size of C1, the second element is the size of C11,
    #         the third element is the size of the T cell population.
    def get_init(self, group, mid):
        total_init = self.growth_df[(self.growth_df["id"]==mid) & (self.growth_df["day"]==7)]["size"].tolist()[0]
        if "A1" in group or "1C" in group: return [total_init, 0.0, self.t_init*total_init]
        if "A2" in group or "2C" in group: return [0.8*total_init, 0.2*total_init, self.t_init*total_init]
        if "A3" in group or "3C" in group: return [0.5*total_init, 0.5*total_init, self.t_init*total_init]
        if "A4" in group or "4C" in group: return [0.2*total_init, 0.8*total_init, self.t_init*total_init]
        if "A5" in group or "5C" in group: return [0.0, total_init, self.t_init*total_init]
        if "B1" in group: return [total_init, 0.0, 0.0]
        if "B2" in group: return [0.8*total_init, 0.2*total_init, 0.0]
        if "B3" in group: return [0.5*total_init, 0.5*total_init, 0.0]
        if "B4" in group: return [0.2*total_init, 0.8*total_init, 0.0]
        if "B5" in group: return [0.0, total_init, 0.0]

    # For all mice in either B1 or B5 and all growth rates passed, get the results of fitting 
    # the simulation to the the data.
    # @param group The name of the group to evaluate
    # @param grs The list of growth rates to evaluate
    # @return A pandas dataframe with the results of each growth rate evaluation for each mouse
    def nude_pure_group(self, group, grs):
        group_df = self.growth_df[self.growth_df["group"] == group]
        results = [[0]]*(len(group_df["id"].unique())*len(grs))
        idx = 0
        for g_idx in range(len(grs)):
            # eq_params = g1, g11, m, k, c, a, f, b, l -- all but one of the growth rates is 0 (l=0.01 to avoid divide by zero)
            if "B1" in group: # C1
                eq_params = [grs[g_idx], 0, 0, 0, 0, 0, 0, 0, 0.01]
                should_win = 0
            else: # C11
                eq_params = [0, grs[g_idx], 0, 0, 0, 0, 0, 0, 0.01]
                should_win = 1
            for mid in group_df["id"].unique():
                init = self.get_init(group, mid)
                end_time = max(group_df[group_df["id"]==mid]["day"])

                sol = pyxfunc.run(init[0], init[1], init[2], 7.0, end_time, np.inf, 
                                    self.cutoff_0_raw, self.cutoff_0_percent, init[2], *eq_params)

                day_idxs = [(d-7)*(1/self.dt) for d in self.growth_df[self.growth_df["id"] == mid]["day"].tolist()]
                sizes = self.growth_df[self.growth_df["id"] == mid]["size"].tolist()
                err_win = pyxfunc.get_error(day_idxs, sizes, sol, should_win, self.max_loser_subline_percent)
                results[idx] = [group, mid, eq_params[0], eq_params[1], err_win[0], err_win[1]]
                idx += 1
        results = pd.DataFrame(results, columns=["group", "id", "g1", "g11", "error", "winner"])
        return results

    # Evaluate the passed growth rates on all mice in groups B1 and B5.
    # @param grs The growth rates to evaluate
    # @return A pandas dataframe with the results of evaluating the growth rates.
    def nude_pure_all(self, grs):
        results = []
        for group in ["Grp. B1 nude (100% C1)", "Grp. B5 nude (100% C11)"]:
            results += [self.nude_pure_group(group, grs)]
        results = pd.concat(results)
        return results

    # For all mice in groups B2, B3, or B4 and all growth rates, k values, and m values, 
    # get the results of fitting the simulation to the the data.
    # @param group The name of the group to evaluate
    # @param g1_list The list of C1 growth rates to evaluate
    # @param g11_list The list of C11 growth rates to evaluate
    # @param m_list The list of m values to evaluate
    # @param k_list The list of k values to evaluate
    # @return A pandas dataframe with the results of each parameter evaluation for each mouse
    def nude_admix_group(self, group, g1_list, g11_list, m_list, k_list):
        group_df = self.growth_df[self.growth_df["group"] == group]
        end_time = max(group_df["day"])
        should_win = 0 # C1 should always win nude admixes
        results = [[0]]*len(group_df["id"].unique())*len(g1_list)*len(g11_list)*len(m_list)*len(k_list)
        idx = 0
        
        for params in product(g1_list, g11_list, m_list, k_list):
            eq_params = [*params, 0, 0, 0, 0, 0.01]
            for mid in group_df["id"].unique():
                init = self.get_init(group, mid)
                sol = pyxfunc.run(init[0], init[1], init[2], 7.0, end_time, np.inf, 
                                    self.cutoff_0_raw, self.cutoff_0_percent, init[2], *eq_params)
                idxs = [(d-7)*(1/self.dt) for d in self.growth_df[self.growth_df["id"] == mid]["day"].tolist()]
                sizes = self.growth_df[self.growth_df["id"] == mid]["size"].tolist()
                err_win = pyxfunc.get_error(idxs, sizes, sol, should_win, self.max_loser_subline_percent)
                results[idx] = [group, mid, *params, err_win[0], err_win[1]]
                idx += 1
        results = pd.DataFrame(results, columns=["group", "id", "g1", "g11", "m", "k", "error", "winner"])

        return results


    # Evaluate the passed growth rates on all mice in groups B2, B3, and B4.
    # @param g1_list The list of C1 growth rates to evaluate
    # @param g11_list The list of C11 growth rates to evaluate
    # @param m_list The list of m values to evaluate
    # @param k_list The list of k values to evaluate
    # @return A pandas dataframe with the results of evaluating the parameter combinations.
    def nude_admix_all(self, g1_list, g11_list, m_list, k_list):
        results = []
        for group in ["Grp. B2 nude (80% C1; 20% C11)", "Grp. B3 nude (50% C1; 50% C11)", "Grp. B4 nude (20% C1; 80% C11)"]:
            results += [self.nude_admix_group(group, g1_list, g11_list, m_list, k_list)]
        results = pd.concat(results)
        return results

    # For any group in B6 mice, and parameters passed, 
    # get the results of fitting the simulation to the the data.
    # @param group The name of the group to evaluate
    # @param gkm_pairs The list of C1 growth rate, C11 growth rate, m, and k tuples to evaluate
    # @param c_list The list of c values to evaluate
    # @param a_list The list of a values to evaluate
    # @param f_list The list of f values to evaluate
    # @param b_list The list of b values to evaluate
    # @param l_list The list of l values to evaluate
    # @return A pandas dataframe with the results of each parameter evaluation for each mouse
    def b6_group(self, group, gkm_pairs, c_list, a_list, f_list, b_list, l_list):
        group_df = self.growth_df[self.growth_df["group"] == group]
        results = [[0]*10]*(len(gkm_pairs)*len(c_list)*len(a_list)*len(f_list)*len(b_list)*len(l_list)*len(group_df["id"].unique()))
        idx=0
        if ("A1" in group or "A2" in group or "1C" in group or "2C" in group):
            should_win = 0 # C1
        else: 
            should_win = 1 # C11
        for params in product(gkm_pairs, c_list, a_list, f_list, b_list, l_list):
            eq_params = [params[0][0], params[0][1], params[0][2], params[0][3], *params[1:]]            
            for mid in group_df["id"].unique():
                init = self.get_init(group, mid)
                end_time = max(group_df[group_df["id"]==mid]["day"])
                sol = pyxfunc.run(init[0], init[1], init[2], 7.0, end_time, self.delay, 
                                    self.cutoff_0_raw, self.cutoff_0_percent, init[2], *eq_params)
                idxs = [(d-7)*(1/self.dt) for d in self.growth_df[self.growth_df["id"] == mid]["day"].tolist()]
                sizes = self.growth_df[self.growth_df["id"] == mid]["size"].tolist()
                err_win = pyxfunc.get_error(idxs, sizes, sol, should_win, self.max_loser_subline_percent)
                results[idx] = [group, mid, params[0][0], params[0][1], params[0][2], params[0][3], *params[1:], err_win[0], err_win[1]]
                idx += 1
        results = pd.DataFrame(results, columns=["group", "id", "g1", "g11", "m", "k", "c", "a", "f", "b", "l", "error", "winner"])
        return results  


    # Evaluate the passed growth rates on B6 mice.
    # @param gkm_pairs The list of C1 growth rate, C11 growth rate, m, and k tuples to evaluate
    # @param c_list The list of c values to evaluate
    # @param a_list The list of a values to evaluate
    # @param f_list The list of f values to evaluate
    # @param b_list The list of b values to evaluate
    # @param l_list The list of l values to evaluate
    # @return A pandas dataframe with the results of evaluating the parameter combinations.
    def b6_all(self, gkm_pairs, c_list, a_list, f_list, b_list, l_list):
        results = []
        for group in ["Grp. A1 B6 (100% C1)", "Grp. A2 B6 (80% C1; 20% C11)", "Grp. A3 B6 (50% C1; 50% C11)", "Grp. A4 B6 (20% C1; 80% C11)", "Grp. A5 B6 (100% C11)", 
                      "1C", "2C", "3C", "4C", "5C"]:
            results += [self.b6_group(group, gkm_pairs, c_list, a_list, f_list, b_list, l_list)]
        results = pd.concat(results)
        return results