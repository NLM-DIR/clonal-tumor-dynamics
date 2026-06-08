#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <cmath>
using namespace std;

// Compile with python setup.py build_ext --inplace
// Any changes to function definitions here need to be reflected in pxdfunc and pyxfunc

/*
    Take one step in the function
    @param double c1 The size of subclone C1
    @param double c11 The size of subclone C11
    @param double tcell The size of the T cell population
    @param bool delay_flag Whether or not the T cell population is changing. 
                           True means the T cell population will not change. 
                           False means the T cell population will change following the dynamics.
    @param double end_time The end time of the simulation
    @param init_tcell The initial T cell population size
    @param double g1 The growth rate of C1
    @param double g11 The growth rate of C11
    @param double m Effect of C11 on C1
    @param double k Effect of C1 on C11
    @param double c Effect of T cells on C1
    @param double a Recruitment rate of T cells by C1
    @param double f Turnover rate of T cells
    @param double b The suppression of T cells by C1
    @param double l Michelis-Menten parameter
*/ 
vector<double> _one_step(double c1, double c11, double tcell, bool delay_flag, double end_time, double init_tcell,
                    double g1, double g11, double m, double k, double c, double a, double f, double b, double l) {
    double s = init_tcell*f;
    double tcell_next;
    if (delay_flag) {
        tcell_next = tcell;
    }
    else {
        tcell_next = tcell + (s+ (a/(l+sqrt(c1)))*sqrt(c1)*tcell - b*sqrt(c1)*tcell - f*tcell)*0.01;
    }
    double c1_next = c1 + (g1*c1 + k*sqrt(c1)*sqrt(c11) - c*c1*tcell)*0.01;
    double c11_next = c11 + (g11*c11 + m*sqrt(c11)*sqrt(c1))*0.01;
    return {c1_next, c11_next, tcell_next};
}


/*
    Run the simulation
    @param double init_c1 The initial size of subclone C1
    @param double init_c11 The initial size of subclone C11
    @param double init_tcell The inital size of the T cell population
    @param double start_time The start day of the simulation
    @param double end_time The end day of the simulation
    @param double delay The day that the T cell population begins to change
    @param double cutoff_0_raw The raw size of the subclone population at which to set the population size to 0
    @param double cutoff_0_percentage The percentage of the subclone in the tumor population at which to set the population size to 0
    @param double base_tcell The minimum size of the T cell population
    @param double g1 The growth rate of C1
    @param double g11 The growth rate of C11
    @param double m Effect of C11 on C1
    @param double k Effect of C1 on C11
    @param double c Effect of T cells on C1
    @param double a Recruitment rate of T cells by C1
    @param double f Turnover rate of T cells
    @param double b The suppression of T cells by C1
    @param double l Michelis-Menten parameter
    @return vector<vector<double>> results Vectors of the sizes of the populations
                                           results[0] is the size of C1
                                           results[1] is the size of C11
                                           results[2] is the size of the T cell population
                                           results[3] is the time values
*/ 
vector<vector<double>> _run(double init_c1, double init_c11, double init_tcell, double start_time, double end_time, 
                           double delay, double cutoff_0_raw, double cutoff_0_percent, double base_tcell,
                           double g1, double g11, double m, double k, double c, double a, double f, double b, double l) {
    double time = start_time;
    bool delay_flag = true;

    vector<vector<double>> results = {vector<double>{init_c1}, vector<double>{init_c11}, vector<double>{init_tcell}, vector<double>{time}}; 
    int idx = 0;

    while (time < end_time) {
        vector<double> curr_results;

        if (time >= delay) {
            delay_flag = false;
        }
        curr_results = _one_step(results.at(0).at(idx), results.at(1).at(idx), results.at(2).at(idx), delay_flag, end_time, base_tcell, g1, g11, m, k, c, a, f, b, l);

        if (curr_results[0] <= cutoff_0_raw || curr_results[0]/(curr_results[0]+curr_results[1]) <= cutoff_0_percent) {
            curr_results[0] = 0;
        }
        if (curr_results[1] <= cutoff_0_raw || curr_results[1]/(curr_results[0]+curr_results[1]) <= cutoff_0_percent) {
            curr_results[1] = 0;
        }
        if (curr_results[2] <= base_tcell) {
            curr_results[2] = base_tcell;
        }

        time = time + 0.01;
        idx = idx + 1;

        results.at(0).push_back(curr_results.at(0));
        results.at(1).push_back(curr_results.at(1));
        results.at(2).push_back(curr_results.at(2));
        results.at(3).push_back(time);
    }
    return results;
}

/*
    Get the error between the simulated growth and the true data
    @param vector<int> idxs The indices of the predicted growth vectors to use to calculate the error
    @param vector<double> true_growth A vector of the sizes of the true tumor growth
    @param vector<vector<double>> pred_growth A vector of the predicted growth of:
                                                pred_growth[0] C1
                                                pred_growth[1] C11
                                                pred_growth[2] T cells
                                                pred_growth[3] time
    @param int should_win Which subclone should dominate the tumor at the end point; 0 = C1, 1 = C11
    @param double max_loser_percent The maximum percent of the tumor that the non-dominating tumor can be at the end point
    @return vector<double> results
                                results[0] double mean squared error
                                results[1] double which subclone dominates; 0 = C1, 1 = C11, 2 = neither
*/
vector<double> _get_error(vector<int> idxs, vector<double> true_growth, vector<vector<double>> pred_growth, int should_win, double max_loser_perc) {

    vector<double> total_pred (pred_growth.at(0).size());
    std::transform(pred_growth.at(0).cbegin(), pred_growth.at(0).cend(), pred_growth.at(1).cbegin(), total_pred.begin(), std::plus<>{});

    vector<double> pred;
    for (unsigned long int i = 0; i < idxs.size(); i++) {
            pred.push_back(total_pred.at(idxs.at(i)));
    }

    double mse = 0;
    for (unsigned long int i=0; i<true_growth.size(); i++) {
        mse += (true_growth.at(i) - pred.at(i)) * (true_growth.at(i) - pred.at(i)) ;
    }
    mse = mse / true_growth.size();

    double end_c1 = pred_growth.at(0).at(pred_growth.at(0).size()-1);
    double end_c11 = pred_growth.at(1).at(pred_growth.at(1).size()-1);
    double winner = 2;
    if ((end_c1 > end_c11) & (end_c11 / (end_c1+end_c11) < max_loser_perc)) { 
        winner = 0;
    }
    else if ((end_c11 > end_c1) & (end_c1 / (end_c1+end_c11) < max_loser_perc)) {
        winner = 1;
    }
    if ((int)winner != should_win) {
        mse = INFINITY;
    }

    return {mse, winner};
}
