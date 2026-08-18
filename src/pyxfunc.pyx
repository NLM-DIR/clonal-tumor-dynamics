# distutils: language=c++

from libcpp cimport bool
from libcpp.vector cimport vector
from pxdfunc cimport _one_step, _run, _get_error, _one_step_it, _run_it

def one_step(double c1, double c11, double tcell, bool delay_flag, double end_time, double init_tcell,
                    double g1, double g11, double m, double k, double d, double a, double f, double b, double l):
    cdef vector[double] result = _one_step(c1, c11, tcell, delay_flag, end_time, init_tcell,
                    g1, g11, m, k, d, a, f, b, l)
    return result


def one_step_it(double c1, double c11, double tcell, bool delay_flag, double end_time, double init_tcell,
                    double g1, double g11, double m, double k, double d, double a, double f, double b, double l, double h):
    cdef vector[double] result = _one_step_it(c1, c11, tcell, delay_flag, end_time, init_tcell,
                    g1, g11, m, k, d, a, f, b, l, h)
    return result

def run(double init_c1, double init_c11, double init_tcell, double start_time, double end_time, 
          double delay, double cutoff_0_raw, double cutoff_0_percent, double base_tcell,
          double g1, double g11, double m, double k, double d, double a, double f, double b, double l):
    cdef vector[vector[double]] result = _run(init_c1, init_c11, init_tcell, 
                                             start_time, end_time, delay, cutoff_0_raw, cutoff_0_percent, base_tcell,
                                             g1, g11, m, k, d, a, f, b, l)
    return result

def run_it(double init_c1, double init_c11, double init_tcell, double start_time, double end_time, 
          double delay, double cutoff_0_raw, double cutoff_0_percent, double base_tcell,
          double g1, double g11, double m, double k, double d, double a, double f, double b, double l, double h):
    cdef vector[vector[double]] result = _run_it(init_c1, init_c11, init_tcell, 
                                             start_time, end_time, delay, cutoff_0_raw, cutoff_0_percent, base_tcell,
                                             g1, g11, m, k, d, a, f, b, l, h)
    return result

def get_error(vector[int] idxs, vector[double] true_growth, vector[vector[double]] pred_growth, int should_win, double max_loser_perc):
    cdef vector[double] result = _get_error(idxs, true_growth, pred_growth, should_win, max_loser_perc)
    return result