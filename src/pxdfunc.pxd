from libcpp cimport bool
from libcpp.vector cimport vector


cdef extern from "cfunc.cpp":
    cdef vector[double] _one_step(double c1, double c11, double tcell, bool delay_flag, double end_time, double init_tcell,
                    double g1, double g11, double m, double k, double d, double a, double f, double b, double l)
    cdef vector[double] _one_step_it(double c1, double c11, double tcell, bool delay_flag, double end_time, double init_tcell,
                    double g1, double g11, double m, double k, double d, double a, double f, double b, double l, double h)
    cdef vector[vector[double]] _run(double init_c1, double init_c11, double init_tcell, double start_time, double end_time, 
                            double delay, double cutoff_0_raw, double cutoff_0_percent, double base_tcell,
                            double g1, double g11, double m, double k, double d, double a, double f, double b, double l)
    cdef vector[vector[double]] _run_it(double init_c1, double init_c11, double init_tcell, double start_time, double end_time, 
                            double delay, double cutoff_0_raw, double cutoff_0_percent, double base_tcell,
                            double g1, double g11, double m, double k, double d, double a, double f, double b, double l, double h)
    cdef vector[double] _get_error(vector[int] idxs, vector[double] true_growth, vector[vector[double]] pred_growth, 
                            int should_win, double max_loser_perc)
