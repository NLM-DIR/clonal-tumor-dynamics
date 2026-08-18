# Mathematical and experimental models of interactions among subclones and immunity in melanoma

Code for the paper Mathematical and experimental models of interactions among subclones and immunity in melanoma, Hirsch et al.

## Reproduce results

### Create the Conda environment

In order to ensure all installations are in place, we have included a Conda environment. In order to create it, run

``conda create --file environment.yml``

Activate the environment using

``conda activate tumor-dynamics``

### Build C++ source:

First, you will need to build the C++ source code.

In src/ folder, run

``python setup.py build_ext --inplace``

### Generate results using grid search

Running the grid search requires a configuration file. Example files are included in ``config_files/`` in this repository.

In order to run a single set of parameters, call the ``runner`` files:

- ``runner_g.py`` will run a grid search for the growth rates of the nude mice with a single subclone.
- ``runner_km.py`` will run a grid search for $k$ and $m$ for the nude mice with an admixture of the subclones. This uses the optimal growth rates estimated from ``runner_g.py``, which must be run first.
- ``runner_t.py`` will run a grid search for $c$, $a$, $f$, and $b$ for the B6 mice. This uses the optimal growth rate and $k$ and $m$ values estimated from ``runner_km.py``, which must be run first.
- ``runner_t_3d.py`` will run a grid search for $c$, $a$, and $f$ for the B6 mice. This is used to generate the 3D figures with the outcome of the simulation when varying $c$, $a$, and $f$ while holding the other parameters constant at the optimal values for each mouse. This uses the optimal parameter estimations estimated from ``runner_t.py``, which must be run first.

The ``runner`` files use ``estimator.py``, which calls the C++ code in the ``src`` folder. See how to compile the C++ code above.

This code was run on the NIH Biowulf cluster, which uses an in house SLURM command SWARM for large batches of jobs. The ``make_km_swarm.sh``, ``make_t_swarm.sh``, and ``make_t_3d_swarm.sh`` files generate the SWARM files (``swarm_km.swarm``, ``swarm_t.swarm``, and ``swarm_t_3d.swarm``) for a range of parameters. These SWARM files contain the python commands used to run the grid search. These commands can then be copied into a standard SBATCH script or run locally. The scripts call the ``runner_km.py``, ``runner_t.py``, and ``runner_t_3d.py`` files, respectively.

### Analyze grid search results

In order to analyze the results of the grid search and generate their respective figures, use the ``analyze`` Jupyter notebooks:

- ``analyze_growth_rates.ipynb`` loads the results from ``runner_g.py``. It finds the optimal growth rate values for each nude mouse with a single subclone and plots the resulting curves (100% proliferative and 100% invasive in Figure 3d). It saves the optimal growth rates for use in ``runner_km.py``.
- ``analyze_km.ipynb`` loads the results from ``runner_km.py``. It finds the optimal growth rate, $k$, and $m$ tuples for the nude mice with an admixture of subclones. It saves the optimal growth rate, $k$, and $m$ tuples for use in ``runner_t.py``.
- ``analyze_t.ipynb`` loads the results from ``runner_t.py``. It finds the optimal growth rate, $k$, $m$, $c$, $a$, $f$, $b$ tuples for each B6 immunocompetent mouse. It saves the optimal tuples for use in ``runner_t_3d.py``.

### Generate figures

- ``Figure_2a_4a_plot_data.ipynb`` generates Figures 2a and 4a
- ``Figure_3a_Table_S1.ipynb`` generates Figure 3a
- ``Figure_3b.ipynb`` generates Figure 3b
- ``Figure_3c_3d.ipynb`` generates Figures 3c and 3d
- ``Figure_5a_5b.ipynb`` generates Figures 5a, 5b, and Supplementary Figure 8
- ``Figure_6_S3_S4_S5_S6.ipynb`` generates Figure 6, Supplementary Figures 3, 4, 5, and 6.
- ``Figure_S9.ipynb`` generates Supplementary Figure 9

(Figures 2b and 4b were generated using Excel Supplementary Table 3).

### Generate supplementary tables

- ``Figure_3a_Table_S1.ipynb`` generates Supplementary Table 1
- ``Table_S2.ipynb`` generates Supplementary Table 2

### Get other results

- ``get_percent_decrease.ipynb`` calculates the percentage decrease between the growth rates of the nude immunocompromised and B6 immunocompetent mice.
