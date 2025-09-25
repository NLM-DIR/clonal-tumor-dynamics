# Mathematical and experimental models of interactions among subclones and immunity in melanoma

This repository contains the code to reproduce the results and figures of Mathematical and experimental models of interactions among subclones and immunity in melanoma.

<b>Study abstract:</b>

Melanoma tumors often contain heterogeneous subclones, typically with proliferative or invasive phenotypes. These subclones interact with each other and the immune system to shape tumor progression. Understanding these interactions between subclones and the tumor microenvironment is a fundamental challenge for melanoma treatment. To address this challenge, we developed dedicated *in vivo* mouse models and use experimental measurements afforded by them to inform mathematical models of these interactions. For the experimental models, admixtures of proliferative and invasive subclones of varying proportions were implanted into immunocompromised and immunocompetent mice, and tumor growth and subclone composition were tracked over time. Based on the results of the experimental study, we built a Lotka–Volterra-based model interaction between clones and extend it with a piecewise equation for T cell dynamics to capture immune–tumor interactions. By identifying key subclone–immune interactions, this study provides insights into how intratumor heterogeneity and immune pressure influence melanoma progression and suggests strategies for optimizing therapeutic timing and adaptive treatment.

This repository contains code to perform exhaustive searches to find $g_P$, $g_I$, $k$, $m$, $d$, $a$, $f$, and $l$, to perform analysis on those results, and plot the results as shown in the manuscript.

The model generates curves based on the following set of differential equations:

$$
\frac{dP}{dt} = g_P P(t) + kP(t)I(t) - dP(t)T(t)
\\
\frac{dI}{dt} = g_I I(t) + mI(t)P(t)
$$
$$
t < t_1: \frac{dT}{dt} = 0 \\
t_1 < t < t_2: \frac{dT}{dt} = \frac{a}{l+P(t)}P(t)T(t) - fT(t) \\
t_2 < t < t_3: \frac{dT}{dt} = -fT(t) \\
t_3 < t: \frac{dT}{dt} = 0 \\
$$

where $t_1$ is the time at which the tumor reaches size 2.6, $t_2$ is the time at which the tumor either reaches 7.25 or 0.01, and $t_3$ is the time at which the T cell population reaches its initial size.

Details of the model can be found in the publication.

## Steps to reproduce the results

Code to reproduce the figures and get values and statistics reported in the paper are in the Python scripts and Jupyter Notebooks in this repository. The notebooks are titled based on the figure(s) they produce.

The user interface to run the exhaustive search (needed to reproduce the results figures), is included in:
- `exhaustive_growth.py` runs the exhaustive search for the subclonal growth rates
- `exhaustive_subclone.py` runs the exhaustive search for the subclone interaction terms: $k$ and $m$
- `exhaustive_tcell.py` runs the exhaustive search for the parameters of the T cell model: $m$, $d$, $a$, $f$, and $l$

The `estimator.py` file includes the `Estimator` class that contains the code that actually performs the exhaustive searches. It also reads and stores the values from the configuration files for ease of access.

The Jupyter notebooks to evaluate and plot the results take a folder for the results and will load all files in that folder. This is to accommadate multiple runs of different parameter ranges which make it easy to use break up the parameter ranges to run in parallel.

### Examples to generate results

Run exhaustive search for the growth rates of the nude mice with growth rates in the range of 0.05 to 0.15 with increments of 0.01, saving to `results/growth/`:

```
python exhaustive_growth.py -c config_subclone.json --gmin 0.05 --gmax 0.15 --ginc 0.01 -s results/growth/
```

Run exhaustive search for $k$ and $m$ evaluating *which subline wins only* (in order to reproduce Figure 3c; note that running with the `-w` flag in order to get results for the winning subline only will result in outputting incorrect error values), with values in the range of -0.2 to 0.2 with increments of 0.05, saving to `results/subclone_winners/`:

```
python exhaustive_subclone.py -c config_subclone.json --kmin -0.2 --kmax 0.2 --kinc 0.05 --mmin -0.2 --mmax 0.2 --minc 0.05 -w -s results/subclone_winners/
```

Run exhaustive search for $k\in[0, 0.2]$ and $m\in[-0.2,0)$ with increments of 0.01, saving to `results/subclone/`:

```
python exhaustive_subclone.py -c config_subclone.json --kmin 0 --kmax 0.2 --kinc 0.01 --mmin -0.2 --mmax -0.01 --minc 0.01 -s results/subclone/
```

Run exhaustive search for $m\in[-0.2, 0)$ with increments of 0.01, $d\in(0, 10]$ with increments of 0.5, $a\in[0.1, 2.5]$ with increments of 0.1, $f\in[0.1, 2.5]$ with increments of 0.1, and $l\in[0.0001, 0.001]$ with increments of 0.0001, with $g_1$, $g_11$, and $k\sim m$ defined in config_tcell.json, saving to `results/tcell/`.

```
python exhaustive_tcell.py -c config_tcell.json --mmin -0.2 --mmax -0.01 --minc 0.01 --dmin 0.5 --dmax 10 --dinc 0.5 --amin 0.1 --amax 2.5 --ainc 0.1 --fmin 0.1 --fmax 2.5 --finc 0.1 --lmin 0.0001 --lmax 0.001 --linc 0.0001 -s results/tcell/
```

## Configuration files

Both the code to run the exhaustive search and the code to reproduce the figures utilizes the Estimator class. This class takes a configuration file for user input parameter values. The configuration file should be in JSON format and must include the the following fields for any instantiation of an Estimator object:

- `growth_file`: (string) The path to the growth data. 
- `t_init`: (float) The initial size of the T cell population. 
- `t_incr_size`: (float) The normalized size of the tumor that triggers the change in the T cell dynamics from $dT/dt=0$ to $dT/dt=\frac{a}{l+P(t)}P(t)T(t)-f(t)$. 
- `dt`: (float) The increment size for generating the curve defined by the ODE
- `cutoff_0`: (float) The size at which a subclone is considered as having reached a population of 0
- `error_type`: (string) What error calculation to use when performing the exhaustive search: either "MSE" or "R2"

Additional fields can be included:
- `g1`: (float) The growth rate of subclone C1 to use when calling one_mouse_nude_admix, all_mice_nude_admix, one_mouse_b6, or all_mice_b6 if g1 is not provided as a parameter to the method.
- `g11`: (float) The growth rate of subclone C11 to use when calling one_mouse_nude_admix, all_mice_nude_admix, one_mouse_b6, or all_mice_b6 if g11 is not provided as a parameter to the method.
- `km_intercept`: (float) The intercept of the linear relationship between $k$ and $m$ ($y$ in $k=xm+y) to use when calling one_mouse_b6 or all_mice_b6 if km_eq is not provided as a parameter to the method.
- `km_slope`: (float) The slope of the linear relationship between $k$ and $m$ ($x$ in $k=xm+y) to use when calling one_mouse_b6 or all_mice_b6 if km_eq is not provided as a parameter to the method.

Example configuration files are provided: `config_subclone.json` includes only the required fields; `config_tcell.json` includes optional fields that define the parameter space used in our experiments to estimate the T cell parameters.

Note that the JSON format does not have an option for comments.

## Reproducing the Figures

- Figure 1: Created in Illustrator
- Figure 2:
    - 2a: `Figure_2a_4a.ipynb`
    - 2b: `data/VAF_medians.xlsx`
- Figure 3:
    - 3a: `Figure_3a.ipynb`
    - 3b: `Figure_3b.ipynb`
    - 3c: `Figure_3c.ipynb`
    - 3d: `Figure_3d.ipynb`
    - 3e: `Figure_3e.ipynb`
- Figure 4:
    - 4a: `Figure_2a_4a.ipynb`
    - 4b: `data/VAF_medians.xlsx`
- Figure 5:
    - 5a: Created in Illustrator
    - 5b: `Figure_5b.ipynb`
- Figure 6: `Figure_3e_6.ipynb`
- Figure 7: `Figure_7_S3.ipynb`
- Figure S1: `Figure_S1_S2.py`
- Figure S2: `Figure_S1_S2.py`
- Figure S3: `Figure_7_S3.ipynb`

## Data

Data files needed to reproduce our results are included in the data/ folder. This includes:
- `growth_data_df.csv`: The raw growth measurements of each tumor.
- `growth_data_scaled.csv`: The normalized growth measurements of each tumor. Each measurement for each tumor is divided by the first measurement of that tumor so that every tumor has a first measurement equal to 1.
- `VAF_medians`: The medians of the VAF measurements (used for Figure 2b and 4b)
- `VAF/{Group ID}_VAF.csv`: The VAF measurements for that Group ID.

This data can also be found in the Supplementary Material of the publication.

## Citation
