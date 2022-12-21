Given the output from Galform runs, here I have gathered the codes to reproduce the observable functions used in the Elliott et al. (2021) paper. Apart from the observable functions, their specific bins are also taken from this paper.

- calibration.py: given the Galform output of a single model, it generates the observable functions. It also save the bins and the values of the observable functions (data to be used by the emulator) we consider for calibration.
- calibration_plots.py: make the calibration plots with the data obtained from calibration.py

- calibration_models.py: compute and graphic the calibration plots in the same way as calibration.py and calibration_plots.py, but for 1 or more than 1 models to compare different runs.

- calibration_dispersion.py: compute and graphic the calibration plots for 250Mpc/h subboxes of the UNIT simulation to analyse the dispersion of these plots.
- calibration_dispersion_ratio.py: given the output from calibration.dispersion.py it computes the ratio between the subboxes and the "real value" of the whole 1000Mpc/h volume of the UNIT simulation. This is a visible way to quantify the dispersion, which could be useful to define the range of these plots where we are going to train the emulator (where the dispersion is lower than 10% f.e.).