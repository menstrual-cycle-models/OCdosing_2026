MATLAB Code for "Modeling the effects of hormonal oral contraceptive dosing on contraceptive efficacy" 

Authors: Ruby Kim, Lihong Zhao, Lucy S. Oremland, Natasha Wozniak, Heather Z. Brooks, Elissa J. Schwartz, Erica J. Graham, and Lisette de Pillis

Getting Started
    Run startup.m to add necessary file paths.

Repository Structure
- main_driver.m
	- Script for running model with specified OC dosing conditions
- EE_NE_max_movstd_heatmaps.m
	- Script for visualizing long-term behavior (max & max of movstd) with specified treatment using heatmaps
- EE_NE_period_heatmaps.m
	- Script for visualizing long-term behavior (period) with specified treatment using heatmaps
- EE_NE_fits.m
	- Script for fitting pharmacokinetic parameters to EE and NE data
- helper_functions/
	- rhs.m
		- Function returning model for use with ODE solver
	- IC_rel_to_LH_peak.m
		- Function returning initial conditions on day relative to LH peak 
	- steadyTreatment.m
		- Function returning long-term behavior on specified treatment; For use with heatmaps or other analyses of long-term behavior
- data/ (Note: Data from simulations are saved to this folder.)
	- pars.xlsx
		- Contains parameters for model (order matters)
	- fda_ee_ne_data.mat
		- Contains ee_data (18x2), ee_error (18x1), pp_data (18x2), pp_error (18x1)
- figures/ (Note: Figures from simulations are saved to this folder. The folder will be created if it doesn't already exist.)
