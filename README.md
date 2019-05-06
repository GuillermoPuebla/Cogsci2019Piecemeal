# Cogsci2019Piecemeal
Model fitting code for the Cogsci proceeding paper "A Piecemeal Processing Strategy Model for Causal-Based Categorization". All data included.

Contents:

models.R: all models used in the paper are here. For each model there is a cuadratic loss version (for model fitting) and a standard version (for prediction).

utlities.R: functions to loop over participants and estimate best-fitting parameter values and correlations of model predictions with data.

chain.compl.fitting.R: model fitting results for experiment 1, complete information condition.

chain.inc.fitting.R: model fitting results for experiment 1, incomplete information condition.

chain.compl.tcl.fitting.R: model fitting results for experiment 2 (traditional category labels), complete information condition.

chain.inc.tcl.fitting.R: model fitting results for experiment 2 (traditional category labels), incomplete information condition.
