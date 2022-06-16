# Multiobjective_Optimization

This repository contains the scripts necessary for:
- Logistic regression
- Parameter acquisition for bisphosphine ligands
- Chemical space visualization for bisphosphine ligands

## Logistic Regression

The logistic regression script is designed to work with two excel spreadsheets: one with ligand descriptors and one reaction data. For examples, 'Bisphosphine_descriptors' and 'Bisphosphine_rxn_data'. 

The Jupyter Notebook is designed to search for the best one-parameter and two-parameter logistic regression models. This version of the script does not search for 3+ parameter logistic regression models. 

Notebook Sections 1-3 prepare the script and read in all data.

Notebook Section 4 runs univariate logistic regression model search.

Notebook Section 5-6 runs bivariate logistic regression model search.

Notebook Section 7 is used for virtual screening of virtual ligands.

It should be noted that the script only tests 2-term models for which the terms are not colinear with eachother (per a user-defined colinearity criteria). 

Logistic_Regression_notebook.ipynb: code for 2-parameter model search is contained in Logistic_Regression.py. Code for validation statistics in contained in logreg_stats.py.

DELETE KENNARDSTONEALGORIGTHM.PY

