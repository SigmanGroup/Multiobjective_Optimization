# Multiobjective_Optimization

This repository contains the scripts necessary for:
- Logistic regression
- Parameter acquisition for bisphosphine ligands
- Chemical space visualization for bisphosphine ligands

## Logistic Regression 

The logistic regression script is designed to work with two excel spreadsheets: one with ligand descriptors and one reaction data. For examples, 'Bisphosphine_descriptors' and 'Bisphosphine_rxn_data'.

Notes: 
- The Jupyter Notebook is designed to search for the best one-parameter and two-parameter logistic regression models. This version of the script does not search for 3+ parameter logistic regression models. 
- Notebook Sections 1-3 prepare the script and read in all data.
- Notebook Section 4 runs univariate logistic regression model search.
- Notebook Section 5-6 runs bivariate logistic regression model search.
- Notebook Section 7 is used for virtual screening of virtual ligands.
- It should be noted that the script only tests 2-term models for which the terms are not colinear with eachother (per a user-defined colinearity criteria). 
- Code for 2-parameter model search is contained in Logistic_Regression.py
- Code for validation statistics in contained in logreg_stats.py.

## Chemical Space
This Jupyter notebook is designed to be run with one excel spreadsheet which contains our bisphosphine ligand parameters, with all ligands treated as though they posses C2v symmetry. Therefore, the parameters used in this representation of chemical space do not account for absolute stereochemistry of the ligands.

Notes: 
- Initialization: imports the necessary packages.
- Load data: Reads in the data from the associated excel spreadsheet
- Run PCA analysis: Performs principal componant analysis on the selected parameters
- Plot a subset of ligands in curated chemical space: Here you can define lists of ligands that you wish to plot in this chemical space
- Plot 2D Chemical Space on Choice of Two Dimensions: Plots the various 2D chemical space maps (Steric/Electronic, Geometric/ElectroniC and Geometric/Steric).
- Neighbour Analysis: Use this section to find the distance between test set ligands and their nearest training set ligand.
