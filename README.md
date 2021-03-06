# Multiobjective_Optimization

This repository contains the scripts necessary for:
- Logistic regression
- Parameter acquisition for bisphosphine ligands
- Chemical space visualization for bisphosphine ligands

## Logistic Regression 

The logistic regression script is designed to work with two excel spreadsheets: one with ligand descriptors and one reaction data. For examples, see 'Bisphosphine_descriptors' and 'Bisphosphine_rxn_data'.

Notes: 
- The Jupyter Notebook (Logistic_Regression_notebook.ipynb) is designed to search for the best one-parameter and two-parameter logistic regression models. This version of the script does not search for 3+ parameter logistic regression models. 
- It should be noted that the script only tests 2-term models for which the terms are not collinear with each other (per a user-defined collinearity criteria). 
- Code for 2-parameter model search is contained in Logistic_Regression.py
- Code for validation statistics in contained in logreg_stats.py.

## Chemical Space
This Jupyter notebook is designed to be run with one excel spreadsheet which contains bisphosphine ligand parameters, with all ligands treated as though they possess C2v symmetry. Therefore, the parameters used in this representation of chemical space do not account for absolute stereochemistry of the ligands.

Notes: 
- Initialization: imports the necessary packages.
- Load data: Reads in the data from the associated excel spreadsheet
- Run PCA analysis: Performs principal component analysis on the selected parameters
- Plot a subset of ligands in curated chemical space: Here you can define lists of ligands that you wish to plot in this chemical space
- Plot 2D Chemical Space on Choice of Two Dimensions: Plots the various 2D chemical space maps (Steric/Electronic, Geometric/Electronic, and Geometric/Steric).
- Neighbour Analysis: Use this section to find the distance between test set ligands and their nearest training set ligand.

## GetSum Parameters

- The GetSum parameter acquisition script is designed to pull parameters from bisphosphine ligands. It is critical that the folder system is organized in the same manner as the example included in this repository. In the same folder that the 'Funky_conformers_v3' Jupyter Notebook is run from there must be a folder titled ???DFT_files_atom_nums???. Nested in this folder should be ???Atom_numbers.xlsx??? ??? a spreadsheet containing the atom numbers for all of the ligands and another folder ???Ligand_Calcs???. Nested within the ligand calcs folder is an individual folder for each ligand. The ligand folders must have the following naming format. 
???Ligand ID???_???ligand name???. If there are spaces in the name, these should be converted to underscores.
(e.g. SS Et DuPhos (ID 10) would have the folder name ???10_SS_Et_DuPhos???).

- Each ligand folder must contain the all of the final DFT output files (ending in ???.log???) for all of the conformers for the ligand. Each DFT-output file should have the following naming convention:
- 'Ligand name _ conformer number(1-5) _ job-type.log'
- A max of 5 conformers are computed at the DFT level for each ligand. Therefore, the script is designed for up to 5 conformers per ligand.
- Job-type is blank for optimizations, ???SPE??? for single-point calculations with the PdCl2, and ???SPE_NoPd??? for the single-point calculations without PdCl2. 
(e.g. DFT-output files for conformer #1 for (S,S) Et DuPhos would be ???SS_Et_DuPhos_1.log???, ???SS_Et_DuPhos_1_SPE.log???, and ???SS_Et_DuPhos_1_SPE_NoPd.log???)
- NOTE: See the included 'DFT_files_atom_nums' folder for an example of the necessary filesystem.

Jupyter Notebook (Funky_conformers_v3.ipynb) Notes:
Once the atom numbers excel spreadsheet and the folders containing the log files are set up, the Jupyter Notebook requires minor user input. The user must specify local file paths, specifications of the atom numbers spreadsheet, and desired names for the resultant parameter spreadsheets that are produced (see Jupyter Notebook for more details). 

The Jupyter Notebook produces an excel spreadsheet with two sheets. One sheet treats all ligands according to their symmetry while the second simplifies parameters and treats all ligands as being C2v symmetric (see manuscript and supplementary information for more details).

Atom_numbers.xlsx Notes:
- Every ligand that will be parameterized must have a row in the atom numbers spreadsheet for both the PdCl2 complex and the complex without PdCl2
- In addition to inputing the atom numbers (see SI for more details) the following columns must be filled out
- Symmetry: Point-group of the ligand (C1, C2, C2v, or Cs)
- Renumber: True for C1 and C2 ligands that do not have the Walphos or Josiphos backbone. False for all other ligands
- Chiral: True if the ligand is chiral, otherwise False.
- log_name: name of the single-point file for the lowest energy conformer (e.g. 'SS_DuPhos_1_SPE')
- log_name (NoPd): name of the single-point file for the lowest energy conformer without PdCl2 (e.g. 'SS_DuPhos_1_SPE_NoPd')
- opt_log_name: name of the optimization file for the lowest energy conformer (e.g. 'SS_DuPhos_1')

## Citation

https://kjelljorner.github.io/morfeus/

## Contributors
- Chemical_Space: code by Lucy van Dijk and Jordan Dotson
- GetSum_Parameters: code by Jordan Dotson and Lucy van Dijk. Portions of code were adapted from a script by Tobias Gensch.
- Logistic Regression: code by Jordan Dotson with contributions from Lucy van Dijk. Portions of code were adpated from a script by Cian Kingston
