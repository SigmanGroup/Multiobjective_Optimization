# Note: This code is currently writen exclusively for 2-parameter logistic regression searches
# Candidate logistic regression models are built for every single parameter and all unique combinations of 2 parameters that
# are below the colinearity cutoff.

import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression, LinearRegression
import logreg_stats as lr_stats
import itertools
import time

import multiprocessing
nproc = max([1,multiprocessing.cpu_count()-2])
from joblib import Parallel,delayed

class Model:
    def __init__(self,terms,X,y,reg):  #
        self.terms = terms                           # model terms (e.g. 'x1', 'x2', etc)
        self.n_terms = len(terms)                    # number of terms in model (int)
        self.formula = '1 + ' + ' + '.join(terms)    # mathematic expression for formula (str)
        self.model = reg.fit(X,y)                    # linear regression modeling object (can be called upon for sklearn functions)
        self.accuracy = self.model.score(X,y)        # model's accuracy
        self.mcfad_r2 = mcfad_par(terms, X, y)

def corrmap_par(t1,t2,data):                         #Checks the correlation of two parameters
    i,j = data[t1].values.reshape(-1,1),data[t2].values.reshape(-1,1)
    t1t2corr = LinearRegression().fit(i,j)
    score = t1t2corr.score(i,j)
    return(t1,t2,score)

def mcfad_par(terms, X, y):
    cand_mcfad_r2 = lr_stats.calc_mcfad(X, y)
    return(cand_mcfad_r2)

# step_par takes list of terms, data, response column name, regressor object
# returns terms, an sklearn Model object, accuracy score and response column name
def step_par(terms,data,response,reg):
    terms = tuple(terms)
    model = Model(terms,data.loc[:,terms],data[response],reg) # instantiates model object. Model trains on only terms provided!
    score = model.accuracy # Generates accuracy score
    mcfad_r2 = model.mcfad_r2 # Generates mcFadden R2 score
    return(terms,model,score,response, mcfad_r2)

# ForwardStep_py Arguments:
# data: pd.dataframe with all features and a response column
# response = the name of the column in the data dataframe that stores the reaction output (DDG, CD, etc...)
# reg = regressor object, and collinearity criteria (float from 0.0 to 1.0), colin_criteria = colinearity criteria (float from 0.0 to 1.0)
# ForwardStep_py Returns:
# results: Final dataframe with the following columns: 'Model' - ('x1', ), ('x1', 'x2')..., n_terms, and accuracy
# models: {('x1',), <MODEL OBJECT>, ('x3', 'x25'), <MODEL OBJECT>}
# scores_mcfad_r2: {('x1',), 0.428, ('x3', 'x25'), 0.551, ...} (values = mcfad_r2)
# sortedmodels: [('x5',), ('x1', 'x2'), ('x1', 'x3'), ('x1', 'x2', 'x5'), ... ]
def ForwardStep_py(data,response,reg=LogisticRegression(max_iter=50000),collin_criteria=0.5):
    start_time = time.time() # Time when function was called
    pool = Parallel(n_jobs=nproc,verbose=0) # Allows the job to be run in parallel on multiple processors
    features = list(data.columns) # it is advised to have the column titles as x1...xn
    features.remove(response)

    corrmap = data.corr() # pearson correlation coefficient R: -1 ... 1
    collin_criteria = np.sqrt(collin_criteria) # convert from R2 to R

    models,scores_accuracy,scores_mcfad_r2 = {},{},{}  # Creates 3 empty dictionaries that will store models, accuracy scores and mcfadden R2 scores

    for step in [1,2]:
        print("Step " + str(step))
        if step == 1: # if it is on the first step, creates a list of 1d tuples. one for each feature
            todo = [(feature,) for feature in features]
        if step == 2: # Gets all combinations of two parameters that are below the collin_criteria
            todo = sorted([(t1,t2) for (t1,t2) in itertools.combinations(features,step) if abs(corrmap.loc[t1,t2]) < collin_criteria])

        parall = pool(delayed(step_par)(terms,data,response,reg) for terms in todo) # This creates a list of outputs for the fuction called within pool.
        for results in parall:
            if len(results) == 0: # if step_par returns nothing, the tuple is skipped
                continue
            models[results[0]] = results[1] # populates the models dictionary (step 1): {('x1',): Model object} (step 2): {('x1', 'x6'): Model object}                                                            
            scores_accuracy[results[0]] = results[2] # populates the scores_accuracy dictionary (step 1): {('x1',): accuracy float} (step 2): {('x1', 'x6'): accuracy float}
            scores_mcfad_r2[results[0]] = results[4]

    candidates_accuracy = tuple(sorted(scores_accuracy,key=scores_accuracy.get,reverse=True))
    print('Finished 1 and 2 parameter models. Time taken (sec): %0.4f' %((time.time()-start_time)))
    sortedmodels = sorted(scores_accuracy,key=scores_accuracy.get,reverse=True)

    results_d = {
        'Model': sortedmodels,
        'n_terms': [models[terms].n_terms for terms in sortedmodels],
        'Accuracy': [models[terms].accuracy for terms in sortedmodels],
        'McFadden_R2': [models[terms].mcfad_r2 for terms in sortedmodels]
    }
    results = pd.DataFrame(results_d)
    print('Done. Time taken (minutes): %0.2f' %((time.time()-start_time)/60))
    return(results,models,scores_mcfad_r2,sortedmodels,candidates_accuracy)
