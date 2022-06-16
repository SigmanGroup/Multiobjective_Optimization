import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import LeaveOneOut, RepeatedKFold, cross_val_score 
from sklearn import metrics
from sklearn.metrics import f1_score, recall_score, precision_score, accuracy_score

# calc_mcfad only takes one X/y dataset
def calc_mcfad(X, y, reg=LogisticRegression()):
    lr = reg.fit(X,y)
    y_pred_prob = lr.predict_proba(X)[:,1] #List of predicted probabilities of output being "1" or above y-cut
    null_probability = np.count_nonzero(y)/len(y) # overall probability of being "1" in dataset

    # Calculate log likelihood and null log likelihood
    null_log_likelihood = 0
    for i in y:
        if i == 1:
            null_log_likelihood += np.log(null_probability)
        elif i == 0:
            null_log_likelihood += np.log(1-null_probability)
        else:
            print("ERROR!!! Values input into logistic regressor are not 1's and 0's.")

    log_likelihood = 0
    for i in range(len(y)):
        if y[i] == 1:
            #calculate the log likelihood where likelihood = probability
            log_likelihood_i = np.log(y_pred_prob[i])
            log_likelihood += log_likelihood_i
        elif y[i] == 0:
            #calculate the log likelihood of 1-probability
            log_likelihood_i = np.log(1-y_pred_prob[i])
            log_likelihood += log_likelihood_i
        else:
            print("ERROR Values input into logistic regressor are not 1's and 0's.")

    Mcfadden_R2 = 1 - log_likelihood/null_log_likelihood
    return(Mcfadden_R2)

# calc_mcfadden_R2 takes both training and test set data
def calc_mcfadden_R2(X_train, y_train, X_test, y_test):

    train_test_split = True
    if len(X_test) == 0:
        train_test_split = False
    lr = LogisticRegression().fit(X_train,y_train)
    y_pred_prob = lr.predict_proba(X_train)[:,1] #List of predicted probabilities of output being "1" or above y-cut
    if train_test_split:
        y_pred_prob_test = lr.predict_proba(X_test)[:,1]
    null_probability = np.count_nonzero(y_train)/len(y_train) # overall probability of being "1" in train dataset
    if train_test_split:
        null_probability_test = np.count_nonzero(y_test)/len(y_test) # overall probability of being "1" in test dataset

    # Calculate log likelihood and null log likelihood
    null_log_likelihood = 0 #The sum of the log likelyhoods of the null

    for i in y_train:
        if i == 1:
            null_log_likelihood += np.log(null_probability)
        elif i == 0:
            null_log_likelihood += np.log(1-null_probability)
        else:
            print("ERROR!!!")

    log_likelihood = 0 #The sum of the log likelyhoods of the null

    for i in range(len(y_train)):
        if y_train[i] == 1:
            #calculate the log likelihood where likelihood = probability
            log_likelihood_i = np.log(y_pred_prob[i])
            log_likelihood += log_likelihood_i
        elif y_train[i] == 0:
            #calculate the log likelihood of 1-probability
            log_likelihood_i = np.log(1-y_pred_prob[i])
            log_likelihood += log_likelihood_i
        else:
            print("ERROR")

    Mcfadden_R2_test = 'N/A'
    if train_test_split:
        null_log_likelihood_test = 0 #The sum of the log likelyhoods of the null

        for i in y_test:
            if i == 1:
                null_log_likelihood_test += np.log(null_probability_test)
            elif i == 0:
                null_log_likelihood_test += np.log(1-null_probability_test)
            else:
                print("ERROR!!!")

        log_likelihood_test = 0 #The sum of the log likelyhoods of the null

        for i in range(len(y_test)):
            if y_test[i] == 1:
                #calculate the log likelihood where likelihood = probability
                log_likelihood_i = np.log(y_pred_prob_test[i])
                log_likelihood_test += log_likelihood_i
            elif y_test[i] == 0:
                #calculate the log likelihood of 1-probability
                log_likelihood_i = np.log(1-y_pred_prob_test[i])
                log_likelihood_test += log_likelihood_i
            else:
                print("ERROR")
        Mcfadden_R2_test = 1 - log_likelihood_test/null_log_likelihood_test

    Mcfadden_R2 = 1 - log_likelihood/null_log_likelihood
    return(Mcfadden_R2, Mcfadden_R2_test)

def precision_recall_f1_score(X_train_sel, y_train):
    lr = LogisticRegression(max_iter=5000).fit(X_train_sel, y_train)
    y_pred = lr.predict(X_train_sel)
    precision = precision_score(y_train, y_pred, average='macro')
    f1 = f1_score(y_train, y_pred, average='macro')
    recall = recall_score(y_train, y_pred, average='macro')
    return(precision, recall, f1)

def test_accuracy_score(X_test_sel, y_test, X_train_sel, y_train):
    lr = LogisticRegression(max_iter=5000).fit(X_train_sel, y_train)
    y_pred = lr.predict(X_test_sel)
    test_accuracy = accuracy_score(y_test, y_pred)
    return(test_accuracy)

def kfold_logreg(X_train, y_train, k=5):
    # prepare the cross-validation procedure
    cv = RepeatedKFold(n_splits=k, n_repeats=3, random_state=1)
    # create model
    model = LogisticRegression(max_iter=5000)
    # evaluate model
    scores = cross_val_score(model, X_train, y_train, scoring='accuracy', cv=cv, n_jobs=-1)
    # report performance
    kfold_score, kfold_stdev = np.mean(scores), np.std(scores)
    #print('Accuracy: %.3f (%.3f)' % (kfold_score, kfold_stdev))
    return(kfold_score, kfold_stdev)
