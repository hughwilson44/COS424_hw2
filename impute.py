#!/usr/bin/python

"""
Script to produce a model for imputing missing values
in the sample partial file. Will return the imputed values
and the most significant predictors with their coefficients
"""

# Import statements and data ReadDataing taken from the precept code
# import statements
#import math
import numpy as np
from sklearn import linear_model
from sklearn.datasets import load_boston
#from sklearn import tree
#from IPython.display import Image
#from sklearn.externals.six import StringIO
#import pydot
#from collections import OrderedDict
#import matplotlib.pyplot as plt
import pandas
#import os
import ReadData as RD
# %matplotlib inline

# Enter your directory here
#os.chdir('/Users/hugh/Google Drive/Hugh/PhD/princeton/COS424/assignments/2')

def lassolars(data,target):
    '''Fit a lasso model to the data and target using the LARS
    algorithm'''
    clf = linear_model.LassoLarsCV(cv=5,fit_intercept=False)
    model = clf.fit(data,target)
    return(model)

def ridgeregression(data,target):
    '''Fit a model to the data using ridge regression'''
    clf = linear_model.RidgeCV(alphas=(0.1,1,10),fit_intercept=False, \
                               store_cv_values=True)
    model = clf.fit(data,target)
    return(model)

def coefficients2file(coefficients,filename):
    '''write the coefficients to file'''
    np.savetxt(filename+'_coefs.txt',coefficients,delimiter='\t')

def predictions2file(results,filename):
    pandas.dataframe.to_csv(results,'./'+filename+'_predictions.txt', \
                            sep='\t', index=false)

def main():
    # Read in the data files
    Datapath = './data/'
    Path_sfull = Datapath + 'intersected_final_chr1_cutoff_20_sample_full.bed'
    Path_spar = Datapath + 'intersected_final_chr1_cutoff_20_sample_partial.bed'
    Path_train = Datapath + 'Train_NaN_Meaned_without_2627'

    Dat_train = RD.ReadDataHead(Path_train)
    Dat_spar = RD.ReadData(Path_spar)
    Dat_sfull = RD.ReadData(Path_sfull)

    # Produce a copy of the full file to which predictions can be added
    lasso_results = Dat_sfull
    rr_results = Dat_sfull

    # METHOD 1
    # Train the estimator across specimens at fixed genome locus
    # prepare the data for passing to the regressors
    datind = map(str, range(1,Dat_train.shape[1]-4))
    Dat_Xtrn1 = Dat_train[Dat_spar['1'].notnull()][datind]
    Dat_Ytrn1 = Dat_spar['1'][Dat_spar['1'].notnull()]

    # lasso LARS
    lasso_model = lassolars(Dat_Xtrn1,Dat_Ytrn1)
    lasso_results['M1_pred'] = lasso_model.predict(Dat_train[datind])

    # ridge regression
    rr_model = ridgeregression(Dat_Xtrn1,Dat_Ytrn1)
    rr_results['M1_pred'] = rr_model.predict(Dat_train[datind])

    #-------------------------------------------------------------------

    # METHOD 2a
    # Train the estimator down a specimen for a given specimen
    # the data remains the same, but we loop over the target
    Dat_Xtrn2 = Dat_Xtrn1.transpose()
    Dat_Xpred2 = Dat_Ytrn1

    # initialise holder arrays
    lasso_M2a_coefs = np.zeros([Dat_spar.shape[0],np.shape(datind)[0]])
    lasso_M2a_pred = np.zeros(Dat_spar.shape[0])

    for target in list(Dat_spar.index):
        Dat_Ytrn2 = Dat_train.loc[target,datind]

        # lasso LARS
        lasso_model = lassolars(Dat_Xtrn2,Dat_Ytrn2)
        lasso_M2a_pred[target] = lasso_model.predict(Dat_Xpred2)
        lasso_M2a_coefs[target] = lasso_model.coef_

    lasso_results['M2a_pred'] = lasso_M2a_pred

    #-------------------------------------------------------------------

    # print the predictions and coefficients to file
    predictions2file(lasso_results,'lasso')
    coefficients2file(lasso_M1_coefs,'lasso_M1')
    coefficients2file(lasso_M2a_coefs,'lasso_M2a')

    if __name__ == '__main__':
#    boston = load_boston()
#    data = boston.data
#    target = boston.target
#    print(lassolars(data,target))
        main()
