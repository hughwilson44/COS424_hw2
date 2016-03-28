
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
import time

# Enter your directory here
#os.chdir('/Users/hugh/Google Drive/Hugh/PhD/princeton/COS424/assignments/

def lassolars(data,target):
    '''Fit a lasso model to the data and target using the LARS
    algorithm'''
    clf = linear_model.LassoLars( alpha = 0.1, fit_intercept = False )
    #clf = linear_model.LassoLarsCV(cv=5,fit_intercept=False)
    model = clf.fit(data,target)
    return(model)

def ridgeregression(data,target):
    '''Fit a model to the data using ridge regression'''
    #clf = linear_model.RidgeCV(alphas=(1e-5,1e-4,1e-3,1e-2,0.1,1,10),fit_intercept=False, \
    #                           store_cv_values=True)
    clf = linear_model.Ridge(alpha=0.1,fit_intercept=False)
    model = clf.fit(data,target)
    return(model)

def elasticnet(data,target):
    '''Fit a model to the data using elastic nets'''
    #clf = linear_model.ElasticNetCV(l1_ratio=0.75,fit_intercept=False)
    clf = linear_model.ElasticNet(alpha=0.1, l1_ratio=0.75,fit_intercept=False)
    model = clf.fit(data,target)
    return(model)

def coefficients2file(coefficients,filename):
    '''write the coefficients to file'''
    np.savetxt(filename+'_coefs.txt',coefficients,delimiter='\t')

def predictions2file(results,filename):
    pandas.DataFrame.to_csv(results,'./'+filename+'_predictions.txt', \
                            sep='\t', index=False)


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
    el_results = Dat_sfull

    # METHOD 1
    # Train the estimator across specimens at fixed genome locus
    # prepare the data for passing to the regressors
    datind = map(str, range(1,Dat_train.shape[1]-4))
    Dat_Xtrn1 = Dat_train[Dat_spar['1'].notnull()][datind]
    Dat_Ytrn1 = Dat_spar['1'][Dat_spar['1'].notnull()]

    # establish a start time
    start_time = time.time()

    # lasso LARS
    lasso_model = lassolars(Dat_Xtrn1,Dat_Ytrn1)
    lasso_results['M1_pred'] = lasso_model.predict(Dat_train[datind])
    lasso_M1_coefs = lasso_model.coef_
    print("--- lasso %s seconds ---" % (time.time() - start_time))
    running_time = time.time()

    # ridge regression
    rr_model = ridgeregression(Dat_Xtrn1,Dat_Ytrn1)
    rr_results['M1_pred'] = rr_model.predict(Dat_train[datind])
    rr_M1_coefs = rr_model.coef_
    print("--- ridge %s seconds ---" % (time.time() - running_time))
    running_time = time.time()

    # Elastic nets
    el_model = elasticnet(Dat_Xtrn1,Dat_Ytrn1)
    el_results['M1_pred'] = el_model.predict(Dat_train[datind])
    el_M1_coefs = el_model.coef_
    print("--- elastic %s seconds ---" % (time.time() - running_time))
    running_time = time.time()

    # Print the alpha values
#    print('lasso',lasso_model.alpha_)
#    print('ridge',rr_model.alpha_)
#    print('elastic',el_model.alpha_)

    #-------------------------------------------------------------------

    # METHOD 2a
    # Train the estimator down a specimen for a given specimen
    # the data remains the same, but we loop over the target
    Dat_Xtrn2 = Dat_train[Dat_spar['1'].notnull()]
    Dat_train.set_index('start', drop=False, inplace=True, verify_integrity=True)
    Dat_Xpred2 = Dat_spar[Dat_spar['1'].notnull()]
    Dat_Xpred2.set_index('start', drop=False, inplace=True, verify_integrity=True)

    # Read in the top predictor values
    TopPrd = pandas.read_csv('./TopPredictors',sep='\t')
    TopPrd.set_index('start', drop=False, inplace=True, verify_integrity=True)

    # initialise holder arrays
    lasso_M2a_pred = np.zeros(TopPrd.shape[0])
    rr_M2a_pred = np.zeros(TopPrd.shape[0])
    el_M2a_pred = np.zeros(TopPrd.shape[0])

    # initialise a counter
    count = 0
    running_time = time.time()

    # for target in list(Dat_spar.index[:5]):
    for target in list(TopPrd.loc[:,'start']):
        Dat_Ytrn2 = Dat_train.loc[target,datind]
        Dat_Xtrn2fs = Dat_train.loc[TopPrd.loc[target,'predictor1':]][datind].transpose()
        Dat_Xpred2fs = Dat_Xpred2.loc[TopPrd.loc[target,'predictor1':],'1']
        # lasso LARS
        lasso_model = lassolars(Dat_Xtrn2fs,Dat_Ytrn2)
        lasso_M2a_pred[count] = \
            lasso_model.predict(Dat_Xpred2fs.reshape(1,-1))    #.reshape(1,-1))
        # rr
        rr_model = ridgeregression(Dat_Xtrn2fs,Dat_Ytrn2)
        rr_M2a_pred[count] = \
            rr_model.predict(Dat_Xpred2fs.reshape(1,-1))    #.reshape(1,-1))
        # el
        el_model = elasticnet(Dat_Xtrn2fs,Dat_Ytrn2)
        el_M2a_pred[count] = \
            el_model.predict(Dat_Xpred2fs.reshape(1,-1))    #.reshape(1,-1))

        count = count + 1

    lasso_results_multi_method = lasso_results.set_index('start',drop=False, verify_integrity=True)
    lasso_results_multi_method = lasso_results_multi_method.loc[TopPrd.loc[:,'start']]
    lasso_results_multi_method['M2a_pred'] = lasso_M2a_pred

    rr_results_multi_method = rr_results.set_index('start',drop=False, verify_integrity=True)
    rr_results_multi_method = rr_results_multi_method.loc[TopPrd.loc[:,'start']]
    rr_results_multi_method['M2a_pred'] = rr_M2a_pred

    el_results_multi_method = el_results.set_index('start',drop=False, verify_integrity=True)
    el_results_multi_method = el_results_multi_method.loc[TopPrd.loc[:,'start']]
    el_results_multi_method['M2a_pred'] = el_M2a_pred
    print("--- method2 %s seconds ---" % (time.time() - running_time))

    #-------------------------------------------------------------------

    # print the predictions and coefficients to file
    predictions2file(lasso_results,'lasso_M1')
    predictions2file(lasso_results_multi_method, 'lasso_methods')
    predictions2file(rr_results,'rr_M1')
    predictions2file(rr_results_multi_method, 'rr_methods')
    predictions2file(el_results,'el_M1')
    predictions2file(el_results_multi_method, 'el_methods')
    coefficients2file(lasso_M1_coefs,'lasso_M1')
    coefficients2file(rr_M1_coefs,'rr_M1')
    coefficients2file(el_M1_coefs,'el_M1')
#    coefficients2file(lasso_M2a_coefs,'lasso_M2a')

if __name__ == '__main__':
    main()
#    boston = load_boston()
#    data = boston.data
#    target = boston.target
#    print(lassolars(data,target))

