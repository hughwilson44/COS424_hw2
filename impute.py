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
    clf = linear_model.LassoLars(alpha=.1,fit_intercept=False)
    model = clf.fit(data,target)
    return(model)

def model2predictions(model,data,sample_full,filename):
    '''Use a fitted regression model to produce a table with:
        the genome start locus as row numbers, then: the genome start
        locus; the strand the site is on; and the imputed
        value for the sample.'''
    sample_full['predicted'] = model.predict(data)
    pandas.DataFrame.to_csv(sample_full,'./' + filename, sep='\t', index=False)


def main():
    # Read in the data files
    Datapath = './data/'
    Path_sfull = Datapath + 'intersected_final_chr1_cutoff_20_sample_full.bed'
    Path_spar = Datapath + 'intersected_final_chr1_cutoff_20_sample_partial.bed'
    Path_train = Datapath + 'Train_NaN_Meaned_without_2627'

    Dat_train = RD.ReadDataHead(Path_train)
    Dat_spar = RD.ReadData(Path_spar)
    Dat_sfull = RD.ReadData(Path_sfull)

    # Call the lassolars regressor and print the result to file
    datind = map(str, range(1,Dat_train.shape[1]-4))
    Dat_Xtrn = Dat_train[Dat_spar['1'].notnull()][datind]
    Dat_Ytrn = Dat_spar['1'][Dat_spar['1'].notnull()]
    lasso_model = lassolars(Dat_Xtrn,Dat_Ytrn)
    model2predictions(lasso_model,Dat_train[datind],Dat_sfull, \
                      'lasso_predictions.csv')


if __name__ == '__main__':
    boston = load_boston()
    data = boston.data
    target = boston.target
    print(lassolars(data,target))
    main()
