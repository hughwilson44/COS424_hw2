import pandas as pd
import numpy as np
from ReadData import ReadData
from ConvertResult import ConvertResultImputation

# ----------
# import original data
# ----------
Datapath = '../../methylation_imputation/data/'
DataSample_full = Datapath + 'intersected_final_chr1_cutoff_20_sample_full.bed'
DataSample_partial = Datapath + 'intersected_final_chr1_cutoff_20_sample_partial.bed'
DataTrain = Datapath + 'intersected_final_chr1_cutoff_20_train.bed'

# Training data
# use either one of the following three:
# DataTrain = ReadData(DataTrain)
# DataTrain = pd.read_csv('../result/Train_NaN_Meaned', sep = '\t')
DataTrain = pd.read_csv('../result/Train_NaN_Meaned_without_2627', sep = '\t')

# Sample data
DataSample_full = ReadData(DataSample_full)
DataSample_partial = ReadData(DataSample_partial)


# ----------
# Format Data so they are indexed by start position
# ----------
DataTrain.set_index('start', drop=False, inplace=True, verify_integrity=True)
DataSample_full.set_index('start', drop=False, inplace=True, verify_integrity=True)
DataSample_partial.set_index('start', drop=False, inplace=True, verify_integrity=True)


# ----------
# Read the Imputation result
# ----------
Resultpath = '../result/'
method = 'lasso_M1'
filename = Resultpath + 'Imputation_smry_' + method

# This script reads data already processed by Result_Error_CV.py,
# SI the data table contains conlumn 'true' 'CV' and 'error'
Distances = pd.read_csv(filename, sep='\t')
# Result = ConvertResult(Result)
Distances.set_index('start', drop=False, inplace=True, verify_integrity=True)

# The column names for the predictors
Predictors = Distances.columns.tolist()[:]

# ----------
# get the distance
# ----------
for r in Distances.index.tolist():
	for pred in Predictors:
		Distances.loc[r,pred] = abs(int(r) - int(pred))		# the distance between predictor and site
		if abs(int(r) - int(pred)) == 0:
			print("site at same position as predictor, means traning data included.\n Result_Error_CV.py doesn't work well in: Delete rows that were used as training data")
			print(r,pred)

# ----------
# Get the min max distance
# ----------
Mins = []
Maxs = []
for r in Distances.index.tolist():
	Mins.append(np.min(Distances.loc[r,Predictors]))
	Maxs.append(np.max(Distances.loc[r,Predictors]))
Distances.insert(Distances.shape[1], 'disMin', Mins)
Distances.insert(Distances.shape[1], 'disMax', Maxs)

# ----------
# Save to file
# ----------
Distances.to_csv('../result/' + 'DistanceError_smry_' + method, index=False, header=True, sep='\t')

