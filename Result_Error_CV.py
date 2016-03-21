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
Resultpath = '../result/raw/'
method = 'lasso_M1'
# method = 'rr_M1'
filename = Resultpath + method + '_predictions.txt'
print(method)
Result = pd.read_csv(filename, sep='\t')
# Result = ConvertResult(Result)
Result.set_index('start', drop=False, inplace=True, verify_integrity=True)


# ----------
# Add true value
# ----------
# Result is in the form: ['start', 'strand', 'predictor1', 'predictor2', ..., 'predictorN', 'predicted', 'true'] indexed by 'start'
Result.insert(Result.shape[1], 'true', Result.loc[:,'1'])


# ----------
# Delete rows that were used as training data
# ----------
# delete rows that have value in Sample_partial
Result = Result[np.isnan(DataSample_partial['1'])]


# ----------
# Delete NaN rows in Result (because Sample_full still contains ~3000 NaNs)
# ----------
Result = Result[np.isfinite(Result['true'])]


# -----------
# Calculate CV (coefficient of variance for each site)
# -----------
# CV is calculated from the original data, not the data with NaNs replaced, NaNs are ignored when calculating CV
temp = []
for r in Result.index.tolist():
	temp.append(DataTrain.loc[r,'1':'33'].std(skipna=True) / DataTrain.loc[r,'1':'33'].mean(skipna=True))
Result.insert(Result.shape[1], 'CV', temp)


# ----------
# Calculate Errors
# ----------
# Error is defined as abs(predict - true)
temp = []
for r in Result.index.tolist():
	# should be either one of the following. If Result contains true value, can use second one. first one is always right
	temp.append(abs(Result.loc[r,'predicted'] - DataSample_full.loc[r,'1']))
	# temp.append((Result.loc[r,'predict'] - Result.loc[r,'true']) ** 2)
Result.insert(Result.shape[1], 'error', temp)

# When write to file, do not write index, because index will just be wrote as the first column
# Even use index_label=False, which makes index first column with no column name.
# when read in R, the index will still be read as first column with colum name 'True'
# Result.to_csv(Resultpath + '/Imputation', index=True, index_label=False, header=True, sep='\t')
Result.to_csv('../result/' + 'Imputation_smry_' + method, index=False, header=True, sep='\t')