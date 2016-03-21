import pandas as pd
import numpy as np
from ReadData import ReadData
# from ReadData import CosineSim_nan
# from ReadData import SimilarityMatrix
# from ReadData import EliminateNaN
from ReadData import EliminateNaN_mean


# ----------
# import data
# ----------
Datapath = '../../methylation_imputation/data/'
DataSample_full = Datapath + 'intersected_final_chr1_cutoff_20_sample_full.bed'
DataSample_partial = Datapath + 'intersected_final_chr1_cutoff_20_sample_partial.bed'
DataTrain = Datapath + 'intersected_final_chr1_cutoff_20_train.bed'

DataTrain = ReadData(DataTrain)
DataSample_full = ReadData(DataSample_full)
DataSample_partial = ReadData(DataSample_partial)

# Inds1 are columns indices in DataTrain without '26' and '27'
Inds1 = [x for x in range(0,25+4)]
Inds1.extend([x for x in range(27+4,DataTrain.shape[1])])

# Eliminating NaNs with 26 27
print('working on with 26 27')
Train = EliminateNaN_mean(DataTrain)
Train.to_csv('../result/Train_NaN_Meaned', sep = '\t', header = True, index = False)

# Eliminating NaNs without 26 27
print('working on without 26 27')
Train = EliminateNaN_mean(DataTrain.iloc[:,Inds1])
Train.to_csv('../result/Train_NaN_Meaned_without_2627', sep = '\t', header = True, index = False)


Train = pd.read_csv('../result/Train_NaN_Meaned', sep = '\t')
Train.set_index('start', drop=False, inplace=True, verify_integrity=True)
Train.to_csv('../result/test', sep = '\t', header = True, index = True, index_label = True)