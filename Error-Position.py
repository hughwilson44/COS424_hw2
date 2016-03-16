import pandas as pd
import numpy as np
from ReadData import ReadData
from ConvertResult import ConvertResultImputation

# ----------
# ----------
Resultpath = '../result/'
Result = ConvertResult(Result)

# ----------
# import data
# ----------
Datapath = '../../methylation_imputation/data/'
DataTrain = Datapath + 'intersected_final_chr1_cutoff_20_train.bed'
DataTrain = ReadData(DataTrain)

temp = []
for i in Result.shape[0]:
	r = DataTrain[DataTrain['start'] == Result.iloc[i,'start']].index.tolist()	
	temp.append(DataTrain.loc[r,'1':'33'].std(skipna=True) / DataTrain.loc[r,'1':'33'].mean(skipna=True))	# coefficient of variance
Result.insert(Result.shape[1], 'CV', temp)

Result.to_csv(Resultpath + '/Imputation', index=None, sep='\t')