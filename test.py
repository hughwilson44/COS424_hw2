import pandas as pd
from ReadData import ReadData

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

# ----------
# ----------
	
