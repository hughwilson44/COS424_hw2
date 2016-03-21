import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from scipy import spatial


def ReadData(filename):
<<<<<<< HEAD
    Data = pd.read_csv(filename, sep = '\t', header = None, names = None)
    colnames = ['chr', 'start', 'end', 'strand']
    colnames.extend([str(i) for i in range(1, Data.shape[1]-4)])
    colnames.append('450K')
    Data.columns = colnames
    return(Data)

def ReadDataHead(filename):
    '''read tab separated file into data frame with/without header'''
    Data = pd.read_csv(filename, sep = '\t', names = None)
    colnames = ['chr', 'start', 'end', 'strand']
    colnames.extend([str(i) for i in range(1, Data.shape[1]-4)])
    colnames.append('450K')
    Data.columns = colnames
    return(Data)
=======
	Data = pd.read_csv(filename, sep = '\t', header = None, names = None)
	colnames = ['chr', 'start', 'end', 'strand']
	colnames.extend([str(i) for i in range(1, Data.shape[1]-4)])
	colnames.append('450K')
	Data.columns = colnames
	return(Data)


# Calculating cosine similarity, ignoring NaNs
def CosineSim_nan(X,Y):
	X1 = []
	Y1 = []
	for i in range(0,len(X)):
		if not (np.isnan(X[i]) or np.isnan(Y[i])):
			X1.append(X[i])
			Y1.append(Y[i])
	return (1 - spatial.distance.cosine(X1,Y1))


# Calculating cosine similarities
def SimilarityMatrix(Data):
	SampleNum = Data.shape[1] - 5
	Similarities = np.zeros([SampleNum,SampleNum])

	print('calculating correlation')
	for i in range(0, SampleNum):
		print(i)
		for j in range(i, SampleNum):
			if i == j:
				Similarities[i,j] = 1
			else:
				Similarities[i,j] = CosineSim_nan(Data.iloc[:,i+4], Data.iloc[:,j+4])
		print(Similarities[i,0:6])

	for i in range(1, SampleNum):
		for j in range(0, i):
			if Similarities[i,j] == 0:
				Similarities[i,j] = Similarities[j,i]

	return Similarities


# Eliminating NaNs in training data by taking the weighted mean
# Data is a dataframe generated from ReadData
def EliminateNaN(Data, Similarities):
	SampleNum = Data.shape[1] - 5

	print('Eliminating NaNs')
	for i in range(0, SampleNum):
		print(i)
		NumNaNs = 0
		for r in range(0, Data.shape[0]):
			if np.isnan(Data.iloc[r,i+4]):
				Data.iloc[r,i+4] = 0
				temp = 0
				for j in range(0, SampleNum):
					if not (np.isnan(Data.iloc[r,j+4]) or i == j):
						temp = temp + Similarities[i,j]
				for j in range(0, SampleNum):
					if not (np.isnan(Data.iloc[r,j+4]) or i == j):
						Data.iloc[r,i+4] = Data.iloc[r,i+4] + (Similarities[i,j] / temp) * Data.iloc[r,j+4]
				NumNaNs = NumNaNs + 1
		print('in sample', str(i), ',', str(NumNaNs), 'NaNs replaced')

	return Data
	

# Eliminating NaNs in training data by taking the mean
# Data is a dataframe generated from ReadData
def EliminateNaN_mean(Data):
	SampleNum = Data.shape[1] - 5
	
	Values = Data.iloc[:, 4:(SampleNum+4)]

	print('calculating means')
	Means = Values.mean(axis=1, skipna=True)

	print('Replacing NaNs')
	for r in range(0,Data.shape[0]):
		for j in range(4,(SampleNum+4)):
			if np.isnan(Data.iloc[r,j]):
				Data.iloc[r,j] = Means[r]

	return Data
>>>>>>> 11d14bfe2c55d4944b6ff2b9bbd012b8e3659aec

