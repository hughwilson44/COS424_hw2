import pandas as pd

def ReadData(filename):
	Data = pd.read_csv(filename, sep = '\t', header = None, names = None)
	colnames = ['chr', 'start', 'end', 'strand']
	colnames.extend([str(i) for i in range(1, Data.shape[1]-4)])
	colnames.append('450K')
	Data.columns = colnames
	return(Data)
