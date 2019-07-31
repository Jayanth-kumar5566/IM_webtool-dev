import pandas

def ind(x):
	return set(x.index)

def re_index(x):
	'''x is a list containing the dataframes to be merged '''
	ind_list=map(ind,x)
	ind_f=set.itersection(ind_list)
	
