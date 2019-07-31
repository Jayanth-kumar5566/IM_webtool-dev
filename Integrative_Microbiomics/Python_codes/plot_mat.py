#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#import seaborn as sns
import pandas

def bplot(M,lab):
	lab.sort_values('labels',inplace=True)
	M=M.reindex(lab.index,axis='index')
	M=M.reindex(lab.index,axis='columns')
	return M
	#ax=sns.heatmap(M,xticklabels=False,yticklabels=False,cmap=plt.get_cmap("Blues"))
	#return plt.show()
	
