
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
for N in [64,128,256,512]:
	plt.clf()
	fi=open('1.1(N='+str(N)+').txt','r')
	while ('seconds' not in fi.readline()):
		pass
	data=np.zeros((N,N))
	for i in range(N):
	    ip=fi.readline().split()
	    for j in range(N):
	        data[j][i]=float(ip[j])
	data=pd.DataFrame(data)
	sns.heatmap(data,xticklabels=[],yticklabels=[]).invert_yaxis()
	plt.savefig('1.1(N='+str(N)+').jpg')
for N in [64,128,256,512]:
	plt.clf()
	fi=open('1.2(N='+str(N)+').txt','r')
	while ('seconds' not in fi.readline()):
		pass
	data=np.zeros((N,N))
	for i in range(N):
	    ip=fi.readline().split()
	    for j in range(N):
	        data[j][i]=float(ip[j])
	data=pd.DataFrame(data)
	sns.heatmap(data,xticklabels=[],yticklabels=[]).invert_yaxis()
	plt.savefig('1.2(N='+str(N)+').jpg')