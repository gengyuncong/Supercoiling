from __future__ import division
import matplotlib
matplotlib.use('Agg')
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import h5py
import math
import sys

fontSize=16
matplotlib.rcParams.update({"axes.formatter.limits": (-4,4), "svg.fonttype" : "none", 'pdf.fonttype':42,'font.family':'sans-serif','font.sans-serif':'Helvetica','font.size': fontSize, "axes.titlesize": fontSize, "xtick.labelsize": fontSize, "ytick.labelsize": fontSize,'text.usetex':True,'text.latex.preamble':[r'\usepackage{sansmath}',r'\sansmath']})
plotStyles={"markersize":8,"markeredgewidth":1.0,"linewidth":3.0}
stepStyles={"markersize":20,"markeredgewidth":3.0,"linewidth":3.0,"where":"post"}
barStyles={"width":0.65, "linewidth":0, "align":"center", "color":"steelblue"}

new_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']

##############################################################
'''
0-1:	DNA
1-2:	RNAP
2-3:	RNAP_tmp
3-4:	Twist
4-5:	Gyrase_unbind
5-6:	Gyrase_bind
6-7:	TopoI_unbind
7-8:	TopoI_bind
8-9:	RNAP_stall
9-10:	mRNA	

'''
##############################################################
num=70;

def get_mRNA(filename):
	fp = h5py.File(filename, "r");
	replicates=fp["/Simulations"].keys();
	mrna=np.array([], dtype=float)
	mymean = 0
	fano = 0
	for i,replicate in enumerate(replicates):
		counts = fp["/Simulations/%s/SpeciesCounts"%replicate]
#		c = fp["/Simulations/%s/SpeciesCounts"%replicate][:,10*num+54]; 
#		c = fp["/Simulations/%s/SpeciesCounts"%replicate][:,range(10*num, 10*num+80)]
#		c = np.sum(counts[-1,range(10*num, 11*num)])/40.0; 
		c = counts[-1, 10*num+54];
		mrna=np.append(mrna, c);
	if (len(mrna)>0):
		mymean = np.mean(mrna)
		fano=np.var(mrna)/np.mean(mrna)
	return (mymean, fano)


lists1 = ['0.001', '0.005', '0.01', '0.02', '0.05', '0.08', '0.1', '0.15', '0.2', '0.3', '0.4', '0.5'];
lists2 = ['0.001', '0.005', '0.01', '0.02', '0.05', '0.08', '0.1', '0.15', '0.2'];

mean = []
fano = []

mean_d1 = []
fano_d1 = []

mean_d2 = []
fano_d2 = []

mean_c = []
fano_c = []

def push(a1,a2,t1,t2):
	a1 = np.append(a1, t1)
	a2 = np.append(a2, t2)
	return (a1, a2)

def re_order(a1, a2):
	index = np.argsort(a1)
	a1 = a1[index]
	a2 = a2[index]
	return (a1, a2)

for i in lists1:
	filename = 'domain.'+str(i)+'.sbml.lm'
	(t1, t2) = get_mRNA(filename);
	(mean_d1, fano_d1) = push(mean_d1, fano_d1, t1, t2);

	filename = 'domain_v2.'+str(i)+'.sbml.lm'	
	(t1, t2) = get_mRNA(filename);
	(mean_d2, fano_d2) = push(mean_d2, fano_d2, t1, t2);

	filename = 'closed.'+str(i)+'.sbml.lm'
	(t1, t2) = get_mRNA(filename);
	(mean_c, fano_c) = push(mean_c, fano_c, t1, t2);


for i in lists2:
	filename = 'normal.'+str(i)+'.sbml.lm'
	(t1, t2) = get_mRNA(filename);
	(mean, fano) = push(mean, fano, t1, t2);	
	
(mean_d1, fano_d1) = re_order(mean_d1, fano_d1)
(mean_d2, fano_d2) = re_order(mean_d2, fano_d2)
(mean_c, fano_c) = re_order(mean_c, fano_c)
(mean, fano) = re_order(mean, fano)


figname = 'S10_fano_all4.pdf'
plt.rcParams["figure.figsize"] = (5,4)
plt.axhline(y=1, color=new_colors[7], linewidth=2)
plt.plot(mean, fano, 'o-', color=new_colors[0], linewidth=2, label='Open')
plt.plot(mean_d2, fano_d2, 'o-', color=new_colors[1], linewidth=2, label='Dynamic1')
plt.plot(mean_d1, fano_d1, 'o-', color=new_colors[2], linewidth=2, label='Dynamic2')
plt.plot(mean_c, fano_c, 'o-', color=new_colors[3], linewidth=2, label='Closed')
plt.ylabel('Fano factor')
plt.xlabel('Mean mRNA copy number')
plt.xscale('log')
plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig(figname)
plt.close()


