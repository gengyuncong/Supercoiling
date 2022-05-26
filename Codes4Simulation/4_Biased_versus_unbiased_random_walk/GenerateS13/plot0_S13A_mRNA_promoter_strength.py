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
num=149;

def get_mRNA(filename):
	fp = h5py.File(filename, "r");
	replicates=fp["/Simulations"].keys();
	mrna=np.array([], dtype=float)
	mymean = 0
	fano = 0
	for i,replicate in enumerate(replicates):
		count = fp["/Simulations/%s/SpeciesCounts"%replicate]
		for t in [-1]:
			c1 = np.sum(count[t,range(10*num, 10*num+80)])/50.0;
#			c1 = count[t,10*num+80];
			mrna=np.append(mrna, c1);
	if (len(mrna)>0):
		mymean = np.mean(mrna)
		fano=np.var(mrna)/np.mean(mrna)
	return (mymean, fano)

def push(a1,a2,t1,t2):
	a1 = np.append(a1, t1)
	a2 = np.append(a2, t2)
	return (a1, a2)

def re_order(a1, a2):
	index = np.argsort(a1)
	a1 = a1[index]
	a2 = a2[index]
	return (a1, a2)


lists1 = ['0.001', '0.005', '0.01', '0.02', '0.05', '0.08', '0.1', '0.15', '0.2'];

figname = 'Fig_7A_mRNA.pdf'

plt.rcParams["figure.figsize"] = (12,3.5)
n_k = 1;
D_dict = {"0.002":"5", "0.02":"50", "0.2":"500"};

for D in ['0.002', '0.02', '0.2']:
	mean_d1 = []
	fano_d1 = []
	mean_d2 = []
	fano_d2 = []
	
	prefix = {"explicit":"expl_sens_"+D_dict[D]+".", "biased":"D"+D_dict[D]+"/sens."}; 
	if D == '0.02':
		prefix = {"explicit":"expl_sens.", "biased":"D"+D_dict[D]+"/sens."};
	suffix = {"explicit":".sbml.lm", "biased":"."+D_dict[D]+".1.sbml.lm"}
	for ini_rate in lists1:
		filename = prefix["explicit"]+str(ini_rate)+suffix["explicit"];
		(t1, t2) = get_mRNA(filename);
		(mean_d1, fano_d1) = push(mean_d1, fano_d1, t1, t2);
		
		filename = prefix["biased"]+str(ini_rate)+suffix["biased"];
		(t1, t2) = get_mRNA(filename);
		(mean_d2, fano_d2) = push(mean_d2, fano_d2, t1, t2);
		
		(mean_d1, fano_d1) = re_order(mean_d1, fano_d1)
		(mean_d2, fano_d2) = re_order(mean_d2, fano_d2)

	x_ = np.array([0.001, 0.005, 0.01, 0.02, 0.05, 0.08, 0.1, 0.15, 0.2])
	plt.subplot(1,3,n_k)
	plt.title('D='+D+' '+r'$um^2/s$')
#	plt.axhline(y=1, color=new_colors[7], linewidth=2)
	plt.plot(x_, mean_d2, 'o-', color=new_colors[0], linewidth=2, label='explicit')
	plt.plot(x_, mean_d1, 'o-', color=new_colors[1], linewidth=2, label='biased')
	if n_k == 1:
		plt.ylabel('Mean mRNA copy number')
	plt.xlabel('$k_{max}$')
	plt.xscale('log')	
	plt.ylim((0.03, 25)); 
	plt.xscale('log')
	plt.yscale('log')
	if n_k == 2:
		plt.legend(loc='lower right')
	n_k = n_k + 1;

plt.tight_layout()
plt.savefig(figname)
plt.close()


