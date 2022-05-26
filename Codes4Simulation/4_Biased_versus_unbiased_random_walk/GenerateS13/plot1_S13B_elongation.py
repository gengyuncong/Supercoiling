from __future__ import division
import matplotlib
matplotlib.use('Agg')
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import h5py
import math
import sys
import os
import re
import getopt
from scipy.optimize import curve_fit
from scipy.misc import factorial
import scipy.stats as stats  


fontSize=16
matplotlib.rcParams.update({"axes.formatter.limits": (-4,4), "svg.fonttype" : "none", 'pdf.fonttype':42,'font.family':'sans-serif','font.sans-serif':'Helvetica','font.size': fontSize, "axes.titlesize": fontSize, "xtick.labelsize": fontSize, "ytick.labelsize": fontSize,'text.usetex':True,'text.latex.preamble':[r'\usepackage{sansmath}',r'\sansmath']})
plotStyles={"markersize":8,"markeredgewidth":1.0,"linewidth":3.0}
stepStyles={"markersize":20,"markeredgewidth":3.0,"linewidth":3.0,"where":"post"}
barStyles={"width":0.65, "linewidth":0, "align":"center", "color":"steelblue"}

##############################################################
import os
from os import listdir
from os.path import isfile, join

num=149;

def get_elongation_rate(myfile):
	fp = h5py.File(myfile, "r");
	replicates=fp["/Simulations"].keys();
	initiation=np.array([])
	elongation=np.array([])
	tt = 750;
	for i,replicate in enumerate(replicates):
		S1 = fp["/Simulations/%s/SpeciesCounts"%replicate][:,13*num+1];
		S2 = fp["/Simulations/%s/SpeciesCounts"%replicate][:,13*num+2];
		S11 = S1[tt:] - S1[tt];
		initiation=np.append(initiation, S11[-1]/(1.0*len(S11)))
		if (S1[-1]<0 or S2[-1]<0):
			next;
		else:
			for ii in range(1,S2[-1]+1):
				t1 = np.where(S1==ii)[0][0]
				t2 = np.where(S2==ii)[0][0]
				if (t1>tt):
					elongation = np.append(elongation, int(3000/(1.0*(t2-t1))))
	return (initiation, elongation)

new_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']

lists3 = ['0.001', '0.005', '0.01', '0.02', '0.05', '0.08', '0.1', '0.15', '0.2'];

figname = 'Fig_7B_elongation_initiation.pdf'

plt.rcParams["figure.figsize"] = (12,3.5)
#ax1 = plt.subplot(1,1,1)
n = len(lists3);
n_k = 1;
D_dict = {"0.002":"5", "0.02":"50", "0.2":"500"};

for D in ['0.002', '0.02', '0.2']:
	j = 0;
	prefix = {"explicit":"expl_sens_"+D_dict[D]+".", "biased":"D"+D_dict[D]+"/sens."}; 
	if D == '0.02':
		prefix = {"explicit":"expl_sens.", "biased":"D"+D_dict[D]+"/sens."};
	suffix = {"explicit":".sbml.lm", "biased":"."+D_dict[D]+".1.sbml.lm"}
	for rw in ['explicit', 'biased']:
		initiation_list = np.zeros(n); 
		elongation_mean = np.zeros(n);
		elongation_sem = np.zeros(n); 
		k = 0;
		for ini_rate in lists3:
			filename = prefix[rw] + ini_rate + suffix[rw];
			(initiation, elongation) = get_elongation_rate(filename); 
			initiation_list[k] = np.mean(initiation);
			elongation_mean[k] = np.mean(elongation);
			elongation_sem[k] = np.std(elongation)/np.sqrt(len(elongation));
			k = k+1;
		print(initiation_list)
		print(elongation_mean)
		subplot(1,3,n_k)
		plt.title('D='+D+' '+r'$um^2/s$')
		plt.plot(initiation_list, elongation_mean,  'o-', color=new_colors[j], linewidth=2, label=rw, alpha = 1)
#		plt.fill_between(initiation_list, elongation_mean-elongation_sem, elongation_mean+elongation_sem, color=new_colors[j], alpha=.25)
		if n_k == 1:
			plt.ylabel('Elongation rate (bp/s)')
		if n_k == 2:
			plt.legend(loc='lower right')
		plt.xlabel('Initiation rate (/s)')
		plt.ylim((0,60))
		j = j+1; 
	n_k = n_k+1;

plt.tight_layout()
plt.savefig(figname)
plt.close()


