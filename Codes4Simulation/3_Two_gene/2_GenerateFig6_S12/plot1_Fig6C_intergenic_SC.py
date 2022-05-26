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

new_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']
##############################################################
def get_mRNA(myfile,num):
	fp = h5py.File(myfile, "r");
	replicates=fp["/Simulations"].keys();
	initiation1 = np.array([])
	initiation2 = np.array([])
	elongation1 = np.array([])
	elongation2 = np.array([])
	SC = np.array([]);
	for i,replicate in enumerate(replicates):
		count = fp["/Simulations/%s/SpeciesCounts"%replicate]
		RNAP1 = count[:, range(num, num+90)]+count[:, range(2*num,2*num+90)]+count[:, range(8*num, 8*num+90)]
		RNAP2 = count[:, range(num+90,2*num)]+count[:, range(2*num+90,3*num)]+count[:, range(8*num+90,9*num)]
		
		SC_all = np.array(count[750:, 3*num+90]);
		SC = np.append(SC, (np.mean(SC_all)-60.0)/60.0)
	return SC
#	return (np.mean(SC), np.std(SC))

def cal_SEM(data, N):
	mean_ = np.zeros(N);
	for i in range(N):
		mean_[i] = np.mean(np.random.choice(data,len(data)))
	return np.std(mean_)

##############################################################
import os
from os import listdir
from os.path import isfile, join

import random

lists = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2];
x_ = np.array([0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]); 
n = len(lists);

num = 180;

figname = 'Fig_6C_intergenic_SC.pdf'
plt.rcParams["figure.figsize"] = (14,3.5)
labels = ['Convergent', 'Divergent' ,'Codirectional-1', 'Codirectional-2'];

j = 1;
for prefix in ['conv', 'div', 'tand', 'tand2']:
	k = 0;
	conv = np.zeros(n);
	conv_err = np.zeros(n);
	for i in lists:
		no_domain = prefix+'.'+str(i)+'.20.sbml.lm';
		tmp = get_mRNA(no_domain,num);
		conv[k] = np.mean(tmp);
		conv_err[k] = cal_SEM(tmp, 100); #np.std(tmp);
		k = k+1;
	subplot(1,4,j)
	plt.plot(x_, conv, 'o-', color=new_colors[0], linewidth=2.0)
	plt.fill_between(x_, conv - conv_err,conv + conv_err, alpha=.25)
	plt.xlabel(r'$k_{max}$')
	plt.xscale('log')
	if j == 1:
		plt.ylabel('Intergenic SC')
		plt.ylim(0, 0.7)
	elif j == 2:
		plt.ylim(-0.7, 0)
	elif j == 3:
		plt.ylim(-0.35, 0.35)
	elif j == 4:
		plt.ylim(-0.35, 0.35)
#	plt.title(labels[j-1])
	j = j+1;

plt.tight_layout()
plt.savefig(figname)
plt.close()


