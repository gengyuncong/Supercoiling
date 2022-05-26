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
def get_elongation(S1, S2, elongation):
	S11 = S1[750:]
	S22 = S2[750:]
	if (S22[-1] < S11[0]+1):
		return elongation
	else:
		for j in range(S11[0]+1, S22[-1]+1):
			t1 = np.where(S11==j)[0][0]
			t2 = np.where(S22==j)[0][0]
			c = 1200.0/((t2-t1)*1.0);
			elongation = np.append(elongation, c)
		return elongation

def cal_SEM(data, N):
	mean_ = np.zeros(N);
	for i in range(N):
		mean_[i] = np.mean(np.random.choice(data,len(data)))
	return (np.nanmean(mean_), np.nanstd(mean_))

def get_mRNA(myfile,num):
	fp = h5py.File(myfile, "r");
	replicates=fp["/Simulations"].keys();
#	initiation = np.array([])
	elongation = np.array([])
	for i,replicate in enumerate(replicates):
		S1 = fp["/Simulations/%s/SpeciesCounts"%replicate][:, 13*num+2];
		S2 = fp["/Simulations/%s/SpeciesCounts"%replicate][:, 13*num+4];
#		initiation = get_initiation(S1, initiation);
		elongation = get_elongation(S1, S2, elongation);
	return (elongation)

def get_mRNA_tand(myfile,num):
	fp = h5py.File(myfile, "r");
	replicates=fp["/Simulations"].keys();
#	initiation1 = np.array([])
#	initiation2 = np.array([])
#	elongation1 = np.array([])
	elongation2 = np.array([])
	for i,replicate in enumerate(replicates):
#		S1_1 = fp["/Simulations/%s/SpeciesCounts"%replicate][:, 13*num+2];
		S1_2 = fp["/Simulations/%s/SpeciesCounts"%replicate][:, 13*num+3];
#		S2_1 = fp["/Simulations/%s/SpeciesCounts"%replicate][:, 13*num+4];
		S2_2 = fp["/Simulations/%s/SpeciesCounts"%replicate][:, 13*num+5];
#		elongation1 = get_elongation(S1_1, S2_1, elongation1)
		elongation2 = get_elongation(S1_2, S2_2, elongation2)
	return (elongation2)


##############################################################
import os
from os import listdir
from os.path import isfile, join

import random

lists = ['0.001', '0.005', '0.01', '0.02', '0.05', '0.1', '0.2'];
x_ = np.array([0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]); 
n = len(lists);

num = 180;

figname = 'Fig_6D_elong.pdf'
plt.rcParams["figure.figsize"] = (14,3.5)
labels = ['Convergent', 'Divergent' ,'Codirectional-1', 'Codirectional-2'];

j = 1;
for prefix in ['conv', 'div', 'tand', 'tand2']:
	k = 0;
	mean = np.zeros(n);
	std = np.zeros(n);
	for i in lists:
		no_domain = prefix+'.'+str(i)+'.20.sbml.lm';
		if prefix == "tand2" or prefix == "div":
			tmp = get_mRNA_tand(no_domain,num);
		else:
			tmp = get_mRNA(no_domain,num);
		(mean[k], std[k]) = cal_SEM(tmp, 100); 
#		mean[k] = np.nanmean(tmp);
#		std[k] = np.nanstd(tmp)/np.sqrt(len(tmp[~np.isnan(tmp)]));
		k = k+1;
	print(mean)
	print(std)
	subplot(1,4,j);
#	plt.axhline(0, color=new_colors[3], linewidth=2.0);
	plt.plot(x_, mean, 'o-', color=new_colors[0], linewidth=2.0)
	plt.fill_between(x_, mean - std, mean + std, color=new_colors[0], alpha=.25);
	plt.xlabel(r'$k_{max}$')
	plt.xscale('log')
#	plt.title(labels[j-1])
	plt.ylim((0, 50))
	if j == 1:
		plt.ylabel('Elongation rate (bp/s)')
	j = j+1;


plt.tight_layout()
plt.savefig(figname)
plt.close()


