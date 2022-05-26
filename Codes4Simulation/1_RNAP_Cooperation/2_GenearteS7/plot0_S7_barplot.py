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
import os
from os import listdir
from os.path import isfile, join

num=149;

def get_elongation_rate(myfile):
	fp = h5py.File(myfile, "r");
	replicates=fp["/Simulations"].keys();
	counts=np.array([])
	elongation=np.array([])
	for i,replicate in enumerate(replicates):
		S1 = fp["/Simulations/%s/SpeciesCounts"%replicate][:,13*num+1];
		S2 = fp["/Simulations/%s/SpeciesCounts"%replicate][:,13*num+2];
		if (S1[-1]<0 or S2[-1]<0):
			next;
		else:
			for ii in range(1,S2[-1]+1):
				t1 = np.where(S1==ii)[0][0]
				t2 = np.where(S2==ii)[0][0]
				if (t1>750):
					elongation = np.append(elongation, int(3000/(1.0*(t2-t1))))
	return (np.mean(elongation), np.std(elongation))


def get_elongation_rate_early(myfile):
	fp = h5py.File(myfile, "r");
	replicates=fp["/Simulations"].keys();
	counts=np.array([])
	elongation=np.array([])
	for i,replicate in enumerate(replicates):
		S1 = fp["/Simulations/%s/SpeciesCounts"%replicate][:,13*num+1];
		S2 = fp["/Simulations/%s/SpeciesCounts"%replicate][:,13*num+2];
		if (S1[-1]<0 or S2[-1]<0):
			next;
		else:
			for ii in range(1,S2[-1]+1):
				t1 = np.where(S1==ii)[0][0]
				t2 = np.where(S2==ii)[0][0]
				if (t1>100):
					elongation = np.append(elongation, int(3000/(1.0*(t2-t1))))
	return (np.mean(elongation), np.std(elongation))


def get_elongation_rate_PI(myfile):
	fp = h5py.File(myfile, "r");
	replicates=fp["/Simulations"].keys();
	counts=np.array([])
	elongation=np.array([])
	for i,replicate in enumerate(replicates):
		S1 = fp["/Simulations/%s/SpeciesCounts"%replicate][:,13*num+1];
		S2 = fp["/Simulations/%s/SpeciesCounts"%replicate][:,13*num+2];
		S4 = fp["/Simulations/%s/SpeciesCounts"%replicate][:,13*num+4];
		tt = np.where(S4==1)[0];
		if (S1[-1]<0 or S2[-1]<0 or len(tt)==0):
			next;
		else:
			for ii in range(1,S2[-1]+1):
				t1 = np.where(S1==ii)[0][0]
				t2 = np.where(S2==ii)[0][0]
				if (t1>750 and tt[0]>750):
					elongation = np.append(elongation, int(3000/(1.0*(t2-t1))))
	return (np.mean(elongation), np.std(elongation))



new_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']

hinac='anta.0.1_1.sbml.lm';
hnormal='sens.0.1.50.1.sbml.lm';
linac='anta.0.001_1.sbml.lm';
lnormal='sens.0.001.50.1.sbml.lm';

mean = np.zeros(8)
std = np.zeros(8)

(mean[0], std[0]) = get_elongation_rate_early(hnormal)
(mean[1], std[1]) = get_elongation_rate_early(hinac)
(mean[2], std[2]) = get_elongation_rate(lnormal)
(mean[3], std[3]) = get_elongation_rate_PI(linac)

figname='S7_barplot.pdf'
plt.rcParams["figure.figsize"] = (5,4)
for i in range(4):
	my_color = new_colors[4];
	if i%2 == 0:
		my_color = new_colors[7];
	if i == 0:
		plt.bar(i, mean[i], yerr=std[i], width=0.8, color=my_color, label='On')
	elif i == 1:
		plt.bar(i, mean[i], yerr=std[i], width=0.8, color=my_color, label='Off')
	else:
		plt.bar(i, mean[i], yerr=std[i], width=0.8, color=my_color);
print(mean)
xticks([0.5, 2.5], ['Group', 'Solo'])
plt.legend()
plt.ylabel('Elongation rate (bp/s)')
plt.tight_layout()
plt.savefig(figname)
plt.close()

