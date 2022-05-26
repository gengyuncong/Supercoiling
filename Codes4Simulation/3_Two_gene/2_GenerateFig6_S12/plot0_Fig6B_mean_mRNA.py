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
	counts = np.zeros((len(replicates),2), dtype=int)
#	initiation = np.array([])
	for i,replicate in enumerate(replicates):
		counts[i,0] = np.sum(fp["/Simulations/%s/SpeciesCounts"%replicate][-1, range(10*num,10*num+90)])/(20.0);#1,2,..,80
		counts[i,1] = np.sum(fp["/Simulations/%s/SpeciesCounts"%replicate][-1, range(10*num+90,11*num)])/(20.0); #81,82,...,160
	return counts

def cal_SEM(data, N):
	mean_ = np.zeros(N);
	for i in range(N):
		mean_[i] = np.mean(np.random.choice(data,len(data)))
	return (np.mean(mean_), np.std(mean_))

##############################################################
import os
from os import listdir
from os.path import isfile, join

import random

lists = ['0.001', '0.005', '0.01', '0.02', '0.05', '0.1', '0.2'];
num = 180;

n = len(lists);

prefix = 'conv';
k = 0;
conv = np.zeros(n);
conv_err = np.zeros(n);
for i in lists:
	no_domain = prefix+'.'+str(i)+'.20.sbml.lm';
	tmp = get_mRNA(no_domain,num);
	data = tmp[:,0];
	(conv[k], conv_err[k]) = cal_SEM(data, 100);
	k = k+1;

prefix = 'div';
k = 0;
div = np.zeros(n); 
div_err = np.zeros(n);
for i in lists:
	no_domain = prefix+'.'+str(i)+'.20.sbml.lm';
	tmp = get_mRNA(no_domain,num);
	data = tmp[:,1];
	(div[k], div_err[k]) = cal_SEM(data, 100);
	k = k+1;

prefix = 'tand';
k = 0;
tand1 = np.zeros(n); tand1_err = np.zeros(n); 
for i in lists:
	no_domain = prefix+'.'+str(i)+'.20.sbml.lm';
	tmp = get_mRNA(no_domain,num);
	data = tmp[:,0];
	(tand1[k], tand1_err[k]) = cal_SEM(data, 100);
	k = k+1;

prefix = 'tand2';
k = 0;
tand2 = np.zeros(n); tand2_err = np.zeros(n);
for i in lists:
	no_domain = prefix+'.'+str(i)+'.20.sbml.lm';
	tmp = get_mRNA(no_domain,num);
	data = tmp[:,1];
	(tand2[k], tand2_err[k]) = cal_SEM(data, 100);
	k = k+1;


x_ = np.array([0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]); 
#x_ = np.array([2, 8, 20, 40, 80]);
#x_ = x_*60;

figname = 'Fig_6B.pdf'
plt.rcParams["figure.figsize"] = (14,3.5)
subplot(1,4,1)
plt.plot(x_, conv, 'o-', color=new_colors[0], linewidth=2.0)
plt.fill_between(x_, conv - conv_err,conv + conv_err, alpha=.25)
plt.xlabel(r'$k_{max}$')
plt.ylabel('mRNA copy number')
#plt.title('Convergent')
plt.xscale('log')
plt.ylim((0,6))

subplot(1,4,2)
plt.plot(x_, div, 'o-', color=new_colors[0], linewidth=2.0)
plt.fill_between(x_, div - div_err, div + div_err, alpha=.25)
plt.xlabel(r'$k_{max}$')
#plt.title('Divergent')
plt.xscale('log')
plt.ylim((0,6))
#'''
subplot(1,4,3)
plt.plot(x_, tand1, 'o-', color=new_colors[0], linewidth=2.0)
plt.fill_between(x_, tand1-tand1_err, tand1+tand1_err, alpha=.25)
plt.xlabel(r'$k_{max}$')
#plt.title('Codirectional-1')
plt.xscale('log')
plt.ylim((0,6))

subplot(1,4,4)
plt.plot(x_, tand2, 'o-', color=new_colors[0], linewidth=2.0)
plt.fill_between(x_, tand2-tand2_err, tand2+tand2_err, alpha=.25)
#plt.legend()
plt.xlabel(r'$k_{max}$')
#plt.title('Codirectional-2')
plt.xscale('log')
plt.ylim((0,6))
#'''
plt.tight_layout()
plt.savefig(figname)
plt.close()


