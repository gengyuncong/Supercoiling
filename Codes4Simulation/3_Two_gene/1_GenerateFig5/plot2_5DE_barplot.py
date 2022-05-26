from __future__ import division
import matplotlib
matplotlib.use('Agg')
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import h5py
import math
import os
import re


fontSize=16
matplotlib.rcParams.update({"axes.formatter.limits": (-4,4), "svg.fonttype" : "none", 'pdf.fonttype':42,'font.family':'sans-serif','font.sans-serif':'Helvetica','font.size': fontSize, "axes.titlesize": fontSize, "xtick.labelsize": fontSize, "ytick.labelsize": fontSize,'text.usetex':True,'text.latex.preamble':[r'\usepackage{sansmath}',r'\sansmath']})
plotStyles={"markersize":8,"markeredgewidth":1.0,"linewidth":3.0}
stepStyles={"markersize":20,"markeredgewidth":3.0,"linewidth":3.0,"where":"post"}
barStyles={"width":0.65, "linewidth":0, "align":"center", "color":"steelblue"}

new_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']
##############################################################
##############################################################
def get_initiation(S1, initiation):
	S11 = S1[750:]
	initiation = np.append(initiation, (S11[-1] - S11[0])/(1.0*len(S11)))
	return initiation

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

def get_mRNA(myfile,num):
	fp = h5py.File(myfile, "r");
	replicates=fp["/Simulations"].keys();
	initiation = np.array([])
	elongation = np.array([])
	for i,replicate in enumerate(replicates):
		S1 = fp["/Simulations/%s/SpeciesCounts"%replicate][:, 13*num+2];
		S2 = fp["/Simulations/%s/SpeciesCounts"%replicate][:, 13*num+4];
		initiation = get_initiation(S1, initiation);
		elongation = get_elongation(S1, S2, elongation);
	return (initiation, elongation)

def get_mRNA_tand(myfile,num):
	fp = h5py.File(myfile, "r");
	replicates=fp["/Simulations"].keys();
	initiation1 = np.array([])
	initiation2 = np.array([])
	elongation1 = np.array([])
	elongation2 = np.array([])
	for i,replicate in enumerate(replicates):
		S1_1 = fp["/Simulations/%s/SpeciesCounts"%replicate][:, 13*num+2];
		S1_2 = fp["/Simulations/%s/SpeciesCounts"%replicate][:, 13*num+3];
		S2_1 = fp["/Simulations/%s/SpeciesCounts"%replicate][:, 13*num+4];
		S2_2 = fp["/Simulations/%s/SpeciesCounts"%replicate][:, 13*num+5];
		initiation1 = get_initiation(S1_1, initiation1)
		initiation2 = get_initiation(S1_2, initiation2)
		elongation1 = get_elongation(S1_1, S2_1, elongation1)
		elongation2 = get_elongation(S1_2, S2_2, elongation2)
	return (initiation1, initiation2, elongation1, elongation2)

##############################################################
import os
from os import listdir
from os.path import isfile, join

import random

num = 180;
i = '20';

prefix = 'conv';
no_domain = prefix+'.0.05.'+str(i)+'.sbml.lm';
no_domain_ori = '../simul4fig6_v2/' + prefix+'.0.05.'+str(i)+'.sbml.lm';
(ini1, elong1) = get_mRNA(no_domain,num);
(ini2, elong2) = get_mRNA(no_domain_ori, num);

prefix = 'div';
no_domain = prefix+'.0.05.'+str(i)+'.sbml.lm';
no_domain_ori = '../simul4fig6_v2/' + prefix+'.0.05.'+str(i)+'.sbml.lm';
(ini3, elong3) = get_mRNA(no_domain,num);
(ini4, elong4) = get_mRNA(no_domain_ori, num);

#'''
prefix = 'tand';
no_domain = prefix+'.0.05.'+str(i)+'.sbml.lm';
no_domain_ori = '../simul4fig6_v2/' + prefix+'.0.05.'+str(i)+'.sbml.lm';
(ini5, ini6, elong5, elong6) = get_mRNA_tand(no_domain,num);
(ini5_, ini6_, elong5_, elong6_) = get_mRNA_tand(no_domain_ori,num);
#'''


figname = 'Fig_5D_barplot.pdf'
plt.rcParams["figure.figsize"] = (14,3.5)
subplot(1,2,1)
plt.ylabel('Initiation rate (/s)')
x_SC = np.array([0, 2, 4, 6]);
y_SC = np.array([np.mean(ini1), np.mean(ini3), np.mean(ini5), np.mean(ini6)])
y_SC_err = np.array([np.std(ini1), np.std(ini3), np.std(ini5), np.std(ini6)])
x_noSC = np.array([1, 3, 5, 7]); 
y_noSC = np.array([np.mean(ini2), np.mean(ini4), np.mean(ini5_), np.mean(ini6_)])
y_noSC_err = np.array([np.std(ini2), np.std(ini4), np.std(ini5_), np.std(ini6_)])
print(y_SC)
print(y_SC_err)
print(y_noSC)
print(y_noSC_err)
plt.bar(x_SC, y_SC, yerr=y_SC_err,width=0.8, color=new_colors[0],  alpha=0.4, label='W/O SC');
plt.bar(x_noSC, y_noSC, yerr=y_noSC_err, width=0.8, color=new_colors[1],  alpha=0.4, label='W/ SC');
plt.xticks([1, 3, 5, 7], ['Convergent','Divergent','Co-(up)','Co-(down)'])
plt.legend()

subplot(1,2,2)
plt.ylabel('Elongation rate (bp/s)')
x_SC = np.array([0, 2, 4, 6]);
y_SC = np.array([np.mean(elong1), np.mean(elong3), np.mean(elong5), np.mean(elong6)])
y_SC_err = np.array([np.std(elong1), np.std(elong3), np.std(elong5), np.std(elong6)])
x_noSC = np.array([1, 3, 5, 7]); 
y_noSC = np.array([np.mean(elong2), np.mean(elong4), np.mean(elong5_), np.mean(elong6_)])
y_noSC_err = np.array([np.std(elong2), np.std(elong4), np.std(elong5_), np.std(elong6_)])
print(y_SC)
print(y_SC_err)
print(y_noSC)
print(y_noSC_err)
plt.bar(x_SC, y_SC, yerr=y_SC_err, width=0.8, color=new_colors[0],  alpha=0.4, label='W/O SC');
plt.bar(x_noSC, y_noSC, yerr=y_noSC_err, width=0.8, color=new_colors[1],  alpha=0.4, label='W/ SC');
plt.xticks([1, 3, 5, 7], ['Convergent','Divergent','Co-(up)','Co-(down)'])
#plt.legend()


plt.tight_layout()
plt.savefig(figname)
plt.close()


