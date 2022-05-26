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
	counts1 = np.zeros(2001);
#	counts2 = np.zeros(2001);
	for i,replicate in enumerate(replicates):
		count = fp["/Simulations/%s/SpeciesCounts"%replicate];
		counts1 = counts1 + np.sum(count[:, range(10*num,10*num+90)], axis=1)/(20.0);#1,2,..,80
#		counts2 = counts2 + np.sum(count[:, range(10*num+90,11*num)], axis=1)/(20.0);
	return counts1/len(replicates)

def get_mRNA2(myfile,num):
	fp = h5py.File(myfile, "r");
	replicates=fp["/Simulations"].keys();
#	counts1 = np.zeros(2001);
	counts2 = np.zeros(2001);
	for i,replicate in enumerate(replicates):
		count = fp["/Simulations/%s/SpeciesCounts"%replicate];
#		counts1 = counts1 + np.sum(count[:, range(10*num,10*num+90)], axis=1)/(20.0);#1,2,..,80
		counts2 = counts2 + np.sum(count[:, range(10*num+90,11*num)], axis=1)/(20.0);
	return counts2/len(replicates)

##############################################################
import os
from os import listdir
from os.path import isfile, join

import random

lists = ['0.001', '0.005', '0.01', '0.02', '0.05', '0.1', '0.2'];
num = 180;

n = len(lists);

title_keys = {'conv':'Convergent', 'div': 'Divergent', 'tand':'Codirectional-1', 'tand2':'Codirectional-2'}; 
plt.rcParams["figure.figsize"] = (12,20)
figname = 'S10_mRNA_over_time.pdf';
k = 1;
for prefix in ['conv', 'div', 'tand', 'tand2']:
	plt.subplot(4,1,k);
	plt.title(title_keys[prefix])
	kj = 0;
	for i in lists:
		no_domain = prefix+'.'+str(i)+'.20.sbml.lm';
		if prefix == "tand2" or prefix == "div":
			tmp = get_mRNA2(no_domain,num);
		else:
			tmp = get_mRNA(no_domain,num);
		xdata = np.arange(0,2001,1);
		plt.plot(xdata, tmp, color=new_colors[kj], label=str(i));
		kj = kj+1;
	plt.xlabel('Time (s)')
	plt.ylabel('mean mRNA copy number')
	plt.legend(title=r'$k_{max}$')
	k = k+1;
plt.tight_layout()
plt.savefig(figname)
plt.close()

