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
#matplotlib.rcParams.update({"axes.formatter.limits": (-4,4), "svg.fonttype" : "none", 'font.family':'MathJax_SansSerif', 'font.size': fontSize, "axes.titlesize": fontSize, "xtick.labelsize": fontSize, "ytick.labelsize": fontSize,'text.usetex':False,'text.latex.preamble':[r'\usepackage{sansmath}',r'\sansmath']})
#fontSize=20
plotStyles={"markersize":8,"markeredgewidth":1.0,"linewidth":3.0}
stepStyles={"markersize":20,"markeredgewidth":3.0,"linewidth":3.0,"where":"post"}
barStyles={"width":0.65, "linewidth":0, "align":"center", "color":"steelblue"}

##############################################################
import os
from os import listdir
from os.path import isfile, join

num=149;

def get_statistics(myfile):
	fp = h5py.File(myfile, "r");
	replicates=fp["/Simulations"].keys();
	counts=np.zeros((len(replicates),), dtype=int)
	mRNA=np.zeros(1001)
	initiation=np.array([])
	for i,replicate in enumerate(replicates):
		count = fp["/Simulations/%s/SpeciesCounts"%replicate]
		c1 = np.sum(count[:,range(10*num, 10*num+80)], axis=1)
		mRNA = mRNA+c1/50.0
		tt = 750;
		S1 = count[:,13*num+1];
		S11 = S1[tt:] - S1[tt];
		initiation=np.append(initiation, S11[-1]/(1.0*len(S11)))
		
	return (mRNA/(1.0*len(replicates)), np.mean(initiation))

new_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']

k = 1;
lists3 = [0.001, 0.005, 0.01, 0.02, 0.05, 0.08, 0.1, 0.15, 0.2]; 
plt.rcParams["figure.figsize"] = (15,10)
figname = 'S4_mRNA_over_time.pdf'
for ini_rate in lists3:
	hnormal = 'sens.'+str(ini_rate)+'.50.1.sbml.lm';
	(mRNA, ini_empirical) = get_statistics(hnormal);
	subplot(3,3,k);
	plt.title('Initiation rate='+str(round(ini_empirical, 3))+'/s')
	step(range(0,len(mRNA)), mRNA, color=new_colors[0], **stepStyles)
	if k > 6:
		xlabel('Time (s)',fontsize=15); 
	if k in [1, 4, 7]:
		ylabel('mean mRNA copy number',fontsize=15); 
	k = k+1; 
plt.tight_layout()
plt.savefig(figname)
plt.close()


