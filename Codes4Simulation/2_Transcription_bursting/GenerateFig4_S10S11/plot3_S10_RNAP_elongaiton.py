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

#np.set_printoptions(threshold=sys.maxsize)

fontSize=16
matplotlib.rcParams.update({"axes.formatter.limits": (-4,4), "svg.fonttype" : "none", 'pdf.fonttype':42,'font.family':'sans-serif','font.sans-serif':'Helvetica','font.size': fontSize, "axes.titlesize": fontSize, "xtick.labelsize": fontSize, "ytick.labelsize": fontSize,'text.usetex':True,'text.latex.preamble':[r'\usepackage{sansmath}',r'\sansmath']})
plotStyles={"markersize":8,"markeredgewidth":1.0,"linewidth":3.0}
stepStyles={"markersize":20,"markeredgewidth":3.0,"linewidth":3.0,"where":"post"}
barStyles={"width":0.65, "linewidth":0, "align":"center", "color":"steelblue"}

##############################################################
import os
from os import listdir
from os.path import isfile, join

num=70;
def get_pdf(mRNA, binwidth):
	(bins,edges) = np.histogram(mRNA,np.arange(-0.5*binwidth,np.max(mRNA)+1.5*binwidth,binwidth))
	centers = (edges[:-1]+edges[1:])/2
	xdata=centers
	ydata=bins/sum(bins)
	pdf = np.zeros((len(bins),2),dtype=float)
	pdf[:,0]=centers
	pdf[:,1] = bins.astype(double)/(np.sum(bins)*(centers[1]-centers[0]))
	return pdf

def get_all_statistics(myfile):
	fp = h5py.File(myfile, "r");
	replicates=fp["/Simulations"].keys();
	initiation=np.array([])
	elongation=np.array([])
	mrna=np.array([])
	density=np.array([])
	tt = 750;
	for i,replicate in enumerate(replicates):
		count = fp["/Simulations/%s/SpeciesCounts"%replicate]
		S1 = count[:,13*num+1];
		S2 = count[:,13*num+2];
		RNAP = count[750:, range(num, 2*num)]+count[750:, range(2*num,3*num)]+count[750:, range(8*num, 9*num)]
		S11 = S1[tt:] - S1[tt];
		initiation=np.append(initiation, S11[-1]/(1.0*len(S11)))
		#cal elongation rate
		if (S1[-1]<0 or S2[-1]<0):
			next;
		else:
			for ii in range(1,S2[-1]+1):
				t1 = np.where(S1==ii)[0][0]
				t2 = np.where(S2==ii)[0][0]
				if (t1>tt):
					elongation = np.append(elongation, int(2400/(1.0*(t2-t1))))

		c = count[-1, 10*num+54];
		mrna=np.append(mrna, c);

		RNAP = np.sum(RNAP, axis=1)
		time_RNAP = np.where(RNAP>0)[0]
		if len(time_RNAP) == 0:
			next;
		else:
			density=np.append(density, RNAP[time_RNAP]); 
	return (get_pdf(initiation,0.005), get_pdf(elongation,2), get_pdf(density,1))

new_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']

figname = 'S10_RNAP_number.pdf'; 
plt.rcParams["figure.figsize"] = (14,4)

filename = 'domain.0.1.sbml.lm'
(initiation, elongation, density) = get_all_statistics(filename)

filename = 'normal.0.1.sbml.lm'
(n_initiation, n_elongation, n_density) = get_all_statistics(filename)

subplot(1,3,1)
plt.bar(initiation[:,0], initiation[:,1], width=0.8*0.005, alpha=0.25, label='Looped')
plt.bar(n_initiation[:,0], n_initiation[:,1], width=0.8*0.005, alpha=0.25, label='Unlooped')
plt.xlabel('Initiation rate (/s)')
plt.ylabel('PDF')
plt.legend()

subplot(1,3,2)
plt.bar(elongation[:,0], elongation[:,1], width=0.8*2, alpha=0.25, label='Looped')
plt.bar(n_elongation[:,0], n_elongation[:,1], width=0.8*2, alpha=0.25, label='Unlooped')
plt.xlabel('Elongation rate (/s)')
plt.ylabel('PDF')
plt.legend()

subplot(1,3,3)
plt.bar(density[:,0], density[:,1], width=0.8, alpha=0.25, label='Looped')
plt.bar(n_density[:,0], n_density[:,1], width=0.8, alpha=0.25, label='Unlooped')
plt.xlabel('Number of co-transcribing RNAP')
plt.ylabel('PDF')
plt.legend()

plt.tight_layout()
plt.savefig(figname)
plt.close()

