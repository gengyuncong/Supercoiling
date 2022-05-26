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
binwidth = 0.1;

def get_pdf(mRNA):
	(bins,edges) = np.histogram(mRNA,np.arange(-0.5*binwidth-1.0, 1.0+1.5*binwidth,binwidth))
	centers = (edges[:-1]+edges[1:])/2
	xdata=centers
	ydata=bins/sum(bins)
	pdf = np.zeros((len(bins),2),dtype=float)
	pdf[:,0]=centers
	pdf[:,1] = bins.astype(double)/(np.sum(bins)*(centers[1]-centers[0]))
	return pdf

def get_mRNA(myfile,num):
	fp = h5py.File(myfile, "r");
	replicates=fp["/Simulations"].keys();
	SC = np.array([]);
	for i,replicate in enumerate(replicates):
		count = fp["/Simulations/%s/SpeciesCounts"%replicate]
		RNAP1 = count[:, range(num, num+90)]+count[:, range(2*num,2*num+90)]+count[:, range(8*num, 8*num+90)]
		RNAP2 = count[:, range(num+90,2*num)]+count[:, range(2*num+90,3*num)]+count[:, range(8*num+90,9*num)]
		SC_all = np.array(count[:, 3*num+90]);
		counts1 = np.sum(RNAP1[750:,:], axis=1)
		counts2 = np.sum(RNAP2[750:,:], axis=1)
		index1 = np.where(counts1>0)[0]
		index2 = np.where(counts2>0)[0]
		numerator = len(np.intersect1d(index1,index2))
		denominator = len(np.union1d(index1,index2))
		if (denominator == 0):
			next;
		else:
			SC_sub = SC_all[denominator];
			SC = np.append(SC, (np.mean(SC_sub)-60.0)/60.0)
	pdf = get_pdf(SC);
	return pdf

##############################################################
import os
from os import listdir
from os.path import isfile, join

import random

num = 180;
i = '20';

prefix = 'conv';
no_domain = prefix+'.0.05.'+str(i)+'.sbml.lm';
no_domain_ori = '../simul4fig6_v2/'+prefix+'.0.05.'+str(i)+'.sbml.lm';
pdf1 = get_mRNA(no_domain,num);
pdf2 = get_mRNA(no_domain_ori, num);

prefix = 'div';
no_domain = prefix+'.0.05.'+str(i)+'.sbml.lm';
no_domain_ori = '../simul4fig6_v2/'+prefix+'.0.05.'+str(i)+'.sbml.lm';
pdf3 = get_mRNA(no_domain,num);
pdf4 = get_mRNA(no_domain_ori, num);

#'''
prefix = 'tand';
no_domain = prefix+'.0.05.'+str(i)+'.sbml.lm';
no_domain_ori = '../simul4fig6_v2/'+prefix+'.0.05.'+str(i)+'.sbml.lm';
pdf5 = get_mRNA(no_domain,num);
pdf6 = get_mRNA(no_domain_ori,num);
#'''

figname = 'Fig_5C.pdf'
plt.rcParams["figure.figsize"] = (14,3)
subplot(1,3,1)
plt.bar(pdf1[:,0],pdf1[:,1], width=0.8*binwidth, alpha=0.4, label='W/O SC')
plt.bar(pdf2[:,0],pdf2[:,1], width=0.8*binwidth, alpha=0.4, label='W/ SC')
plt.ylabel('PDF')
plt.xlabel('Intergenic SC density')
#plt.ylim((0,0.85))
plt.xlim((-0.8,0.8))
plt.legend()
plt.title('Convergent')

subplot(1,3,2)
plt.bar(pdf3[:,0],pdf3[:,1], width=0.8*binwidth, alpha=0.4, label='W/O SC')
plt.bar(pdf4[:,0],pdf4[:,1], width=0.8*binwidth, alpha=0.4, label='W/ SC')
plt.xlabel('Intergenic SC density')
#plt.ylim((0,0.85))
plt.xlim((-0.8,0.8))
#plt.legend()
plt.title('Divergent')

#'''
subplot(1,3,3)
plt.bar(pdf5[:,0],pdf5[:,1], width=0.8*binwidth, alpha=0.4, label='W/O SC')
plt.bar(pdf6[:,0],pdf6[:,1], width=0.8*binwidth, alpha=0.4, label='W/ SC')
#plt.legend()
plt.xlabel('Intergenic SC density')
plt.title('Codirectional')
plt.xlim((-0.8,0.8))
#plt.ylim((0,0.85))

plt.tight_layout()
plt.savefig(figname)
plt.close()


