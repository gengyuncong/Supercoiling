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
def get_pdf(mRNA):
	(bins,edges) = np.histogram(mRNA,np.arange(-0.5,np.max(mRNA)+1.5,1))
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
	counts = np.zeros((len(replicates),2), dtype=int)
	for i,replicate in enumerate(replicates):
		counts[i,0] = np.sum(fp["/Simulations/%s/SpeciesCounts"%replicate][-1, range(10*num, 10*num+90)])/(20.0);#1,2,..,80
		counts[i,1] = np.sum(fp["/Simulations/%s/SpeciesCounts"%replicate][-1, range(10*num+90,11*num)])/(20.0); #81,82,...,160
	pdf = get_pdf(counts[:,0])
	return (pdf, np.mean(counts[:,0]))

def get_mRNA_tand(myfile,num):
	fp = h5py.File(myfile, "r");
	replicates=fp["/Simulations"].keys();
	counts = np.zeros((len(replicates),2), dtype=int)
	for i,replicate in enumerate(replicates):
		counts[i,0] = np.sum(fp["/Simulations/%s/SpeciesCounts"%replicate][-1, range(10*num, 10*num+90)])/(20.0);#1,2,..,80
		counts[i,1] = np.sum(fp["/Simulations/%s/SpeciesCounts"%replicate][-1, range(10*num+90,11*num)])/(20.0); #81,82,...,160
	pdf_1 = get_pdf(counts[:,0])
	pdf_2 = get_pdf(counts[:,1])
	return (pdf_1, pdf_2, np.mean(counts[:,0]), np.mean(counts[:,1]))

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
(pdf1, mean1) = get_mRNA(no_domain,num);
(pdf2, mean2) = get_mRNA(no_domain_ori, num);

prefix = 'div';
no_domain = prefix+'.0.05.'+str(i)+'.sbml.lm';
no_domain_ori = '../simul4fig6_v2/' + prefix+'.0.05.'+str(i)+'.sbml.lm';
(pdf3, mean3) = get_mRNA(no_domain,num);
(pdf4, mean4) = get_mRNA(no_domain_ori, num);

#'''
prefix = 'tand';
no_domain = prefix+'.0.05.'+str(i)+'.sbml.lm';
no_domain_ori = '../simul4fig6_v2/' + prefix+'.0.05.'+str(i)+'.sbml.lm';
(pdf5, pdf6, mean5, mean6) = get_mRNA_tand(no_domain,num);
(pdf5_, pdf6_, mean5_, mean6_) = get_mRNA_tand(no_domain_ori,num);
#'''

figname = 'Fig_5B_mRNA_50_inf.pdf'
plt.rcParams["figure.figsize"] = (14,3)
subplot(1,4,1)
plt.bar(pdf1[:,0],pdf1[:,1], width=0.8, alpha=0.4, label='W/O SC')
plt.bar(pdf2[:,0],pdf2[:,1], width=0.8, alpha=0.4, label='W/ SC')
plt.text(3.5,0.3, 'mean='+str(round(mean1,2)), color=new_colors[0])
plt.text(3.5,0.15, 'mean='+str(round(mean2,2)), color=new_colors[1])
plt.xlabel('mRNA copy number')
plt.ylabel('Probability')
plt.ylim((0,0.85))
plt.xlim((-0.5,10.5))
plt.legend()
plt.title('Convergent')

subplot(1,4,2)
plt.bar(pdf3[:,0],pdf3[:,1], width=0.8, alpha=0.4, label='W/O SC')
plt.bar(pdf4[:,0],pdf4[:,1], width=0.8, alpha=0.4, label='W/ SC')
plt.text(3.5,0.45, 'mean='+str(round(mean3,2)), color=new_colors[0])
plt.text(3.5,0.3, 'mean='+str(round(mean4,2)), color=new_colors[1])
plt.xlabel('mRNA copy number')
plt.ylim((0,0.85))
plt.xlim((-0.5,10.5))
plt.title('Divergent')

#'''
subplot(1,4,3)
plt.bar(pdf5[:,0],pdf5[:,1], width=0.8, alpha=0.4, label='W/O SC')
plt.bar(pdf5_[:,0],pdf5_[:,1], width=0.8, alpha=0.4, label='W/ SC')
plt.text(3.5,0.45, 'mean='+str(round(mean5,2)), color=new_colors[0])
plt.text(3.5,0.3, 'mean='+str(round(mean5_,2)), color=new_colors[1])
plt.xlabel('mRNA copy number')
plt.title('Codirectional (upstream)')
plt.xlim((-0.5,10.5))
plt.ylim((0,0.85))

subplot(1,4,4)
plt.bar(pdf6[:,0],pdf6[:,1], width=0.8, alpha=0.4, label='W/O SC')
plt.bar(pdf6_[:,0],pdf6_[:,1], width=0.8, alpha=0.4, label='W/ SC')
plt.text(3.5,0.45, 'mean='+str(round(mean6,2)), color=new_colors[0])
plt.text(3.5,0.3, 'mean='+str(round(mean6_,2)), color=new_colors[1])
plt.xlabel('mRNA copy number')
plt.title('Codirectional (downstream)')
plt.xlim((-0.5,10.5))
plt.ylim((0,0.85))
#'''

plt.tight_layout()
plt.savefig(figname)
plt.close()


