from __future__ import division
import matplotlib
matplotlib.use('Agg')
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import h5py
import math
import sys
from scipy.optimize import curve_fit
from scipy.misc import factorial

fontSize=16
matplotlib.rcParams.update({"axes.formatter.limits": (-4,4), "svg.fonttype" : "none", 'pdf.fonttype':42,'font.family':'sans-serif','font.sans-serif':'Helvetica','font.size': fontSize, "axes.titlesize": fontSize, "xtick.labelsize": fontSize, "ytick.labelsize": fontSize,'text.usetex':True,'text.latex.preamble':[r'\usepackage{sansmath}',r'\sansmath']})
plotStyles={"markersize":8,"markeredgewidth":1.0,"linewidth":3.0}
stepStyles={"markersize":20,"markeredgewidth":0.0,"linewidth":1.0,"where":"post"}
barStyles={"width":0.65, "linewidth":0, "align":"center", "color":"steelblue"}

new_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']

##############################################################
def my_func(x,upper):
	if x >= 0:
		return 0.001
	elif x <= -0.06:
		return upper
	else:
		return ((-x)/0.06*(upper-0.001)+0.001)

def poisson(k, lamb):
	return (lamb**k/factorial(k)) * np.exp(-lamb)

def get_pdf(mRNA):
	(bins,edges) = np.histogram(mRNA,np.arange(-0.5,np.max(mRNA)+1.5,1))
	centers = (edges[:-1]+edges[1:])/2
	xdata=centers
	ydata=bins/sum(bins)
	popt, pcov = curve_fit(poisson, xdata, ydata)
	pdf = np.zeros((len(bins),2),dtype=float)
	pdf[:,0]=centers
	pdf[:,1] = bins.astype(double)/(np.sum(bins)*(centers[1]-centers[0]))
	return (pdf, popt[0])

def get_mRNA(filename):
	fp = h5py.File(filename, "r");
	replicates=fp["/Simulations"].keys();
	counts=np.array([], dtype=float)
	for i,replicate in enumerate(replicates):
		c = sum(fp["/Simulations/%s/SpeciesCounts"%replicate][-1,range(10*num, 11*num)])
		counts=np.append(counts, c/20.0)
	mymean = np.mean(counts)
	fano=np.var(counts)/np.mean(counts)
	(pdf, lamda) = get_pdf(counts)
	return (pdf, lamda, mymean, fano)

def get_traj(filename, replicate):
	fp = h5py.File(filename, "r");
	barrier1 = 2; #3,4,...,68 -> 2,3,...,67 -> range(2,68)
	barrier2 = 68; #
	counts=fp["/Simulations/%07d/SpeciesCounts"%replicate];
	times=fp["/Simulations/%07d/SpeciesCountTimes"%replicate];
	ini_count=[0]*len(times);
	upstream_SC=[0]*len(times);
	downstream_SC=[0]*len(times);
	SC_count=[0]*len(times);
	ini_event = 0
	t = 0;
	while (t<=times[-1]):
		effective_len=sum(counts[t, range(barrier1, barrier2)]);
		SC_count[t] = (sum(counts[t, range(3*num+barrier1, 3*num+barrier2)]) - lk0*effective_len)/(lk0*effective_len);
		if (counts[t, 13*num+1] > ini_event):
			ini_count[t] = 1;
			ini_event = counts[t, 13*num+1]; 
		L1 = counts[t, range(barrier1, barrier2)]; #6,7,8,...,54
		L1 = L1.tolist();
		if (len(L1) == len([i for i in L1 if i > 0])):
			upstream_SC[t] = SC_count[t];
			downstream_SC[t] = upstream_SC[t];
		else:
			L2 = L1[::-1];
			index1 = L1.index(0);
			index2 = L2.index(0);
			upstream_SC[t] = (sum(counts[t, range(3*num+barrier1, 3*num+barrier1+index1)]) - lk0*index1)/(lk0*index1);
			downstream_SC[t] = (sum(counts[t, range(3*num+barrier2-index2, 3*num+barrier2)]) - lk0*index2)/(lk0*index2);
		t=t+1;
	ini_rate=round(sum(ini_count)/len(times),3)
	return (times, counts[:,1], ini_count, upstream_SC, downstream_SC, SC_count)

##############################################################
lk0 = 60;
num = 70;

x=np.linspace(-0.10,0.05,1000)
my_func=np.vectorize(my_func)

plt.rcParams["figure.figsize"] = (12,4)

filename = 'normal.0.1.sbml.lm'
(times, loop, ini_count, upstream_SC, downstream_SC, SC_count) = get_traj(filename, 1);
ax = plt.subplot(2,2,1)
ax.set_title('Strong promoter (open)')
ax.step(times, ini_count, label='initiation', color=new_colors[0], **stepStyles)
ax.step(times, loop, label="loop open", color=new_colors[3], **stepStyles)
ax.set_xlabel('Time (s)'); 
ax.set_ylabel('Number');
#ax.legend(loc="upper right", ncol=2, framealpha=0.7, fontsize=12)
ax = plt.subplot(2,2,3)
ax.step(times, SC_count, label='average', color=new_colors[0], **stepStyles)
ax.step(times, upstream_SC, label='upstream', color=new_colors[2], **stepStyles)
ax.step(times, downstream_SC, label='downstream', color=new_colors[3], **stepStyles)
ax.set_ylabel(r'$\sigma$')
ax.set_xlabel('Time (s)')
ax.set_ylabel(r'$\sigma$')
ax.set_ylim((-0.65,0.78)) 
ax.legend(loc="upper right", ncol=3, framealpha=0.7, fontsize=12)

filename = 'domain.0.1.sbml.lm'
(times, loop, ini_count, upstream_SC, downstream_SC, SC_count) = get_traj(filename, 1);
ax = plt.subplot(2,2,2)
ax.set_title('Strong promoter (dynamic loop)')
ax.step(times, ini_count, label='initiation', color=new_colors[0], **stepStyles)
ax.step(times, loop, label="loop open", color=new_colors[3], **stepStyles)
ax.set_xlabel('Time (s)');
ax.set_ylabel('Number');
ax.legend(loc="upper right", ncol=1, framealpha=0.5, fontsize=10)
ax = plt.subplot(2,2,4)
ax.step(times, SC_count, label='average', color=new_colors[0], **stepStyles)
ax.step(times, upstream_SC, label='upstream', color=new_colors[2], **stepStyles)
ax.step(times, downstream_SC, label='downstream', color=new_colors[3], **stepStyles)
ax.set_xlabel('Time (s)'); 
ax.set_ylabel(r'$\sigma$')
ax.set_ylim((-0.65,0.78))
#ax.legend(loc="upper right", ncol=1, framealpha=0.5, fontsize=10)

figname='Fig_4GH.pdf'
plt.tight_layout()
plt.savefig(figname)
plt.close()

