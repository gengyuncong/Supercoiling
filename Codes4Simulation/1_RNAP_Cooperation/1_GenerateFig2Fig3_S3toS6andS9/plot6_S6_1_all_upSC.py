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
	density=np.array([])
	for i,replicate in enumerate(replicates):
		S1 = fp["/Simulations/%s/SpeciesCounts"%replicate][:,13*num+1];
		S2 = fp["/Simulations/%s/SpeciesCounts"%replicate][:,13*num+2];
		t1 = np.where(S1[:-1] != S1[1:])[0] + 1;
		t2 = np.where(S2[:-1] != S2[1:])[0] + 1;
		t1 = t1.tolist();
		t2 = t2.tolist();
		if (len(t1)<1 or len(t2)<1):
			next;
		else:
			length = min(len(t1), len(t2))
			for ii in range(0, length):
				elongation = np.append(elongation, int(3000/(t2[ii] - t1[ii])));

		counts=np.append(counts, S1[-1]/(1.0*len(S1)))

		count = fp["/Simulations/%s/SpeciesCounts"%replicate]
		RNAP = count[:, range(num, num+80)]+count[:, range(2*num,2*num+80)]+count[:, range(8*num, 8*num+80)]
		for t in range(0,len(RNAP)):
			if (np.sum(RNAP[t,:])>0):
				RNAP_index = np.where(RNAP[t,:]>0)[0]
				if (len(RNAP_index)>1):
					density=np.append(density, np.diff(RNAP_index))

	return (np.mean(counts), np.mean(elongation), np.std(elongation), np.mean(density), np.std(density))

lk0 = 60.0;
def cal_torque(lk):
	SC = (lk-lk0)/lk0;
	return SC

binwidth = 0.1

def get_hist(data):
	data = np.array(data)
	H, xedges, yedges = np.histogram2d(data[:,0], data[:,1], (np.arange(-0.5*binwidth+np.min(data[:,0]), np.max(data[:,0])+1.5*binwidth, binwidth), np.arange(-0.5,2.5,1)))
	centers = (xedges[:-1]+xedges[1:])/2
	xdata=centers
	pdf1 = np.zeros((len(H),2),dtype=float)
	pdf1[:,0]=centers
	pdf1[:,1] = H[:,0].astype(double)/sum(sum(H))*(centers[1]-centers[0])

	pdf2 = np.zeros((len(H),2),dtype=float)
	pdf2[:,0]=centers
	pdf2[:,1] = H[:,1].astype(double)/sum(sum(H))*(centers[1]-centers[0])
	return pdf1, pdf2

def get_torque_dist(filename,n):
	fp = h5py.File(filename, "r");
	replicates=fp["/Simulations"].keys();
	all_torque=np.empty((0,2), dtype=float)
	for i,replicate in enumerate(replicates):
		print(i)
		counts = fp["/Simulations/%s/SpeciesCounts"%replicate]
		RNAP = counts[:, range(num, num+80)]
		RNAP_stall = counts[:, range(2*num,2*num+80)]
		RNAP_tmp = counts[:, range(8*num, 8*num+80)]
		Lk = counts[:, range(3*num, 3*num+80)]
		all_RNAP = np.array(RNAP) + np.array(RNAP_stall) + np.array(RNAP_tmp);
		all_RNAP_t = np.sum(all_RNAP, axis=1)

		index_t = np.where(all_RNAP_t>0)[0]
		if (len(index_t)==0):
			next;
		else:
			for t in index_t:
				all_RNAP_p = all_RNAP[t,:]
				index_p = np.where(all_RNAP_p>0)[0]
				for p in index_p:
					if (p >= 79):
						next;
					elif (Lk[t,p-1] >0 and Lk[t,p+1] >0):
						up1 = cal_torque(Lk[t,p-1])
#						down1 = cal_torque(Lk[t,p+1])
						if (RNAP[t, p]==0):
							all_torque=np.vstack((all_torque, [up1, 0]))
						else:
							all_torque=np.vstack((all_torque, [up1, 1]))
		if (i>n):
			(pdf_stall, pdf_non_stall) = get_hist(all_torque)
			return (pdf_stall, pdf_non_stall)

def get_initiaiton_rate(filename):
	fp = h5py.File(filename, "r");
	replicates=fp["/Simulations"].keys();
	initiation=np.array([])
	for i,replicate in enumerate(replicates):
		counts = fp["/Simulations/%s/SpeciesCounts"%replicate]
		S1 = counts[:,13*num+1];
		initiation=np.append(initiation, S1[-1]/(1.0*len(S1)))
	return np.mean(initiation)

lists3 = ['0.001', '0.005', '0.01', '0.02', '0.05', '0.08', '0.1', '0.15', '0.2'];

figname = 'S6_1_upSC.pdf'
#plt.rcParams["figure.figsize"] = (25,10)
plt.rcParams["figure.figsize"] = (13,5)

index_list = range(1,10);
num_list = [98,50,20,5,5,5,5,5,5,5];
k = 0;
for ini_rate in lists3:
	hnormal = 'sens.'+str(ini_rate)+'.50.1.sbml.lm';
	(pdf_stall, pdf_non_stall)=get_torque_dist(hnormal,num_list[k])
	rate1=get_initiaiton_rate(hnormal)
	rate1=round(rate1,3)
	ax = plt.subplot(2,5,index_list[k])
	ax.set_title('Initiation rate='+str(rate1)+'/s')
	ax.bar(pdf_non_stall[:,0], pdf_non_stall[:,1], width=binwidth, color=new_colors[2], label='Normal')
	ax.bar(pdf_stall[:,0], pdf_stall[:,1], width=binwidth, bottom=pdf_non_stall[:,1], color=new_colors[3], label='Stalled')
	ax.set_xlim((-0.9,0.9))
	ax.set_xlabel(r'$\sigma$'+' (upstream)')
	if k == 0 or k == 5:
		ax.set_ylabel('PDF')
	k = k+1;

plt.tight_layout()
plt.savefig(figname)
plt.close()


