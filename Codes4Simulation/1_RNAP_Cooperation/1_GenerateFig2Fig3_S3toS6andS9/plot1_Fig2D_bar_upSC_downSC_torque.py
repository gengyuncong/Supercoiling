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

lk0 = 60.0;

def cal_SC(lk):
	SC = (lk-lk0)/lk0;
	return SC

def cal_torque(lk):
	SC = (lk-lk0)/lk0;
	torque = 0
	if abs(SC) < 0.0126:
		torque = 516.33*SC
	elif abs(SC) <= 0.0375:
		torque = 6.51*np.sign(SC)
	else:
		torque = 173.6*SC
	if torque < -10.5:
		return -10.5
	else:
		return torque

def get_hist2(data, binwidth):
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


def get_hist(data, binwidth):
	data = np.array(data)
	H, xedges, yedges = np.histogram2d(data[:,0], data[:,1], (np.arange(-0.5+np.min(data[:,0]), np.max(data[:,0])+1.5, binwidth), np.arange(-0.5,2.5,1)))
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
	all_down=np.empty((0,2), dtype=float)
	all_up=np.empty((0,2), dtype=float)
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
						up0 = cal_SC(Lk[t,p-1])
						down0 = cal_SC(Lk[t,p+1])
						up1 = cal_torque(Lk[t,p-1])
						down1 = cal_torque(Lk[t,p+1])
						torque1=down1-up1
						if (RNAP[t, p]==0):
							all_up=np.vstack((all_up, [up0, 0]))
							all_down=np.vstack((all_down, [down0, 0]))
							all_torque=np.vstack((all_torque, [torque1, 0]))
						else:
							all_up=np.vstack((all_up, [up0, 1]))
							all_down=np.vstack((all_down, [down0, 1]))
							all_torque=np.vstack((all_torque, [torque1, 1]))
		if (i>n):
			(pdf_stallup, pdf_non_stallup) = get_hist2(all_up, 0.05)
			(pdf_stalldown, pdf_non_stalldown) = get_hist2(all_down, 0.05)
			(pdf_stall, pdf_non_stall) = get_hist(all_torque, 5)
			return (pdf_stallup, pdf_non_stallup, pdf_stalldown, pdf_non_stalldown, pdf_stall, pdf_non_stall)

def get_initiaiton_rate(filename):
	fp = h5py.File(filename, "r");
	replicates=fp["/Simulations"].keys();
	initiation=np.array([])
	for i,replicate in enumerate(replicates):
		counts = fp["/Simulations/%s/SpeciesCounts"%replicate]
		S1 = counts[:,13*num+1];
		initiation=np.append(initiation, S1[-1]/(1.0*len(S1)))
	return np.mean(initiation)

lists3 = ['0.001', '0.01', '0.2'];

figname = 'Fig2D_tor.pdf'
plt.rcParams["figure.figsize"] = (12,7)

index_list = range(1,10);
index_list = [1,4,7,2,5,8,3,6,9];
num_list = [98,20,5];
k = 0; j = 0;
for ini_rate in lists3:
	hnormal = 'sens.'+str(ini_rate)+'.50.1.sbml.lm';
	(pdf_stallup, pdf_non_stallup, pdf_stalldown, pdf_non_stalldown, pdf_stall, pdf_non_stall)=get_torque_dist(hnormal,num_list[j])
	j = j+1;

	#plot up SC
	binwidth = 0.05
	ax = plt.subplot(3,3,index_list[k])
	ax.bar(pdf_non_stallup[:,0], pdf_non_stallup[:,1], width=binwidth, color=new_colors[2], label='Normal')
	ax.bar(pdf_stallup[:,0], pdf_stallup[:,1],width=binwidth, bottom=pdf_non_stallup[:,1], color=new_colors[3], label='Stalled')
	ax.set_xlim((-0.9,0.9))
	ax.set_xlabel(r'$\sigma$ (upstream)')
	rate1=get_initiaiton_rate(hnormal)
	rate1=round(rate1,3)
	ax.set_title('Initiation rate='+str(rate1)+'/s')
	if index_list[k] == 1:
		ax.set_ylabel('PDF')
	k = k+1;
	
	#plot down SC
	binwidth = 0.05
	ax = plt.subplot(3,3,index_list[k])
	ax.bar(pdf_non_stalldown[:,0], pdf_non_stalldown[:,1], width=binwidth, color=new_colors[2], label='Normal')
	ax.bar(pdf_stalldown[:,0], pdf_stalldown[:,1], width=binwidth, bottom=pdf_non_stalldown[:,1], color=new_colors[3], label='Stalled')
	ax.set_xlim((-0.9,0.9))
	ax.set_xlabel('$\sigma$ (downstream)')
	if (index_list[k] ==4):
		ax.set_ylabel('PDF')
	k = k+1;
	
	#plot torque
	binwidth = 5
	ax = plt.subplot(3,3,index_list[k])
	ax.bar(pdf_non_stall[:,0], pdf_non_stall[:,1], width=binwidth, color=new_colors[2], label='Normal')
	ax.bar(pdf_stall[:,0], pdf_stall[:,1], width=binwidth, bottom=pdf_non_stall[:,1], color=new_colors[3], label='Stalled')
	ax.axvspan(10.5, 120, facecolor='black', alpha=0.1)
	ax.set_xlim((-120,120))
	ax.set_xlabel(r'Torque (pN$\cdot$nm)')
	if index_list[k]==7:
		ax.set_ylabel('PDF')
	k = k+1;

plt.tight_layout()
plt.savefig(figname)
plt.close()


