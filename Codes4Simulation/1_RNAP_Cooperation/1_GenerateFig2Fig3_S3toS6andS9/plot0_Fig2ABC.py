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

num=149;

def get_traj(S1, RNAP, t_begin, t_end):
	S1_tmp = S1[t_begin]
	index1 = 0;
	pos = np.array([])
	for t in range(t_begin, t_end+1):
		RNAP_t = RNAP[t,:];
		S1_t = S1[t];
		if (S1_t > S1_tmp):
			index1 += S1_t-S1_tmp;
			S1_tmp = S1_t;
		if (len(np.where(RNAP_t==1)[0])-1 < index1):
			break;
		pos_tmp = np.where(RNAP_t==1)[0][index1]
		pos = np.append(pos, pos_tmp)
	return pos

def interpolate(traj, index):
	if (int(index) == index):
		return traj[int(index)]
	else:
		up = int(math.ceil(index))
		down = int(math.floor(index))
		x = traj[down]+(traj[up]-traj[down])*(index-down)
		return x

def stretch(traj, my_len):
	tmp = np.zeros(my_len);   
	for i in range(0, my_len):  #0,1,2  ->  #0,1,2,3 ,new_len-1=3, ori_len-1=2, 1/3*2=2/3, 2/3*2=4/3, ...
		index_tmp = i/(1.0*(my_len-1))*(len(traj)-1)
		tmp[i] = interpolate(traj, index_tmp);
	return tmp	


def cal_corr(traj_pre, traj_now):
	len1 = len(traj_pre)
	len2 = len(traj_now)
	if (len1 > len2):
		traj_pre = traj_pre[range(0, len2)];
	elif (len1 < len2):
		traj_now = traj_now[range(0, len1)];
	if (len(traj_now)!=len(traj_pre)):
		print("not match")
	else:
		corr = np.corrcoef([traj_pre,traj_now])[0][1]
		if (np.isnan(corr)):
			print("Couldn't calculate correlation")
		else:
			return corr

def cal_correlation(S1, S2, RNAP, correlation, tt):
	traj_pre = np.array([])
	S1_pre = S1[tt];
	while (S1_pre<S1[-1]):
		if (len(np.where(S1>S1_pre)[0]) == 0):
			break;
		t_begin = int(np.where(S1>S1_pre)[0][0])
		RNAP_begin = sum(RNAP[t_begin,:])
		S2_end = S2[t_begin]+RNAP_begin;
		if (len(np.where(S2==S2_end)[0])==0):
			break;
		t_end = np.where(S2==S2_end)[0][0] - 1;
		traj_now = get_traj(S1, RNAP, t_begin, t_end);
		if (len(traj_pre) == 0):
			traj_pre = traj_now;
		else:
			corr = cal_corr(traj_pre, traj_now);
			correlation = np.append(correlation, corr)
			traj_pre = traj_now;
		S1_pre = S1_pre+1;
	return correlation

def get_all_statistics(myfile):
	fp = h5py.File(myfile, "r");
	replicates=fp["/Simulations"].keys();
	initiation=np.array([])
	elongation=np.array([])
	density=np.array([])
	correlation=np.array([])
	tt = 750;
	for i,replicate in enumerate(replicates):
		count = fp["/Simulations/%s/SpeciesCounts"%replicate]
		S1 = count[:,13*num+1];
		S2 = count[:,13*num+2];
		RNAP = count[:, range(num, 2*num)]+count[:, range(2*num,3*num)]+count[:, range(8*num, 9*num)]
#		RNAP = count[:, range(num, num+80)]+count[:, range(2*num,2*num+80)]+count[:, range(8*num, 8*num+80)]
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
					elongation = np.append(elongation, int(3000/(1.0*(t2-t1))))
		#cal RNAP space
		for t in range(0,len(RNAP)):
			if (np.sum(RNAP[t,:])>0):
				RNAP_index = np.where(RNAP[t,:]>0)[0]
				if (len(RNAP_index)>1):
					density=np.append(density, np.diff(RNAP_index)-1)
		#cal correlation
		correlation=cal_correlation(S1, S2, RNAP, correlation, tt);
	return (np.mean(initiation), np.mean(elongation), np.std(elongation), density*60, np.mean(correlation), np.std(correlation))


new_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']

lists3 = ['0.001', '0.005', '0.01', '0.02', '0.05', '0.08', '0.1', '0.15', '0.2'];
n = 9;
initiation_list = np.zeros(n); #
elong_list = np.zeros(n); #elongation
elong_error = np.zeros(n); #elongation
correlation_list = np.zeros(n); #
correlation_error = np.zeros(n); 
d_mean = np.zeros(n);
d_error = np.zeros(n);
d_median = np.zeros(n);
d_first = np.zeros(n);
d_third = np.zeros(n);
d_IQR = np.zeros(n);

k = 0;
for ini_rate in lists3:
	hnormal = 'sens.'+str(ini_rate)+'.50.1.sbml.lm';
	(initiation_list[k], elong_list[k], elong_error[k], density_tmp, correlation_list[k], correlation_error[k]) = get_all_statistics(hnormal);
	d_mean[k] = np.mean(density_tmp)
	d_error[k] = np.std(density_tmp)
	d_median[k] = np.median(density_tmp)
	d_first[k] = np.quantile(density_tmp, 0.25, interpolation='nearest')
	d_third[k] = np.quantile(density_tmp, 0.75, interpolation='nearest')
	d_IQR[k] = d_third[k]-d_first[k]
	k = k+1;

figname = 'Fig2ABC_v2.pdf'
plt.rcParams["figure.figsize"] = (12,4)

ax1 = plt.subplot(1,3,1)
ax1.errorbar(initiation_list, elong_list, yerr=elong_error, fmt='-o', color=new_colors[0], linewidth=2)
ax1.axvspan(initiation_list[-5], initiation_list[-1], facecolor=new_colors[0], alpha=0.1)
ax1.set_ylabel('Apparent elongation rate (bp/s)')
ax1.set_xlabel('Initiation rate (/s)')
ax1.set_ylim((0,40))

ax3 = plt.subplot(1,3,2)
corr_min = min(correlation_list)
corr_max = max(correlation_list)
for i in range(0,len(correlation_list)):
	correlation_list[i] = (correlation_list[i]-corr_min)/(corr_max-corr_min);
ax3.plot(initiation_list, correlation_list, '-o', color=new_colors[0], linewidth=2)
#ax3.errorbar(initiation_list, correlation_list, yerr=correlation_error, fmt='-o', color=new_colors[0], linewidth=2)
ax3.axvspan(initiation_list[-5], initiation_list[-1], facecolor=new_colors[0], alpha=0.1)
ax3.set_ylabel('Normalized correlation of motions')
ax3.set_xlabel('Initiation rate (/s)')

ax2 = plt.subplot(1,3,3)
ax2.errorbar(initiation_list, d_mean, yerr=d_error, fmt='-o', color=new_colors[0], linewidth=2)
ax2.bar(initiation_list, d_IQR, bottom=d_first, width=0.004, fill=False, edgecolor=new_colors[0], linewidth=1.5)
for i in range(0,len(d_median)):
	y_tmp = d_median[i];
	x_min = initiation_list[i] - 0.004*0.5;
	x_max = initiation_list[i] + 0.004*0.5;
	ax2.hlines(y_tmp, x_min, x_max, colors=new_colors[0], linewidth=1.5)

ax2.axvspan(initiation_list[-5], initiation_list[-1], facecolor=new_colors[0], alpha=0.1)
ax2.set_ylabel('Inter RNAP distance (bp)')	
ax2.set_xlabel('Initiation rate (/s)')
ax2.set_ylim((0,2000))

plt.tight_layout()
plt.savefig(figname)
plt.close()


