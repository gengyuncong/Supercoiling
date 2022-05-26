from __future__ import division
import matplotlib
matplotlib.use('Agg')
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import h5py
import math
import sys

fontSize=16
matplotlib.rcParams.update({"axes.formatter.limits": (-4,4), "svg.fonttype" : "none", 'pdf.fonttype':42,'font.family':'sans-serif','font.sans-serif':'Helvetica','font.size': fontSize, "axes.titlesize": fontSize, "xtick.labelsize": fontSize, "ytick.labelsize": fontSize,'text.usetex':True,'text.latex.preamble':[r'\usepackage{sansmath}',r'\sansmath']})
plotStyles={"markersize":8,"markeredgewidth":1.0,"linewidth":3.0}
stepStyles={"markersize":20,"markeredgewidth":3.0,"linewidth":3.0,"where":"post"}
barStyles={"width":0.65, "linewidth":0, "align":"center", "color":"steelblue"}

new_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']

##############################################################
'''
0-1:	DNA
1-2:	RNAP
2-3:	RNAP_tmp
3-4:	Twist
4-5:	Gyrase_unbind
5-6:	Gyrase_bind
6-7:	TopoI_unbind
7-8:	TopoI_bind
8-9:	RNAP_stall
9-10:	mRNA	

'''
##############################################################
def cal_domain_ratio(tmp, ini_count, ratio1, ratio2, time):
	S1_total=np.sum(ini_count)
	index = np.where(tmp>0)[0] #when domain is open/when gyrase binds
	S1_open = np.sum(ini_count[index]) #the initiation events when domain is open/when gyrase binds
	S1_closed = S1_total-S1_open; #the initiation events when domain is closed/when gyrase unbinds
	t_open = len(index);   #the time of domain open/gyrase binds
#	print(t_open)
	t_closed = time-t_open; #the time of domain closed/gyrase unbinds
	if (t_open!=0):
		ratio1=np.append(ratio1,S1_open/(1.0*t_open))  #initiation rate when domain is open/gyrase binds
	if (t_closed!=0):
		ratio2=np.append(ratio2,S1_closed/(1.0*t_closed))  #initiation rate when domain is closed/gyrase unbinds
	return (ratio1, ratio2)


num=70;
#filename = sys.argv[1]
def get_mRNA4unlooped(filename):
	fp = h5py.File(filename, "r");
	replicates=fp["/Simulations"].keys();
	mrna=np.array([], dtype=float)
	for i,replicate in enumerate(replicates):
		counts = fp["/Simulations/%s/SpeciesCounts"%replicate]
		c = counts[-1, 10*num+54];
		mrna=np.append(mrna, c);
	mymean = np.mean(mrna)
	fano=np.var(mrna)/np.mean(mrna)
	return (mymean, fano)

def new_array(gyrase, binsize, time):
	gyrase = np.array(gyrase)
	index = np.where(gyrase>0)[0]
	(bins,edges) = np.histogram(index, np.arange(0,time+1,binsize))
	centers = (edges[:-1]+edges[1:])/2
	pdf = np.zeros((len(bins),2),dtype=float)
	pdf[:,0]=centers
	pdf[:,1] = bins.astype(double)/(np.sum(bins)*(centers[1]-centers[0]))
	index2 = np.where(pdf[:,1]>0)[0]
	if (len(index2)==0):
		print("no binding")
	else:
		new_array = np.zeros(time)
		for i in index2:
			start = int(edges[i])
			end = int(edges[i+1])
			new_array[range(start, end)] = 1;
		return new_array

def get_mRNA(filename):
	fp = h5py.File(filename, "r");
	replicates=fp["/Simulations"].keys();
	mrna=np.array([], dtype=float)
	ratio_open=np.array([]);
	ratio_closed=np.array([]);
	ratio_bound=np.array([]);
	ratio_unbound=np.array([]);
	ratio_bound1=np.array([]);
	ratio_unbound1=np.array([]);
	for i,replicate in enumerate(replicates):
		counts = fp["/Simulations/%s/SpeciesCounts"%replicate]
		counts = counts[750:,:];
		time=len(counts[:,0]);
		# calculate initiation
		ini_count=np.zeros(time);
		ini_event = counts[0, 13*num+1];
		for t in range(0,time):
			if (counts[t, 13*num+1] > ini_event):
				ini_count[t] = counts[t, 13*num+1]-ini_event;
				ini_event = counts[t, 13*num+1];
		# domain open and closed
		domain = np.array(counts[:,1]);
		(ratio_open, ratio_closed) = cal_domain_ratio(domain, ini_count, ratio_open, ratio_closed, time)
		# gyrase bound
#		gyrase = np.sum(counts[:,range(5*num+55,5*num+68)], axis=1) #55,56,...,67 -> 56,57,...,68 
		gyrase = np.sum(counts[:,range(5*num+2, 5*num+15)], axis=1) #2,3,..,14 -> 2,3,...,15
#		gyrase = np.sum(counts[:,range(5*num+15,5*num+55)], axis=1) #15,16,...,54 -> 16,17,...,55
		(ratio_bound, ratio_unbound) = cal_domain_ratio(gyrase, ini_count, ratio_bound, ratio_unbound, time)
		# topoI bound
		topoI = np.sum(counts[:,range(7*num+2, 7*num+15)], axis=1)
		(ratio_bound1, ratio_unbound1) = cal_domain_ratio(topoI, ini_count, ratio_bound1, ratio_unbound1, time)
		# mRNA
		c = counts[-1, 10*num+54];
		mrna=np.append(mrna, c);
	mymean = np.mean(mrna)
	fano=np.var(mrna)/np.mean(mrna)
	return (mymean, fano, np.mean(ratio_open), np.std(ratio_open), np.mean(ratio_closed), np.std(ratio_closed), np.mean(ratio_bound), np.std(ratio_bound), np.mean(ratio_unbound), np.std(ratio_unbound), np.mean(ratio_bound1), np.std(ratio_bound1), np.mean(ratio_unbound1), np.std(ratio_unbound1))

def append_array(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14, a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14):
	a1 = np.append(a1,t1)
	a2 = np.append(a2,t2)
	a3 = np.append(a3,t3)
	a4 = np.append(a4,t4)
	a5 = np.append(a5,t5)
	a6 = np.append(a6,t6)
	a7 = np.append(a7,t7)
	a8 = np.append(a8,t8)
	a9 = np.append(a9,t9)
	a10 = np.append(a10,t10)
	a11 = np.append(a11,t11)
	a12 = np.append(a12,t12)
	a13 = np.append(a13,t13)
	a14 = np.append(a14,t14)
	return (a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14);

def swap_order(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14, index0):
	a1 = a1[index0]
	a2 = a2[index0]
	a3 = a3[index0]
	a4 = a4[index0]
	a5 = a5[index0]
	a6 = a6[index0]
	a7 = a7[index0]
	a8 = a8[index0]
	a9 = a9[index0]
	a10 = a10[index0]
	a11 = a11[index0]
	a12 = a12[index0]
	a13 = a13[index0]
	a14 = a14[index0]
	return (a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14);

lists1 = ['0.001', '0.005', '0.01', '0.02', '0.05', '0.08', '0.1', '0.15', '0.2', '0.3', '0.4', '0.5'];
#lists1 = ['0.001', '0.005', '0.01', '0.02', '0.05', '0.08', '0.1', '0.15', '0.2'];
lists2 = ['0.001', '0.005', '0.01', '0.02', '0.05', '0.08', '0.1', '0.15', '0.2'];

mean = []
fano = []

mean_d = []
fano_d = []
ratio_open_d = []
ratio_open_d_error = []
ratio_closed_d = []
ratio_closed_d_error = []
ratio_bound_d = []
ratio_bound_d_error = []
ratio_unbound_d = []
ratio_unbound_d_error = []
ratio_bound_t = []
ratio_bound_t_error = []
ratio_unbound_t = []
ratio_unbound_t_error = []

for i in lists1:
	filenamed = 'domain.'+str(i)+'.sbml.lm'
	(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14) = get_mRNA(filenamed);
	(mean_d, fano_d, ratio_open_d, ratio_open_d_error, ratio_closed_d, ratio_closed_d_error, ratio_bound_d, ratio_bound_d_error, ratio_unbound_d, ratio_unbound_d_error, ratio_bound_t, ratio_bound_t_error, ratio_unbound_t, ratio_unbound_t_error) = append_array(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,mean_d, fano_d, ratio_open_d, ratio_open_d_error, ratio_closed_d, ratio_closed_d_error, ratio_bound_d, ratio_bound_d_error, ratio_unbound_d, ratio_unbound_d_error, ratio_bound_t, ratio_bound_t_error, ratio_unbound_t, ratio_unbound_t_error);

for i in lists2:
	filename = 'normal.'+str(i)+'.sbml.lm'
	(mymean, myfano) = get_mRNA4unlooped(filename)
	mean = np.append(mean, mymean)
	fano = np.append(fano, myfano)
	
index = np.argsort(mean);
mean = mean[index]
fano = fano[index]

index = np.argsort(mean_d);
(mean_d, fano_d, ratio_open_d, ratio_open_d_error, ratio_closed_d, ratio_closed_d_error, ratio_bound_d, ratio_bound_d_error, ratio_unbound_d, ratio_unbound_d_error, ratio_bound_t, ratio_bound_t_error, ratio_unbound_t, ratio_unbound_t_error)=swap_order(mean_d, fano_d, ratio_open_d, ratio_open_d_error, ratio_closed_d, ratio_closed_d_error, ratio_bound_d, ratio_bound_d_error, ratio_unbound_d, ratio_unbound_d_error, ratio_bound_t, ratio_bound_t_error, ratio_unbound_t, ratio_unbound_t_error, index)

#index = np.argsort(mean_d);
#(mean_d, fano_d, ratio_open_d, ratio_open_d_error, ratio_closed_d, ratio_closed_d_error, ratio_bound_d, ratio_bound_d_error, ratio_unbound_d, ratio_unbound_d_error)=swap_order(mean_d, fano_d, ratio_open_d, ratio_open_d_error, ratio_closed_d, ratio_closed_d_error, ratio_bound_d, ratio_bound_d_error, ratio_unbound_d, ratio_unbound_d_error, index)


figname = 'Fig_4IJ.pdf'
plt.rcParams["figure.figsize"] = (10,3.5)

ax1 = plt.subplot(1,2,1)
ax1.axhline(y=1, color=new_colors[7], linewidth=2)
ax1.plot(mean, fano, 'o-', color=new_colors[0], linewidth=2, label='Unlooped')
ax1.plot(mean_d, fano_d, 'o-', color=new_colors[1], linewidth=2, label='Looped')
ax1.set_ylabel('Fano factor')	
ax1.set_xlabel('Mean mRNA copy number')
ax1.set_xscale('log')
#ax1.set_ylim((0,5.2))
ax1.legend(loc='upper left')

ax3 = plt.subplot(1,2,2)
ax3.errorbar(mean_d, ratio_open_d, yerr=ratio_open_d_error, fmt='-o', color=new_colors[2], linewidth=2, label='Loop open')
ax3.errorbar(mean_d, ratio_closed_d, yerr=ratio_closed_d_error, fmt='-o', color=new_colors[3], linewidth=2, label='Loop closed')
ax3.set_ylabel('Initiation rate (/s)')
ax3.set_xlabel('Mean mRNA copy number')
ax3.set_xscale('log')
ax3.legend(loc='upper left')

'''
ax2 = plt.subplot(2,2,3)
ax2.errorbar(mean_d, ratio_bound_d, yerr=ratio_bound_d_error, fmt='-o', color=new_colors[2], linewidth=2, label='Gyrase bound')
ax2.errorbar(mean_d, ratio_unbound_d, yerr=ratio_unbound_d_error, fmt='-o', color=new_colors[3], linewidth=2, label='Gyrase unbound')
ax2.set_ylabel('Initiation rate (/s)')
ax2.set_xlabel('Mean mRNA copy number')
ax2.set_xscale('log')
ax2.legend(loc='upper left')

ax4 = plt.subplot(2,2,4)
ax4.errorbar(mean_d, ratio_bound_t, yerr=ratio_bound_t_error, fmt='-o', color=new_colors[2], linewidth=2, label='TopoI bound')
ax4.errorbar(mean_d, ratio_unbound_t, yerr=ratio_unbound_t_error, fmt='-o', color=new_colors[3], linewidth=2, label='TopoI unbound')
ax4.set_ylabel('Initiation rate (/s)')
ax4.set_xlabel('Mean mRNA copy number')
ax4.set_xscale('log')
ax4.legend(loc='upper left')
'''
plt.tight_layout()
plt.savefig(figname)
plt.close()


