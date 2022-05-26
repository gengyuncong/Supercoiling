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
import pandas as pd

fontSize=16
matplotlib.rcParams.update({"axes.formatter.limits": (-4,4), "svg.fonttype" : "none", 'font.family':'MathJax_SansSerif', 'font.size': fontSize, "axes.titlesize": fontSize, "xtick.labelsize": fontSize, "ytick.labelsize": fontSize,'text.usetex':False,'text.latex.preamble':[r'\usepackage{sansmath}',r'\sansmath']})
#fontSize=20
#matplotlib.rcParams.update({"axes.formatter.limits": (-4,4), "svg.fonttype" : "none", 'pdf.fonttype':42,'font.family':'DejaVu Sans','font.sans-serif':'Helvetica','font.size': fontSize, "axes.titlesize": fontSize, "xtick.labelsize": fontSize, "ytick.labelsize": fontSize,'text.usetex':False,'text.latex.preamble':[r'\usepackage{sansmath}',r'\sansmath']})
plotStyles={"markersize":8,"markeredgewidth":1.0,"linewidth":3.0}
stepStyles={"markersize":20,"markeredgewidth":3.0,"linewidth":3.0,"where":"post"}
barStyles={"width":0.65, "linewidth":0, "align":"center", "color":"steelblue"}

##############################################################
import os
from os import listdir
from os.path import isfile, join

num=149;

def get_elongation_rate(myfile):
	fp = h5py.File(myfile, "r");
	replicates=fp["/Simulations"].keys();
	initiation=np.array([])
	elongation=np.array([])
	tt = 750;
	for i,replicate in enumerate(replicates):
		S1 = fp["/Simulations/%s/SpeciesCounts"%replicate][:,13*num+1];
		S2 = fp["/Simulations/%s/SpeciesCounts"%replicate][:,13*num+2];
		S11 = S1[tt:] - S1[tt];
		initiation=np.append(initiation, S11[-1]/(1.0*len(S11)))
		if (S1[-1]<0 or S2[-1]<0):
			next;
		else:
			for ii in range(1,S2[-1]+1):
				t1 = np.where(S1==ii)[0][0]
				t2 = np.where(S2==ii)[0][0]
				if (t1>tt):
					elongation = np.append(elongation, int(3000/(1.0*(t2-t1))))
	return (initiation, elongation)

new_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']

ini_lists = ['0.001', '0.005', '0.01', '0.02', '0.05', '0.08', '0.1', '0.15', '0.2'];
D_lists = ['5', '15.8', '50', '158', '500', '1580', '5000'];
topo_lists = ['0.001', '0.01', '0.1',  '1', '10', '100']; 

prefix = {"5":"sens.", "15.8":"newD/sens.", "50":"sens.", "158":"newD/sens.", "500":"sens.", "1580":"newD/sens.", "5000":"sens."};
D_dict = {"5":"0.002","15.8":"0.006", "50": "0.02", "158":"0.06", "500":"0.2", "1580":"0.6","5000":"2"};

df = pd.DataFrame(columns = ['TopoI', 'Diffusion','initiaiton','elongation']);
for topo in topo_lists:
	figname = 'elongation_meta_topo'+topo+'.pdf'
	plt.rcParams["figure.figsize"] = (8,6)
	ax1 = plt.subplot(1,1,1)
	j = 0;
	n = len(ini_lists);
	for D in D_lists:
		initiation_list = np.zeros(n); #
		elongation_mean = np.zeros(n);
		elongation_std = np.zeros(n); 
		k = 0;
		for rot_rate in ini_lists:
			filename = prefix[D]+rot_rate+'.'+D +'.'+topo + '.sbml.lm';
			#sens.0.05.5000.0.01.sbml.lm
			(initiation, elongation) = get_elongation_rate(filename); 
			initiation_list[k] = np.mean(initiation);
			elongation_mean[k] = np.mean(elongation);
			elongation_std[k] = np.std(elongation);
			k = k+1;
		my_arr = np.stack((np.repeat(topo, n),np.repeat(D, n),initiation_list, elongation_mean), axis=-1);
		df_tmp = pd.DataFrame(my_arr, columns = ['TopoI', 'Diffusion','initiaiton','elongation']);
		df = pd.concat([df, df_tmp], ignore_index=True)

		plt.errorbar(initiation_list, elongation_mean, yerr = elongation_std, color=new_colors[j], linewidth=2, label='D='+D_dict[D]+' '+r'$um^2/s$', alpha = 1)
	#ax1.axvspan(initiation_list[-5], initiation_list[-1], facecolor=new_colors[0], alpha=0.1)
		j = j+1; 
	
	plt.ylabel('Elongation rate (bp/s)')
	plt.xlabel('Initiation rate (/s)')
	plt.legend(ncol=2, fontsize=12)
	plt.ylim((0,60))
	plt.title('TopoI unbinding rate='+topo+' /s')
	plt.tight_layout()
	plt.savefig(figname)
	plt.close()

df.to_csv('elongation_meta.csv')
