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

fontSize=16
matplotlib.rcParams.update({"axes.formatter.limits": (-4,4), "svg.fonttype" : "none", 'pdf.fonttype':42,'font.family':'sans-serif','font.sans-serif':'Helvetica','font.size': fontSize, "axes.titlesize": fontSize, "xtick.labelsize": fontSize, "ytick.labelsize": fontSize,'text.usetex':True,'text.latex.preamble':[r'\usepackage{sansmath}',r'\sansmath']})
plotStyles={"markersize":8,"markeredgewidth":1.0,"linewidth":3.0}
stepStyles={"markersize":20,"markeredgewidth":3.0,"linewidth":3.0,"where":"post"}
barStyles={"width":0.65, "linewidth":0, "align":"center"}

new_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']

##############################################################
'''
0-1:	DNA			*
1-2:	RNAP		*
2-3: 	RNAP_stall  *
3-4:	Twist		*
4-5:	Gyrase_unbind
5-6:	Gyrase_bind	*
6-7:	TopoI_unbind
7-8:	TopoI_bind	*
8-9:	RNAP_tmp
9-10:	step
10+0:	mRNA1	
10+1:	protein1
'''
##############################################################
'''
num=160;
barrier=[10,100];
promoter_index=[30,60,120];
terminator_index=[50,80,140];
'''
'''
num=79;
barrier=[10,70];
promoter_index=[30];
terminator_index=[50];
'''
num=149;
barrier1=0;
barrier2=149;
#effective_len=149;

filename = sys.argv[1];
fp = h5py.File(filename, "r");
replicate=int(sys.argv[2]);
counts=fp["/Simulations/%07d/SpeciesCounts"%replicate];

t=0;
#times=100;
#SC_count=np.zeros((79,times));
times=200;
SC_count=np.zeros((149,times));
lk0 = 60.0;
ini_event = 0

j = 0;
while (t<=2000):
	if (counts[t, 13*num+1] > 0):
		ini_event = 1;
	if (ini_event == 1 and j < times and t >= 750):
		print(t)
#		SC = np.array(counts[t, range(3*num, 3*num+79)]);
		SC = np.array(counts[t, range(3*num, 3*num+149)]);
		SC_count[:,j] = (SC - lk0)/(1.0*lk0);
		j = j+1;
	elif (j >= times):
		break;
	t=t+1;

SC_count[(SC_count == -1).astype(bool)]=np.nan
SC_count=np.transpose(SC_count)

matplotlib.rcParams.update({"figure.figsize": (6,6)})
cmap = matplotlib.cm.coolwarm
cmap.set_bad(color='black')
plt.imshow(SC_count, cmap=cmap, interpolation='nearest', vmin=-0.5, vmax=0.5, extent=[1,149,950,750])
plt.xlabel('Position (bp)')
#plt.xlabel('Position (segment)')
plt.xticks([50,100], ['3,000','6,000'])

plt.ylabel('Time (s)')
plt.title('Initiation rate='+str(sys.argv[3]))
clb = plt.colorbar()
clb.set_label(r'$\sigma$')

plt.tight_layout()
plt.savefig(sys.argv[1]+'_'+str(sys.argv[2])+'.pdf')
plt.close()
#'''

