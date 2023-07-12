# Imports
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
import scipy.interpolate
from matplotlib import cm
import scipy.stats
from matplotlib import gridspec
import matplotlib.lines as mlines
import scipy.stats
from matplotlib import gridspec
import matplotlib.lines as mlines
from scipy.stats import chi2
from scipy import interpolate
import pandas as pd
import math
import matplotlib.colors as colors
from os.path import exists
import pickle as pkl
from PIL import Image
import glob
import time
import seaborn as sns


# Set custom colors and linestyles
data_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
fit_color = '#e31a1c'   # red
quad_color = '#1f78b4'  # blue
data_style = 'o'
fit_style = '-'
quad_style = '-'

# Set custom font family
font = {'family' : 'Avenir',
        'weight' : 'normal',
        'size'   : 12}
plt.rc('font', **font)

def makeHeatMap(x, y, arrays, xlabel, ylabel, title, out ):
        fig, ax = plt.subplots(1, 1)
        ax = sns.heatmap(arrays, linewidth=0.001)
        plt.xlabel(xlabel,fontsize=20)
        plt.ylabel(ylabel,fontsize=20)
        plt.title(title,fontsize=20)
        ax.invert_yaxis()
        #plt.xscale('log')
        #plt.yscale('log')
        plt.grid()
        #leg=ax.legend(lines,labels,loc="lower left",fontsize=8,title='Sensitivity, 90% CL');
        #plt.legend(loc="lower right")
        fig.savefig(out+'.pdf',transparent=True)

def makeContourStandard(x, y, globfix, globfity, z0, z1, z2, z3, xlabel, ylabel, title, out):
        #Make contour plots
        fig, ax = plt.subplots(1, 1)
        CS1 = plt.contour(x, y, z0, [6.25], colors='darkgray', linestyles='-.');
        CS2 = plt.contour(x, y, z1, [0.9], colors=data_colors[0]);
        CS4 = plt.contour(x, y, z3, [0.9], colors=data_colors[2], linestyles='-.');
        CS3 = plt.contour(x, y, z2, [0.9], colors=data_colors[1], linestyles='--');
        #By design, labels for contour levels are expected to be put inline, instead of in the legend box, so this needs to be add manually
        h1,_ = CS1.legend_elements()
        h2,_ = CS2.legend_elements()
        h4,_ = CS4.legend_elements()
        h3,_ = CS3.legend_elements()
        lines = [ h1[0], h2[0], h4[0], h3[0] ]#, h6[0] ]
        labels = ['Wilks, 3 d.o.f', 'Grid-Scan (inner grid 25x25x25)', 'Grid-Scan (inner grid 13x13x13)', 'Minimizer' ]#, 'Katie\'s Minimizer']
        plt.xlabel(xlabel,fontsize=20)
        plt.ylabel(ylabel,fontsize=20)
        plt.title(title,fontsize=20)
        #if( ylabel=="$\Delta m^2_{41}$"):
        #        plt.ylim(0.01, 100.)
        #else: 
        #        plt.ylim(0.01, 1.0)
        #plt.xlim(0.01, 1.)
        plt.xscale('log')
        plt.yscale('log')
        plt.grid()
        # Set the padding for the bottom of the plot
        plt.subplots_adjust(bottom=0.15)
        plt.subplots_adjust(left=0.15)
        leg=ax.legend(lines,labels,loc="lower left",fontsize=8,title='Sensitivity, 90% CL');
        #plt.legend(loc="lower right")
        # add labels to the left and right of the plot
        ax.text(0.96, 0.96, 'excluded\nregion', transform=ax.transAxes, fontsize=12, va='top', ha='right')# bbox=dict(boxstyle='square', facecolor='white', alpha=0.8))
        ax.text(0.38, 0.55, 'allowed\n region', transform=ax.transAxes, fontsize=12, va='top', ha='left')# bbox=dict(boxstyle='square', facecolor='white', alpha=0.8))
        fig.savefig(out+'.pdf',transparent=True)
        plt.close()

start = time.time()
# initialize arr to fill
# set number of universes
# look at each index individually
num=0
# also initialize 3d array:
cls_wilks_vec=np.zeros((25, 25, 25))
cls_min_vec=np.zeros((25, 25, 25))
cls_minimizer_vec=np.zeros((25, 25, 25))
cls_minimizer_katies_vec=np.zeros((25, 25, 25))

cl_wilks=[]
cl_min=[]
pvals_min=[]
cl_minimizer=[]
pvals_minimizer=[]
cl_minimizer_katies=[]
pvals_minimizer_katies=[]

with open('cl_wilks_25x25x25.pkl', 'rb') as f4:
    cl_wilks = pkl.load(f4)
with open('cl_minimizergridscanHPC15625pts.pkl', 'rb') as f5:
    cl_min = pkl.load(f5)
with open('pvals_mingridscanHPC15625pts.pkl', 'rb') as f6:
    pvals_min = pkl.load(f6)
with open('cl_minimizertestminimizernewgrad_iter50_0to15625.pkl', 'rb') as f7:
    cl_minimizer = pkl.load(f7)
with open('pvals_mintestminimizernewgrad_iter50_0to15625.pkl', 'rb') as f8:
    pvals_min_minimizer = pkl.load(f8)
#with open('cl_minimizertestminimizernewgrad_100x25x25_0to15625_iter50.pkl', 'rb') as f7:
#    cl_minimizer = pkl.load(f7)
#with open('pvals_mintestminimizernewgrad_100x25x25_0to15625_iter50.pkl', 'rb') as f8:
#    pvals_min_minimizer = pkl.load(f8)
#with open('cl_minimizertestminimizernewgrad_iter50_0tO15625_TRY2.pkl', 'rb') as f7:
#    cl_minimizer = pkl.load(f7)
#with open('pvals_mintestminimizernewgrad_iter50_0tO15625_TRY2.pkl', 'rb') as f8:
#    pvals_min_minimizer = pkl.load(f8)
#with open('cl_minimizertestminimizernewgrad_iter50_0tO15625_TRY2.pkl', 'rb') as f7:
#    cl_minimizer = pkl.load(f7)
#with open('pvals_mintestminimizernewgrad_iter50_0tO15625_TRY2.pkl', 'rb') as f8:
#    pvals_min_minimizer = pkl.load(f8)
with open('cl_minimizercoarseinnerscan.pkl', 'rb') as f7:
    cl_minimizer_katies = pkl.load(f7)
with open('pvals_mincoarseinnerscan.pkl', 'rb') as f8:
    pvals_min_minimizer_katies = pkl.load(f8)

cls_wilks_vec=np.array(cl_wilks, dtype=float).reshape((25, 25, 25))
cls_min_vec=np.array(cl_min, dtype=float).reshape((25, 25, 25))
cls_minimizer_vec=np.array(cl_minimizer, dtype=float).reshape((25, 25, 25))
cls_minimizer_katies_vec=np.array(cl_minimizer_katies, dtype=float).reshape((25, 25, 25))

print("cls wilks um4: ", len([ a for v in cls_wilks_vec[:,:,17] for a in v if a > 6.25 ]))
print("cls grid scan um4: ", len([ a for v in cls_min_vec[:,:,17] for a in v if a > 0.9 ]))
print("cls minimizer um4: ", len([ a for v in cls_minimizer_vec[:,:,17] for a in v if a > 0.9 ]))
print("cls minimizer katies um4: ", len([ a for v in cls_minimizer_katies_vec[:,:,17] for a in v if a > 0.9 ]))

print("cls wilks ue4: ", len([ a for v in cls_wilks_vec[:,16,:] for a in v if a > 6.25 ]))
print("cls grid scan ue4: ", len([ a for v in cls_min_vec[:,16,:] for a in v if a > 0.9 ]))
print("cls minimizer ue4: ", len([ a for v in cls_minimizer_vec[:,16,:] for a in v if a > 0.9 ]))
print("cls minimizer katies ue4: ", len([ a for v in cls_minimizer_katies_vec[:,16,:] for a in v if a > 0.9 ]))

print("cls wilks delta m2: ", len([ a for v in cls_wilks_vec[0,:,:] for a in v if a > 6.25 ]))
print("cls grid scan delta m2: ", len([ a for v in cls_min_vec[0,:,:] for a in v if a > 0.9 ]))
print("cls minimizer delta m2: ", len([ a for v in cls_minimizer_vec[0,:,:] for a in v if a > 0.9 ]))
print("cls minimizer katies delta m2: ", len([ a for v in cls_minimizer_katies_vec[0,:,:] for a in v if a > 0.9 ]))

#get binning
#load in the bins for the grid points
with open("SBNfitoutputs/bins_sens.txt","r") as bins:
    lines = bins.readlines()
    lines[0].strip()
    lines[1].strip()
    lines[2].strip()
    vals0 = lines[0].split(" ")
    print("vals0: ",vals0)
    vals1 = lines[1].split(" ")
    vals2 = lines[2].split(" ")
    m41_list= []
    ue4_list= []
    um4_list= []
    for v in range(len(vals0)-1):
        m41_list.append(float(vals0[v]))
    for v in range(len(vals1)-1):
        ue4_list.append(float(vals1[v]))
    for v in range(len(vals2)-1):
        um4_list.append(float(vals2[v]))
print("m values:", m41_list )
print("ue4 values:",ue4_list)
print("um4 values:",um4_list)

# switch to bin centers for plotting
def getbincenters( fulllist):
    newlist =[]
    for i in range(len(fulllist)-1):
        topval = fulllist[i+1]
        botval = fulllist[i]
        center = (topval+botval)/float(2.0)
        newlist.append(center)
    return newlist

m41_center = getbincenters(m41_list)
ue4_center = getbincenters(ue4_list)
um4_center = getbincenters(um4_list)

stop = time.time()
print("elapsed time: ", stop-start)

## draw slices for best fit point
## Umu4
y,x = np.meshgrid(m41_center, ue4_center)
makeContourStandard(x, y, [.116],[1.16], cls_wilks_vec[:,:,17].transpose(), cls_min_vec[:,:,17].transpose(), cls_minimizer_vec[:,:,17].transpose(), cls_minimizer_katies_vec[:,:,17].transpose(), "$U_{e4}$","$\Delta m^2_{41}$", "$U_{\mu4}$ ="+str(um4_list[17]), "FCcontours_minimizer_vs_gridscan_umu4_"+str(um4_list[17])+"_25x25x25")
makeHeatMap(x, y, cls_min_vec[:,:,17],  "$U_{e4}$","$\Delta m^2_{41}$", "$U_{\mu4}$ ="+str(um4_list[17]), "FCpvalsheatmap_gridscan_umu4_"+str(um4_list[17])+"_25x25x25")
makeHeatMap(x, y, cls_minimizer_vec[:,:,17],  "$U_{e4}$","$\Delta m^2_{41}$", "$U_{\mu4}$ ="+str(um4_list[17]), "FCpvalsheatmap_minimizer_umu4_"+str(um4_list[17])+"25x25x25")
makeHeatMap(x, y, cls_minimizer_katies_vec[:,:,17],  "$U_{e4}$","$\Delta m^2_{41}$", "$U_{\mu4}$ ="+str(um4_list[17]), "FCpvalsheatmap_minimizer_katies_umu4_"+str(um4_list[17])+"25x25x25")

### Ue4
y,x = np.meshgrid(m41_center, um4_center)
makeContourStandard(x, y, [.135],[1.16], cls_wilks_vec[:,16,:].transpose(), cls_min_vec[:,16,:].transpose(), cls_minimizer_vec[:,16,:].transpose(), cls_minimizer_katies_vec[:,16,:].transpose(), "$U_{\mu 4}$","$\Delta m^2_{41}$","$U_{e4}$ ="+str(ue4_list[16]), "FCcontours_minimizer_vs_gridscan_ue4_"+str(ue4_list[16])+"_25x25x25")
makeHeatMap(x, y, cls_min_vec[:,16,:], "$U_{\mu 4}$","$\Delta m^2_{41}$","$U_{e4}$ ="+str(ue4_list[16]), "FCpvalsheatmap_gridscan_ue4_"+str(ue4_list[16])+"_25x25x25")
makeHeatMap(x, y, cls_minimizer_vec[:,16,:], "$U_{\mu 4}$","$\Delta m^2_{41}$","$U_{e4}$ ="+str(ue4_list[16]), "FCpvalsheatmap_minimizer_ue4_"+str(ue4_list[16])+"_25x25x25")
makeHeatMap(x, y, cls_minimizer_katies_vec[:,16,:], "$U_{\mu 4}$","$\Delta m^2_{41}$","$U_{e4}$ ="+str(ue4_list[16]), "FCpvalsheatmap_minimizer_katies_ue4_"+str(ue4_list[16])+"_25x25x25")

for idx in [13]:
    ## delta msq
    x,y = np.meshgrid(ue4_center, um4_center)
    makeContourStandard(x, y, [.116],[.135], cls_wilks_vec[idx,:,:], cls_min_vec[idx,:,:], cls_minimizer_vec[idx,:,:], cls_minimizer_katies_vec[idx,:,:], "$U_{\mu 4}$", "$U_{e4}$", "$\Delta m^2_{41}$ ="+str(m41_list[idx]), "FCcontours_minimizer_vs_gridscan_m41_"+str(m41_list[idx])+"_25x25x25")
    makeHeatMap(x, y, cls_min_vec[idx,:,:], "$U_{\mu 4}$", "$U_{e4}$", "$\Delta m^2_{41}$ ="+str(m41_list[idx]), "FCpvalsheatmap_gridscan_m41_"+str(m41_list[idx])+"_25x25x25")
    makeHeatMap(x, y, cls_minimizer_vec[idx,:,:], "$U_{\mu 4}$", "$U_{e4}$", "$\Delta m^2_{41}$ ="+str(m41_list[idx]), "FCpvalsheatmap_minimizer_MFAmodel_100x25x25_iter50_m41_"+str(m41_list[idx])+"_25x25x25")
    makeHeatMap(x, y, cls_minimizer_katies_vec[idx,:,:], "$U_{\mu 4}$", "$U_{e4}$", "$\Delta m^2_{41}$ ="+str(m41_list[idx]), "FCpvalsheatmap_minimizer_katies_m41_"+str(m41_list[idx])+"_25x25x25")
