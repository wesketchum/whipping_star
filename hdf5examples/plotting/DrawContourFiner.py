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
totpoints = 2500

# Set custom font family
font = {'family' : 'Avenir',
        'weight' : 'normal',
        'size'   : 12}
plt.rc('font', **font)

def makeHeatMap(x, y, arrays, xlabel, ylabel, title, out ):
        fig, ax = plt.subplots(1, 1)
        ax = sns.heatmap(arrays, linewidth=0.001)#, vmax=100, vmin = 0.23)
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
        CS3 = plt.contour(x, y, z2, [0.9], colors=data_colors[1], linestyles='--');
        CS4 = plt.contour(x, y, z3, [0.9], colors=data_colors[2], linestyles='-.');
        #By design, labels for contour levels are expected to be put inline, instead of in the legend box, so this needs to be add manually
        h1,_ = CS1.legend_elements()
        h2,_ = CS2.legend_elements()
        h3,_ = CS3.legend_elements()
        h4,_ = CS4.legend_elements()
        lines = [ h1[0], h2[0], h3[0], h4[0] ]#, h6[0] ]
        #lines = [ h1[0], h2[0]]#, h4[0], h3[0] ]#, h6[0] ]
        #labels = ['Wilks, 3 d.o.f', 'Minimizer' ]#, 'Katie\'s Minimizer']
        labels = ['Wilks', '100x25x25','200x25x25','400x25x25']#, 'Katie\'s Minimizer']
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

#open files
f0 = h5.File("test_scan.h5",mode="r") 
f1 = h5.File("testmin_0to2500_1000univs.h5",mode="r") 
#f2 = h5.File("testnewmin_200x25x25_0to15625_iter50_1000univs.h5",mode="r") 
#f3 = h5.File("testnewmin_400x25x25_0to15625_iter50_1000univs.h5",mode="r") 
#f = h5.File("testnewmin_MFAmodel_400x50x50_0to2500_iter200.h5",mode="r") 
#f = h5.File("testminimizernewgrad_800x50x50_0to2500_iter100.h5",mode="r") 
#f = h5.File("testminimizernewgrad_800x50x50_125000to127500_iter100.h5",mode="r") 
#delta_chi
deltachi0 = f0['delta_chi'][:2500]
deltachi1 = f1['delta_chi'][:]
#deltachi2 = f2['delta_chi'][:]
#deltachi3 = f3['delta_chi'][:]
print("deltachi1.shape(): ", np.array(deltachi1))

delta_chi0 = np.array(deltachi0).reshape(totpoints,1)
delta_chi1 = np.array(deltachi1).reshape(totpoints,1000)
#delta_chi2 = np.array(deltachi2).reshape(15625,1000)
#delta_chi3 = np.array(deltachi3).reshape(15625,1000)
#chisens
chi2_arr=[]
#with open("chi2values100by50by50/chis_sens.txt","r") as chis:
with open("chi2values100by50by50/chis_sens.txt","r") as chis:
    lines = chis.readlines()
    for line in lines:
        chi2_arr.append(float(line))  

minpt=1000000

for pt in range(len(chi2_arr)):
    if chi2_arr[pt]<minpt:
        minpt = chi2_arr[pt]
        bestpt = pt
print("min, pt: ",minpt, bestpt)


cl_wilks=[]
cl_min1=[]
cl_min2=[]
cl_min3=[]
num=0
#pvals
#for i in range(122500,125000):
for i in range(0,totpoints):
    #print("i: ",i)
    #chi2_sens = np.array(chi2_arr[i])-minpt
    # make an R plot for each method
    #print("delta_chi[num].shape()  ", np.array(delta_chi).shape )
    #pvals_min=sum( np.array(delta_chi[num,:])  > chi2_sens )/1000
    pvals_min1=sum( np.array(delta_chi1[i,:])  > deltachi0[i] )/1000
    #pvals_min2=sum( np.array(delta_chi2[i,:])  > deltachi0[i] )/1000
    #pvals_min3=sum( np.array(delta_chi3[i,:])  > deltachi0[i] )/1000
    if 1.0-pvals_min1 < 0.24:
        print("grid: ", i)
    #cl_wilks.append(chi2_sens)
    if i != 0 and (i%49 == 0 or int(i/49) == 49):
        continue;
    cl_wilks.append(deltachi0[i])
    cl_min1.append(1.0-pvals_min1)
    #cl_min2.append(1.0-pvals_min2)
    #cl_min3.append(1.0-pvals_min3)
    #print("chi2_sens,cl_min: ", chi2_sens,1.0-pvals_min)
    num+=1


cls_wilks_vec=np.array(cl_wilks, dtype=float).reshape(1,49,49)
cls_min_vec1=np.array(cl_min1, dtype=float).reshape(1,49,49)
#cls_min_vec2=np.array(cl_min2, dtype=float).reshape(25,25,25)
#cls_min_vec3=np.array(cl_min3, dtype=float).reshape(25,25,25)

print("cls_min_vec1: ", cls_min_vec1[0,:,:] )

#get binning
#load in the bins for the grid points
#with open("chi2values100by50by50/bins_sens.txt","r") as bins:
with open("SBNfitoutputs/bins_sens.txt","r") as bins:
    lines = bins.readlines()
    lines[0].strip()
    lines[1].strip()
    lines[2].strip()
    vals0 = lines[0].split(" ")
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
#print("m values:", m41_list )
#print("ue4 values:",ue4_list)
#print("um4 values:",um4_list)

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
#y,x = np.meshgrid(m41_center, ue4_center)
#makeContourStandard(x, y, [.116],[1.16], cls_wilks_vec[:,:,17].transpose(), cls_min_vec1[:,:,17].transpose(), cls_min_vec2[:,:,17].transpose(), cls_min_vec3[:,:,17].transpose(), "$U_{e4}$","$\Delta m^2_{41}$", "$U_{\mu4}$ ="+str(um4_list[17]), "FCcontours_minimizer_vs_gridscan_umu4_"+str(um4_list[17])+"_25x25x25")
#makeContourStandard(x, y, [.116],[1.16], cls_wilks_vec[:,:,17].transpose(), cls_min_vec[:,:,17].transpose(), cls_minimizer_vec[:,:,17].transpose(), cls_minimizer_katies_vec[:,:,17].transpose(), "$U_{e4}$","$\Delta m^2_{41}$", "$U_{\mu4}$ ="+str(um4_list[17]), "FCcontours_minimizer_vs_gridscan_umu4_"+str(um4_list[17])+"_25x25x25")
#makeHeatMap(x, y, cls_min_vec[:,:,17],  "$U_{e4}$","$\Delta m^2_{41}$", "$U_{\mu4}$ ="+str(um4_list[17]), "FCpvalsheatmap_gridscan_umu4_"+str(um4_list[17])+"_25x25x25")
#makeHeatMap(x, y, cls_minimizer_vec[:,:,17],  "$U_{e4}$","$\Delta m^2_{41}$", "$U_{\mu4}$ ="+str(um4_list[17]), "FCpvalsheatmap_minimizer_umu4_"+str(um4_list[17])+"25x25x25")
#makeHeatMap(x, y, cls_minimizer_katies_vec[:,:,17],  "$U_{e4}$","$\Delta m^2_{41}$", "$U_{\mu4}$ ="+str(um4_list[17]), "FCpvalsheatmap_minimizer_katies_umu4_"+str(um4_list[17])+"25x25x25")
#
#### Ue4
#y,x = np.meshgrid(m41_center, um4_center)
#makeContourStandard(x, y, [.135],[1.16], cls_wilks_vec[:,16,:].transpose(), cls_min_vec1[:,16,:].transpose(), cls_min_vec2[:,16,:].transpose(), cls_min_vec3[:,16,:].transpose(), "$U_{\mu 4}$","$\Delta m^2_{41}$","$U_{e4}$ ="+str(ue4_list[16]), "FCcontours_minimizer_vs_gridscan_ue4_"+str(ue4_list[16])+"_25x25x25")
#makeContourStandard(x, y, [.135],[1.16], cls_wilks_vec[:,16,:].transpose(), cls_min_vec[:,16,:].transpose(), cls_minimizer_vec[:,16,:].transpose(), cls_minimizer_katies_vec[:,16,:].transpose(), "$U_{\mu 4}$","$\Delta m^2_{41}$","$U_{e4}$ ="+str(ue4_list[16]), "FCcontours_minimizer_vs_gridscan_ue4_"+str(ue4_list[16])+"_25x25x25")
#makeHeatMap(x, y, cls_min_vec[:,16,:], "$U_{\mu 4}$","$\Delta m^2_{41}$","$U_{e4}$ ="+str(ue4_list[16]), "FCpvalsheatmap_gridscan_ue4_"+str(ue4_list[16])+"_25x25x25")
#makeHeatMap(x, y, cls_minimizer_vec[:,16,:], "$U_{\mu 4}$","$\Delta m^2_{41}$","$U_{e4}$ ="+str(ue4_list[16]), "FCpvalsheatmap_minimizer_ue4_"+str(ue4_list[16])+"_25x25x25")
#makeHeatMap(x, y, cls_minimizer_katies_vec[:,16,:], "$U_{\mu 4}$","$\Delta m^2_{41}$","$U_{e4}$ ="+str(ue4_list[16]), "FCpvalsheatmap_minimizer_katies_ue4_"+str(ue4_list[16])+"_25x25x25")

#for idx in [0,5,10,15,20,24]:
for idx in [0]:
    ## delta msq
    x,y = np.meshgrid(ue4_center[:49], um4_center[:49])
    #makeContourStandard(x, y, [.116],[.135], cls_wilks_vec[idx,:,:], cls_min_vec1[idx,:,:], cls_min_vec2[idx,:,:], cls_min_vec3[idx,:,:], "$U_{\mu 4}$", "$U_{e4}$", "$\Delta m^2_{41}$ ="+str(m41_list[idx]), "FCcontours_minimizer_vs_gridscan_m41_"+str(m41_list[idx])+"_25x25x25")
    #makeContourStandard(x, y, [.116],[.135], cls_wilks_vec[0,:,:], cls_min_vec[0,:,:], "$U_{\mu 4}$", "$U_{e4}$", "$\Delta m^2_{41}$ ="+str(m41_list[idx]), "FCcontours_minimizer_vs_gridscan_m41_"+str(m41_list[idx])+"_100x50x50")
    makeHeatMap(x, y, cls_min_vec1[0,:,:], "$U_{\mu 4}$", "$U_{e4}$", "$\Delta m^2_{41}$ ="+str(0.01), "FCpvalsheatmap_gridscan_m41_MFAmodel400x50x50_1x49x49_iter50")
    #makeHeatMap(x, y, cls_minimizer_vec[idx,:,:], "$U_{\mu 4}$", "$U_{e4}$", "$\Delta m^2_{41}$ ="+str(m41_list[idx]), "FCpvalsheatmap_minimizer_MFAmodel_100x25x25_iter50_m41_"+str(m41_list[idx])+"_25x25x25")
    #makeHeatMap(x, y, cls_minimizer_katies_vec[idx,:,:], "$U_{\mu 4}$", "$U_{e4}$", "$\Delta m^2_{41}$ ="+str(m41_list[idx]), "FCpvalsheatmap_minimizer_katies_m41_"+str(m41_list[idx])+"_25x25x25")
