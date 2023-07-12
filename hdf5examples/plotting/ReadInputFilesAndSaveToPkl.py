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



folder,nindices,univ,outputlabel = input("Enter folder name, number of gridpoints, number of universes, output name for pkl file: ").split()

# get folder name
folder=str(folder)
# set number of indices
nindices=int(nindices)
# set number of universes
num_uni=int(univ)
# label for output file
outputlabel = str(outputlabel)

start = time.time()

indices=[i for i in range(nindices)]

# reformat indices to string format used to save files
for i in range(nindices): 
    indices[i]=f'{indices[i]:05d}'

# -2lnL dicts for test point, bf from grid, and bf from min
RDict_pt = {}
RDict_min = {}

# loop through text files and save to dict
for k in indices:
    filled=True
    pT_list=[]
    min_list=[]
    #print("folder: ",folder)
    if exists(folder+"/textfiles/chis_seek_"+k+".txt"):
    #if exists(folder+"/chis_seek_"+k+".txt"):
        #print("folder exists")
        with open(folder+"/textfiles/chis_seek_"+k+".txt","r") as chi:
        #with open(folder+"/chis_seek_"+k+".txt","r") as chi:
            lines = chi.readlines()
            if(len(lines)<2000):
                print(folder+"/textfiles/chis_seek_"+k+".txt",len(lines))
                #print(folder+"/chis_seek_"+k+".txt",len(lines))
            for i in range(int(len(lines)/2)):
                lines[i*2].strip() #pt
                lines[i*2+1].strip() #min

                vals0 = lines[i*2].split(" ")
                vals1 = lines[i*2+1].split(" ")
                #print("vals0, vals1: ", vals0, vals1)
                for j in range(len(vals0)):
                    pT_list.append(float(vals0[j]))
                for j in range(len(vals1)):
                    min_list.append(float(vals1[j]))
    else: #default to wilks theorem for now if job failed
        for j in range(num_uni):
            min_list.append(float(6.25))
        for j in range(num_uni):
            pT_list.append(float(6.25))

#   save to dictionary        
    RDict_pt[k] = pT_list
    RDict_min[k] = min_list

with open('Rdict_pt'+outputlabel+'.pkl', 'wb') as f1:
    pkl.dump(RDict_pt, f1)
with open('Rdict_min'+outputlabel+'.pkl', 'wb') as f2:
    pkl.dump(RDict_min, f2)
