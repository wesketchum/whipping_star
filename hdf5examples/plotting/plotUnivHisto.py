import h5py as h5
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import seaborn as sns


startpoint=15021
totpoints=1
nuniv=1000
pt=0
#filename1="testminiminew_400x25x25_0to5_10iters2.h5"
#filename2="testminiminew_400x25x25_0to5_20iters.h5"
#filename3="testminiminew_400x25x25_0to5_50iters2.h5"
#filename4="test3d_FC_1000univ_gp0to5.h5"
filename1="testminiminew_400x25x25_gp"+str(startpoint)+"_10iters.h5"
filename2="testminiminew_400x25x25_gp"+str(startpoint)+"_20iters.h5"
filename3="testminiminew_400x25x25_gp"+str(startpoint)+"_50iters.h5"
filename4="test3d_FC_1000univ_gp"+str(startpoint)+".h5"
f0 = h5.File("test_scan2.h5", mode="r")
f1 = h5.File(filename1, mode="r")
f2 = h5.File(filename2, mode="r")
f3 = h5.File(filename3, mode="r")
f4 = h5.File(filename4, mode="r")

startx1_ = f1['best_start_pointx'][:totpoints*nuniv]
starty1_ = f1['best_start_pointy'][:totpoints*nuniv]
startz1_ = f1['best_start_pointz'][:totpoints*nuniv]
bestx1_ = f1['best_grid_pointx'][:totpoints*nuniv]
besty1_ = f1['best_grid_pointy'][:totpoints*nuniv]
bestz1_ = f1['best_grid_pointz'][:totpoints*nuniv]
delta1_ = f1['delta_chi'][:totpoints*nuniv] 
minimum1_ = f1['last_chi_min'][:totpoints*nuniv] 
data1_ = np.add(delta1_,minimum1_)
delta0 = f0['delta_chi'][:]
minimum0 = f0['last_chi_min'][:]
data0 = np.add(delta0,minimum0)

startx2_ = f2['best_start_pointx'][:totpoints*nuniv]
starty2_ = f2['best_start_pointy'][:totpoints*nuniv]
startz2_ = f2['best_start_pointz'][:totpoints*nuniv]
bestx2_ = f2['best_grid_pointx'][:totpoints*nuniv]
besty2_ = f2['best_grid_pointy'][:totpoints*nuniv]
bestz2_ = f2['best_grid_pointz'][:totpoints*nuniv]
delta2_ = f2['delta_chi'][:totpoints*nuniv] 
minimum2_ = f2['last_chi_min'][:totpoints*nuniv] 
data2_ = np.add(delta2_,minimum2_)

startx3_ = f3['best_start_pointx'][:totpoints*nuniv]
starty3_ = f3['best_start_pointy'][:totpoints*nuniv]
startz3_ = f3['best_start_pointz'][:totpoints*nuniv]
bestx3_ = f3['best_grid_pointx'][:totpoints*nuniv]
besty3_ = f3['best_grid_pointy'][:totpoints*nuniv]
bestz3_ = f3['best_grid_pointz'][:totpoints*nuniv]
delta3_ = f3['delta_chi'][:totpoints*nuniv] 
minimum3_ = f3['last_chi_min'][:totpoints*nuniv] 
data3_ = np.add(delta3_,minimum3_)

delta4_ = f4['delta_chi'][:totpoints*nuniv] 
minimum4_ = f4['last_chi_min'][:totpoints*nuniv] 
data4_ = np.add(delta4_,minimum4_)

#now group based on the gridpoints:
startx1 = startx1_.reshape(totpoints,nuniv)
starty1 = starty1_.reshape(totpoints,nuniv)
startz1 = startz1_.reshape(totpoints,nuniv)
bestx1 = bestx1_.reshape(totpoints,nuniv)
besty1 = besty1_.reshape(totpoints,nuniv)
bestz1 = bestz1_.reshape(totpoints,nuniv)
delta1 = delta1_.reshape(totpoints,nuniv)
minimum1 = minimum1_.reshape(totpoints,nuniv)
data1 = data1_.reshape(totpoints,nuniv)

startx2 = startx2_.reshape(totpoints,nuniv)
starty2 = starty2_.reshape(totpoints,nuniv)
startz2 = startz2_.reshape(totpoints,nuniv)
bestx2 = bestx2_.reshape(totpoints,nuniv)
besty2 = besty2_.reshape(totpoints,nuniv)
bestz2 = bestz2_.reshape(totpoints,nuniv)
delta2 = delta2_.reshape(totpoints,nuniv)
minimum2 = minimum2_.reshape(totpoints,nuniv)
data2 = data2_.reshape(totpoints,nuniv)

startx3 = startx3_.reshape(totpoints,nuniv)
starty3 = starty3_.reshape(totpoints,nuniv)
startz3 = startz3_.reshape(totpoints,nuniv)
bestx3 = bestx3_.reshape(totpoints,nuniv)
besty3 = besty3_.reshape(totpoints,nuniv)
bestz3 = bestz3_.reshape(totpoints,nuniv)
delta3 = delta3_.reshape(totpoints,nuniv)
minimum3 = minimum3_.reshape(totpoints,nuniv)
data3 = data3_.reshape(totpoints,nuniv)

delta4 = delta4_.reshape(totpoints,nuniv)
minimum4 = minimum4_.reshape(totpoints,nuniv)
data4 = data4_.reshape(totpoints,nuniv)

data1_  = []
delta1_ = []
bestx1_ = []
besty1_ = []
bestz1_ = []
startx1_ = []
starty1_ = []
startz1_ = []

data2_  = []
delta2_ = []
bestx2_ = []
besty2_ = []
bestz2_ = []
startx2_ = []
starty2_ = []
startz2_ = []

data3_  = []
delta3_ = []
bestx3_ = []
besty3_ = []
bestz3_ = []
startx3_ = []
starty3_ = []
startz3_ = []

data4_  = []
delta4_ = []

i=0
tot=0
for pt in range(totpoints):
    if bestx1[pt][i] <= 384.0 and besty1[pt][i] <= 48.0 and bestz1[pt][i] <= 48.0 and starty1[pt][i] <= 48.0 and startz1[pt][i] <= 48.0:
        tot += 1
        data1_.append(minimum1[pt][i])
        startx1_.append(startx1[pt][i])
        starty1_.append(starty1[pt][i])
        startz1_.append(startz1[pt][i])
        bestx1_.append(bestx1[pt][i])
        besty1_.append(besty1[pt][i])
        bestz1_.append(bestz1[pt][i])
        delta1_.append(delta1[pt][i])
        data2_.append(minimum2[pt][i])
        startx2_.append(startx2[pt][i])
        starty2_.append(starty2[pt][i])
        startz2_.append(startz2[pt][i])
        bestx2_.append(bestx2[pt][i])
        besty2_.append(besty2[pt][i])
        bestz2_.append(bestz2[pt][i])
        delta2_.append(delta2[pt][i])
        data3_.append(minimum3[pt][i])
        startx3_.append(startx3[pt][i])
        starty3_.append(starty3[pt][i])
        startz3_.append(startz3[pt][i])
        bestx3_.append(bestx3[pt][i])
        besty3_.append(besty3[pt][i])
        bestz3_.append(bestz3[pt][i])
        delta3_.append(delta3[pt][i])
        data4_.append(minimum4[pt][i])
        delta4_.append(delta4[pt][i])
print("finish loop, 400x400, how many minimum,percentage: ", len(bestx1_)/tot)
print("len: ", len(bestx1_))
print("tot: ", tot)
#print("data_: ", data_)


for pt,point in enumerate(range(startpoint,startpoint+totpoints)):
    fig, ax = plt.subplots(1, 1)
    #print(data1[pt])
    plt.hist(np.array(data0[pt]),bins=range(100, 500, 5), histtype='step', edgecolor='black', label='Wilks')
    plt.hist(np.array(data1[pt]),bins=range(100, 500, 5), histtype='step', edgecolor='b', label='10 iters')
    plt.hist(np.array(data2[pt]),bins=range(100, 500, 5), histtype='step', edgecolor='r', label='20 iters')
    plt.hist(np.array(data3[pt]),bins=range(100, 500, 5), histtype='step', edgecolor='g', label='50 iters')
    plt.hist(np.array(data4[pt]),bins=range(100, 500, 5), histtype='step', edgecolor='orange', label='FC grid scan')
    legend = ax.legend(loc='upper right', shadow=True)
    ax.set_xlabel(r'$\chi^2$')
    #ax.set_ylabel('volts')
    ax.set_title(r'$\chi^2$, gp = '+ str(point))
    fig.savefig('histo_this_chi_gp'+str(point)+'.pdf',transparent=True)
    plt.close()

    fig, ax = plt.subplots(1, 1)
    plt.hist(np.array(minimum0[pt]),bins=range(100, 500, 5), histtype='step', edgecolor='black', label='Wilks')
    plt.hist(np.array(minimum1[pt]),bins=range(100, 500, 5), histtype='step', edgecolor='b', label='10 iters')
    plt.hist(np.array(minimum2[pt]),bins=range(100, 500, 5), histtype='step', edgecolor='r', label='20 iters')
    plt.hist(np.array(minimum3[pt]),bins=range(100, 500, 5), histtype='step', edgecolor='r', label='50 iters')
    plt.hist(np.array(minimum4[pt]),bins=range(100, 500, 5), histtype='step', edgecolor='orange', label='FC grid scan')
    legend = ax.legend(loc='upper right', shadow=True)
    ax.set_xlabel(r'min $\chi^2$')
    #ax.set_ylabel('volts')
    ax.set_title(r'minimum $\chi^2$, gp = '+ str(point))
    fig.savefig('histo_minchi_gp'+str(point)+'.pdf',transparent=True)
    plt.close()

    fig, ax = plt.subplots(1, 1)
    #print(delta1[pt]-)
    plt.hist(np.array(delta1[pt]),bins=range(-10, 90, 2), histtype='step', edgecolor='b', label='10 iters, $R_{90\%}$ = ' + str(np.round(np.percentile(np.array(delta1[pt]), 90),2)))
    plt.hist(np.array(delta2[pt]),bins=range(-10, 90, 2), histtype='step', edgecolor='r', label='20 iters, $R_{90\%}$ = ' + str(np.round(np.percentile(np.array(delta2[pt]), 90),2)))
    plt.hist(np.array(delta3[pt]),bins=range(-10, 90, 2), histtype='step', edgecolor='g', label='50 iters, $R_{90\%}$ = ' + str(np.round(np.percentile(np.array(delta3[pt]), 90),2)))
    plt.hist(np.array(delta4[pt]),bins=range(-10, 90, 2), histtype='step', edgecolor='orange', label='FC grid scan, $R_{90\%}$ = ' + str(np.round(np.percentile(np.array(delta4[pt]), 90),2)))
    plt.plot([np.percentile(np.array(delta1[pt]), 90),np.percentile(np.array(delta1[pt]), 90)],[0,300], color='blue')
    plt.plot([np.percentile(np.array(delta2[pt]), 90),np.percentile(np.array(delta2[pt]), 90)],[0,300], color='red')
    plt.plot([np.percentile(np.array(delta3[pt]), 90),np.percentile(np.array(delta3[pt]), 90)],[0,300], color='green')
    plt.plot([np.percentile(np.array(delta4[pt]), 90),np.percentile(np.array(delta4[pt]), 90)],[0,300], color='orange')
    plt.plot([delta0[point],delta0[point]],[0,300], color='black', linestyle='dashed', label='Wilks, $R = $' + str(delta0[point]))
    legend = ax.legend(loc='upper right', shadow=True)
    ax.set_xlabel(r'R')
    #ax.set_ylabel('volts')
    ax.set_title(r'R, gp = '+ str(point))
    fig.savefig('histo_delta_gp'+str(point)+'.pdf',transparent=True)
    plt.close()
