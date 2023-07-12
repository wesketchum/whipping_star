import h5py as h5
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
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

startpoint=0
totpoints=500000
nuniv=1
pt=0
labels=["50x25x25","100x25x25","200x25x25","400x25x25","800x25x25"]
filename4="testmin_MFAresol800x25x25_input800mass.h5"
filename3="testmin_MFAresol400x25x25_input800mass.h5"
filename2="testmin_MFAresol200x25x25_input800mass.h5"
filename1="testmin_MFAresol100x25x25_input800mass.h5"
filename0="testmin_MFAresol50x25x25_input800mass.h5"
files=[]
files.append(h5.File(filename0, mode="r"))
files.append(h5.File(filename1, mode="r"))
files.append(h5.File(filename2, mode="r"))
files.append(h5.File(filename3, mode="r"))
files.append(h5.File(filename4, mode="r"))

chi_raw_=[]
chi_mfa_=[]
chi_gs_=[]
chi_raw=[]
chi_mfa=[]
chi_gs=[]
speccoll_gs_=[]
speccoll_gs=[]
speccoll_mfa_=[]
speccoll_mfa=[]
speccoll_diff=[]
for i in range(len(files)):
    chi_raw_.append(files[i]['this_chi_raw'][:totpoints*nuniv])
    chi_mfa_.append(files[i]['this_chi_mfa'][:totpoints*nuniv])
    chi_gs_.append(files[i]['this_chi_gs'][:totpoints*nuniv])
    speccoll_gs_.append(files[i]['speccoll_gs'][:totpoints*nuniv])
    speccoll_mfa_.append(files[i]['speccoll_mfa'][:totpoints*nuniv])
    #calculate chi2
    for p in range(totpoints*nuniv):
        chi_raw_[i][p] = np.sum(np.true_divide(np.square(np.subtract(np.array(speccoll_gs_[i][p]),np.array(speccoll_mfa_[i][p]))),np.array(speccoll_gs_[i][p])))
    #reshape to 800x25x25
    #chi_raw.append(np.divide(np.array(chi_raw_[i]),31).reshape(800,25,25))
    chi_raw.append(np.divide(np.array(chi_raw_[i]),31).reshape(800,25,25))
    chi_mfa.append(np.divide(np.array(chi_mfa_[i]),31).reshape(800,25,25))
    chi_gs.append(np.divide(np.array(chi_gs_[i]),31).reshape(800,25,25))
    speccoll_gs.append(np.array(speccoll_gs_[i]).reshape(800,25,25,31))
    speccoll_mfa.append(np.array(speccoll_mfa_[i]).reshape(800,25,25,31))
    speccoll_diff.append(np.subtract(np.array(speccoll_gs_[i]),np.array(speccoll_mfa_[i])).reshape(800,25,25,31))
    
print("chi_raw shape: ", chi_raw[0].shape)
print("speccoll_diff shape: ", speccoll_diff[0].shape)

#get binning
f0 = h5.File("testmin_MFAresol800x25x25_input800mass.h5",mode="r") 
gridx_ = f0['gridx'][:]
gridy_ = f0['gridy'][:]
gridz_ = f0['gridz'][:]
gridx = np.unique(gridx_)
gridy = np.unique(gridy_)
gridz = np.unique(gridz_)
#print("gridx = ", gridx)
print("gridy = ", gridy[15], gridy[16], gridz[17])
#print("gridz = ", gridz)
m41_list = [round(10**(2*x),4) for x in gridx]
ue4_list = [round(10**(y),4) for y in gridy]
um4_list = [round(10**(z),4) for z in gridz]
#print("m41_list: ",m41_list)
#print("ue4_list: ",ue4_list)
#print("um4_list: ",um4_list)


#wash out any points that we share with the true points
truebins=[32,16,8,4,2,1]
transparency=[[] for t in range(len(files))]
print("transparency: ", transparency)
for f in range(len(files)):
    for i in range(800):
        if i > 0 and i%truebins[f] == 0:
            transparency[f].append(0.0)
        else:
            transparency[f].append(0.8-f*0.1)
print("transparency[0]: ", transparency[0])

for ue4 in [0,23]:#range(25):
    for um4 in [0,23]:#range(25):
        fig, ax = plt.subplots(1, 1)
        for f in range(len(files)-1):
            plt.scatter(np.arange(800),np.array(chi_raw[f][:,ue4,um4]),alpha=transparency[f][:],marker='.',color=data_colors[f],linestyle='None',label=labels[f])
            #massindex=np.argmax(np.array(chi_raw[f][:,ue4,um4]))
            #print("max index, value in "+labels[f]+" = "+str(massindex)+", "+str(chi_raw[f][massindex,ue4,um4]))
        legend = ax.legend(loc='upper left', shadow=True)
        ax.set_xlabel(r'$\Delta m^2 points$')
        ax.set_ylabel(r'Raw $ \chi^2/31 $')
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.set_title(r'$U_{e4}$ = '+ str(ue4_list[ue4]) +', $U_{\mu4}$ = '+ str(um4_list[um4]))
        fig.savefig('histo_MFAvalid_chi_raw_ue4_'+str(ue4_list[ue4])+'_um4_'+ str(um4_list[um4]) + '.pdf',transparent=True)
        #plt.ylim(-0.1, 10)
        plt.yscale('log')
        fig.savefig('histo_logscale_MFAvalid_chi_raw_ue4_'+str(ue4_list[ue4])+'_um4_'+ str(um4_list[um4]) + '.pdf',transparent=True)
        plt.close()

bin_edges = np.arange(31)   
#for m in [5, 655, 705, 755, 795]:
for ue4 in [0,23]:#range(25):
    for um4 in [23]:#range(25):
        for m in [21, 617, 692, 782, 797]:
            axes = []
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, sharey=True)
            axes.append(ax1)
            axes.append(ax2)
            axes.append(ax3)
            axes.append(ax4)
            fig.suptitle(r'$\Delta m^2$ = '+str( m41_list[m] )+', $U_{e4}$ = ' +str( ue4_list[ue4] )+', $U_{\mu4}$ = ' +str( ue4_list[um4] ) )
            for f in range(len(files)-1):
                axes[f].bar(x=bin_edges[:], height=speccoll_diff[f][m,ue4,um4,:],  alpha=0.99, align='edge', fc=data_colors[f], ec='white' )
                axes[f].set_title( labels[f], fontsize = 9 )
                axes[f].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
                axes[f].tick_params(axis='both', which='major', labelsize=7)
                #axes[f-1].bar(x=bin_edges[:], height=speccoll_mfa[f][m,ue4,um4,:],  alpha=0.65, align='edge', fc=data_colors[f+len(files)], ec='white', label=labels[f]+', MFA' )
            fig.text(0.5, 0.04, 'Bins', fontname = 'Avenir', fontsize = 12, ha='center')
            fig.text(0.03, 0.5, r'$MFA (\vec{q}) - \cal{S} (\vec{q})$', fontname = 'Avenir', fontsize = 12, va='center', rotation='vertical')
            fig.savefig('histo_MFAvalid_specGSvsMFA_m4='+str(m41_list[m])+'_ue4_'+str(ue4_list[ue4])+'_um4_'+ str(um4_list[um4]) + '.pdf',transparent=True)
            plt.close()
            axes = []
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, sharey=True)
            axes.append(ax1)
            axes.append(ax2)
            axes.append(ax3)
            axes.append(ax4)
            fig.suptitle(r'$\Delta m^2$ = '+str( m41_list[m] )+', $U_{e4}$ = ' +str( ue4_list[ue4] )+', $U_{\mu4}$ = ' +str( ue4_list[um4] ) )
            for f in range(len(files)-1):
                axes[f].bar(x=bin_edges[:], height=speccoll_mfa[f][m,ue4,um4,:],  alpha=0.99, align='edge', fc='none', ec=data_colors[f], label='MFA Spectrum' )
                axes[f].bar(x=bin_edges[:], height=speccoll_gs[f][m,ue4,um4,:],  alpha=0.99, align='edge', fc='none', ec=data_colors[f+len(files)], linestyle='dashed',label='Input Spectrum' )
                axes[f].set_title( labels[f], fontsize = 9 )
                axes[f].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
                axes[f].tick_params(axis='both', which='major', labelsize=7)
                legend = axes[f].legend(loc='upper left', shadow=True)
                axes[f].legend(fontsize=5)
		#calculate chi2
                chir=np.divide(np.sum(np.divide(np.multiply((speccoll_diff[f][m,ue4,um4,:]),np.sqrt(speccoll_gs[f][m,ue4,um4,:])),speccoll_diff[f][m,ue4,um4,:])),31)
                print(str(labels[f])+", chi2="+str('{:.2e}'.format(chir))+", m="+str(m41_list[m])+", ue4="+str(ue4_list[ue4])+" um4="+str(um4_list[um4]))
                print("array=", np.multiply(np.multiply((speccoll_diff[f][m,ue4,um4,:]),np.sqrt(speccoll_gs[f][m,ue4,um4,:])),speccoll_diff[f][m,ue4,um4,:]))
                print("stat err =", np.sqrt(speccoll_gs[f][m,ue4,um4,:]))
                print("diff =", speccoll_diff[f][m,ue4,um4,:])
                #axes[f-1].bar(x=bin_edges[:], height=speccoll_mfa[f][m,ue4,um4,:],  alpha=0.65, align='edge', fc=data_colors[f+len(files)], ec='white', label=labels[f]+', MFA' )
            fig.text(0.5, 0.04, 'Bins', fontname = 'Avenir', fontsize = 12, ha='center')
            fig.text(0.03, 0.5, r'$MFA (\vec{q}) - \cal{S} (\vec{q})$', fontname = 'Avenir', fontsize = 12, va='center', rotation='vertical')
            fig.savefig('histo_MFAvalid_specGSvsMFAoverlay_m4='+str(m41_list[m])+'_ue4_'+str(ue4_list[ue4])+'_um4_'+ str(um4_list[um4]) + '.pdf',transparent=True)
            plt.close()
