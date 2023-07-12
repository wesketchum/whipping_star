#import required libraries
import sys, getopt
import h5py as h5
import numpy as np
import math as math
import matplotlib.pyplot as plt 
from matplotlib.colors import hsv_to_rgb
from cycler import cycler
import matplotlib.cm as cm
import seaborn as sns 
from matplotlib.colors import LogNorm
from sympy import pretty_print as pp, latex
from sympy.abc import a, b, n
import scipy.interpolate
from scipy.interpolate import interp2d


def makeHistograms(spec, bin_edges, dmsq, ue4, um4, ftname, out ):
        fig, axes = plt.subplots(6, 3, sharex='col', sharey='row', figsize=(30,20))
        for i in range(6):
                for j in range(3):
                        k = i*3+j;
                        if k >= len(dmsq):
                                continue;
                        mass = dmsq[k]
                        Ue4 = ue4[k]
                        Umu4 = um4[k]
                        #print(i, j, k, bin_edges[:], spec[k])

                        axes[i,j].bar(x=bin_edges[:], height=spec[k], alpha=0.99, align='edge', fc='#1f77b4', ec='white')
                        axes[i,j].set_title( "Bin $\Delta m^2$ = "+str(mass)+", Ue4 = "+str(Ue4)+", Umu4 ="+str(Umu4), fontsize = 20)
        
        fig.text(0.5, 0.04, 'Bins', fontname = ftname, fontsize = 40, ha='center')
        fig.text(0.04, 0.5, r'$MFA (\mathbf{q}) - \cal{S} (\mathbf{q})$', fontname = ftname, fontsize = 40, va='center', rotation='vertical')
        fig.savefig(out+'.pdf',transparent=True)


# Read H5 file
#f0 = h5.File("test_m122.h5", mode="r") #/Users/wospakrk/test_.h5
#f0 = h5.File("/Users/wospakrk/testnewmin_coarsermfa25x25_1universes_0to2500_scan.h5", mode="r") #/Users/wospakrk/test_.h5
f0 = h5.File("/Users/wospakrk/testnewmin_coarsermfa50x50x50_1universes_0to2500_scan.h5", mode="r") #/Users/wospakrk/test_.h5
nuniv=1000.
#: extract reflectance data from the H5 file
spec_mfa = f0['speccoll_mfa']
spec_gp  = f0['speccoll_gs']
thischi_mfa = f0['this_chi_mfa']
thischi_gp  = f0['this_chi_gs']

#calculte the relative errors:
#relerr = np.true_divide(np.subtract(thischi_gp,thischi_mfa),thischi_gp) 
#
##print the indices
#print("relerr: ", relerr )
#largediff=[i for i,v in enumerate(relerr) if v < -0.0005]
#
#diffs=[]
bin_edges = np.arange(31)
#for i in largediff:
#    diffs.append(np.true_divide(np.subtract(spec_gp[i:],spec_mfa[i:]), spec_gp[i:]))
#print("diffs: ", len(diffs[0][0]), len(bin_edges))

#load in grid parameters
xgrid   = f0['gridx'][:]
ygrid   = f0['gridy'][:]
zgrid   = f0['gridz'][:]
totpoints = xgrid.shape[0]

#load in chi2 values
speccoll_gs  = f0['speccoll_gs'][:]
speccoll_mfa = f0['speccoll_mfa'][:]
diff_spec = np.subtract(speccoll_mfa, speccoll_gs);

#prepare grid x-y
x = 10**(xgrid*2.)
x.reshape(x.shape[0],)
dx=np.sort(np.unique(x))
y = 10**(ygrid)
y.reshape(y.shape[0],)
dy=np.sort(np.unique(y))
z = 10**(zgrid)
z.reshape(z.shape[0],)
dz=np.sort(np.unique(z))

print("len dx, dx: ", len(dx), dx)
print("len dy, dy: ", len(dy), dy)
print("len dz, dz: ", len(dz), dz)

indexx=np.arange(len(dx))
indexy=np.arange(len(dy))
indexz=np.arange(len(dz))

#reshape 1d vector to 3d
speccoll_gs_3d  = np.reshape(speccoll_gs, (len(dx), len(dy), len(dz), 31))
speccoll_mfa_3d = np.reshape(speccoll_mfa, (len(dx), len(dy), len(dz), 31))
diff_spec_4d = np.reshape(diff_spec, (len(dx), len(dy), len(dz), 31))

print("gs: ", speccoll_gs_3d[0,47,47,:])
print("mfa: ", speccoll_mfa_3d[0,47,47,:])
print("diff: ", diff_spec_4d[0,47,47,:])

vec_diff_specs = [diff_spec_4d[1,47,47,:], diff_spec_4d[1,43,49,:], diff_spec_4d[1,45,49,:], diff_spec_4d[1,47,49,:], diff_spec_4d[1,49,47,:], diff_spec_4d[1,49,49,:], diff_spec_4d[51,47,47,:], diff_spec_4d[51,43,49,:], diff_spec_4d[51,45,49,:], diff_spec_4d[51,47,49,:], diff_spec_4d[51,49,47,:], diff_spec_4d[51,49,49,:], diff_spec_4d[97,47,47,:], diff_spec_4d[97,43,49,:], diff_spec_4d[97,45,49,:], diff_spec_4d[97,47,49,:], diff_spec_4d[97,49,47,:], diff_spec_4d[97,49,49,:]]

dmsq=[1,1,1,1,1,1,51,51,51,51,51,51,97,97,97,97,97,97]
ue4=[47,43,45,47,49,49,47,43,45,47,49,49,47,43,45,47,49,49]
um4=[47,49,49,49,47,49,47,49,49,49,47,49,47,49,49,49,47,49]

makeHistograms( vec_diff_specs, bin_edges, dmsq, ue4, um4, 'Avenir', 'dmsq_odd_spectracomparisons')
