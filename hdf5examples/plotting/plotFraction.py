import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np

# Define the exponential function to fit
#def func(x, a, b, c):
#    return a * np.exp(-b * x) + c
# Define the quadratic function to fit the data
def func(x, a, b, c):
    return a*x**2 + b*x + c


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



# Sample data
#cases = ['25', '50', '100', '200', '400', '800', '1600']
#cases = ['100', '200', '400', '800', '1600']
cases = ['100', '200', '400', '800']
threshold = 0 #0.0*0.005
data = np.random.rand(len(cases), 100)

# Fraction of data points above threshold
# load in files
#f1 = h5.File("test_scan.h5", mode="r")
#f1 = h5.File("testnewminimizer_"+str(cases[0])+"x25x25_iter50_10000univs2.h5", mode="r")
#f2 = h5.File("testnewminimizer_"+str(cases[1])+"x25x25_iter50_10000univs2.h5", mode="r")
#f3 = h5.File("testnewmin_"+str(cases[2])+"x25x25_gp0_iter50_10000univ_3.h5", mode="r")
#f4 = h5.File("testnewminimizer_"+str(cases[3])+"x25x25_iter50_10000univs2.h5", mode="r")
#f5 = h5.File("testnewminimizer_"+str(cases[4])+"x25x25_iter50_10000univs2.h5", mode="r")
#
##f6 = h5.File("test_scan2.h5", mode="r")
#f6 = h5.File("highres1/testnewmin_"+str(cases[0])+"x50x50_gpt0_iter50_10000univs.h5", mode="r")
#f7 = h5.File("highres1/testnewmin_"+str(cases[1])+"x50x50_gpt0_iter50_10000univs.h5", mode="r")
#f8 = h5.File("highres1/testnewmin_"+str(cases[2])+"x50x50_gpt0_iter50_10000univs.h5", mode="r")
#f9 = h5.File("highres1/testnewmin_"+str(cases[3])+"x50x50_gpt0_iter50_10000univs.h5", mode="r")
#f10 = h5.File("highres1/testnewmin_"+str(cases[4])+"x50x50_gpt0_iter50_10000univs.h5", mode="r")

f1 = h5.File("MFAFiles-perlmutter/testminimizernewgrad_"+str(cases[0])+"x25x25_0to15625_iter1.h5", mode="r")
f2 = h5.File("MFAFiles-perlmutter/testminimizernewgrad_"+str(cases[1])+"x25x25_0to15625_iter1.h5", mode="r")
f3 = h5.File("MFAFiles-perlmutter/testminimizernewgrad_"+str(cases[2])+"x25x25_0to15625_iter1.h5", mode="r")
f4 = h5.File("MFAFiles-perlmutter/testminimizernewgrad_"+str(cases[3])+"x25x25_0to15625_iter1.h5", mode="r")
#f5 = h5.File("MFAFiles-perlmutter/testminimizernewgrad_"+str(cases[3])+"x25x25_0to15625_iter1.h5", mode="r")
#f13 = h5.File("MFAFiles-perlmutter/testminimizernewgrad_"+str(cases[4])+"x25x25_0to15625_iter1.h5", mode="r")
#f14 = h5.File("MFAFiles-perlmutter/testminimizernewgrad_"+str(cases[4])+"x25x25_0to15625_iter1.h5", mode="r")
#
f6  = h5.File("testmin_100x25x25_0to15625.h5", mode="r")
f7  = h5.File("testmin_200x25x25_0to15625.h5", mode="r")
f8  = h5.File("testmin_400x25x25_0to15625.h5", mode="r")
f9  = h5.File("testmin_800x25x25_0to15625.h5", mode="r")
#f10 = h5.File("testmin_400x25x25_0to15625.h5", mode="r")
#f6  = h5.File("testnewmin_"+str(cases[0])+"x"+str(cases[0])+"_0to242500_iter1_1.h5", mode="r")
#f7  = h5.File("testnewmin_"+str(cases[1])+"x"+str(cases[1])+"_0to242500_iter1_1.h5", mode="r")
#f8  = h5.File("testnewmin_"+str(cases[2])+"x"+str(cases[2])+"_0to242500_iter1_1.h5", mode="r")
#f9  = h5.File("testnewmin_"+str(cases[4])+"x"+str(cases[4])+"_0to242500_iter1_1.h5", mode="r")
#f10 = h5.File("testnewmin_"+str(cases[4])+"x"+str(cases[4])+"_0to242500_iter1_1.h5", mode="r")
#f11 = h5.File("testnewmin_"+str(cases[4])+"x"+str(cases[4])+"_0to242500_iter1_1.h5", mode="r")
#f12 = h5.File("testnewmin_"+str(cases[4])+"x"+str(cases[4])+"_0to242500_iter1_1.h5", mode="r")
#f11 = h5.File("MFAFiles-perlmutter/testminimizernewgrad_"+str(cases[4])+"x50x50_0to250000_iter1.h5", mode="r")
#f12 = h5.File("MFAFiles-perlmutter/testminimizernewgrad_"+str(cases[4])+"x50x50_0to250000_iter1.h5", mode="r")

#f6 = h5.File("MFAFiles-perlmutter/testminimizernewgrad_"+str(cases[0])+"x50x50_0to250000_iter1.h5", mode="r")
#f7 = h5.File("MFAFiles-perlmutter/testminimizernewgrad_"+str(cases[1])+"x50x50_0to250000_iter1.h5", mode="r")
#f8 = h5.File("MFAFiles-perlmutter/testminimizernewgrad_"+str(cases[2])+"x50x50_0to250000_iter1.h5", mode="r")
#f9 = h5.File("MFAFiles-perlmutter/testminimizernewgrad_"+str(cases[3])+"x50x50_0to250000_iter1.h5", mode="r")
#f10 = h5.File("MFAFiles-perlmutter/testminimizernewgrad_"+str(cases[4])+"x50x50_0to250000_iter1.h5", mode="r")
#f11 = h5.File("MFAFiles-perlmutter/testminimizernewgrad_"+str(cases[5])+"x50x50_0to250000_iter1.h5", mode="r")
#f12 = h5.File("MFAFiles-perlmutter/testminimizernewgrad_"+str(cases[6])+"x50x50_0to250000_iter1.h5", mode="r")
#
#f1 = h5.File("testminimizernewgrad_gp7499_"+str(cases[0])+"x25x25_0to15625_iter1.h5", mode="r")
#f2 = h5.File("testminimizernewgrad_gp7499_"+str(cases[1])+"x25x25_0to15625_iter1.h5", mode="r")
#f3 = h5.File("testminimizernewgrad_gp7499_"+str(cases[2])+"x25x25_0to15625_iter1.h5", mode="r")
#f4 = h5.File("testminimizernewgrad_gp7499_"+str(cases[3])+"x25x25_0to15625_iter1.h5", mode="r")
#f5 = h5.File("testminimizernewgrad_gp7499_"+str(cases[4])+"x25x25_0to15625_iter1.h5", mode="r")
#
#f6 = h5.File("testminimizernewgrad_gp152349_"+str(cases[0])+"x50x50_0to250000_iter1.h5", mode="r")
#f7 = h5.File("testminimizernewgrad_gp152349_"+str(cases[1])+"x50x50_0to250000_iter1.h5", mode="r")
#f8 = h5.File("testminimizernewgrad_gp152349_"+str(cases[2])+"x50x50_0to250000_iter1.h5", mode="r")
#f9 = h5.File("testminimizernewgrad_gp152349_"+str(cases[3])+"x50x50_0to250000_iter1.h5", mode="r")
#f10 = h5.File("testminimizernewgrad_gp152349_"+str(cases[4])+"x50x50_0to250000_iter1.h5", mode="r")
#
#f11 = h5.File("testminimizernewgrad_gp7499_"+str(cases[0])+"x25x25_0to15625_iter1.h5", mode="r")
#f12 = h5.File("testminimizernewgrad_gp7499_"+str(cases[1])+"x25x25_0to15625_iter1.h5", mode="r")
#f13 = h5.File("testminimizernewgrad_gp7499_"+str(cases[2])+"x25x25_0to15625_iter1.h5", mode="r")
#f14 = h5.File("testminimizernewgrad_gp7499_"+str(cases[3])+"x25x25_0to15625_iter1.h5", mode="r")
#f15 = h5.File("testminimizernewgrad_gp7499_"+str(cases[4])+"x25x25_0to15625_iter1.h5", mode="r")

#load in chi2
#min1   = f1['last_chi_min'][:]
#delta1   = f1['delta_chi'][:]
#data1   = np.add(min1,delta1)
data1   = f1['last_chi_min'][:]
data2   = f2['last_chi_min'][:]
data3   = f3['last_chi_min'][:]
data4   = f4['last_chi_min'][:]
#data5   = f5['last_chi_min'][:]
#data13   = f5['last_chi_min'][:]
#data14   = f5['last_chi_min'][:]
starty1   = f1['best_start_pointy'][:]
starty2   = f2['best_start_pointy'][:]
starty3   = f3['best_start_pointy'][:]
starty4   = f4['best_start_pointy'][:]
#starty5   = f5['best_start_pointy'][:]
startz1   = f1['best_start_pointz'][:]
startz2   = f2['best_start_pointz'][:]
startz3   = f3['best_start_pointz'][:]
startz4   = f4['best_start_pointz'][:]
#startz5   = f5['best_start_pointz'][:]

#min6   = f6['last_chi_min'][:]
#delta6   = f6['delta_chi'][:]
#data6   = np.add(min6,delta6)
data6    = f6['last_chi_min'][:]
data7    = f7['last_chi_min'][:]
data8    = f8['last_chi_min'][:]
data9    = f9['last_chi_min'][:]
#data10   = f10['last_chi_min'][:]
starty6   = f6['best_start_pointy'][:]
starty7   = f7['best_start_pointy'][:]
starty8   = f8['best_start_pointy'][:]
starty9   = f9['best_start_pointy'][:]
#starty10   = f10['best_start_pointy'][:]
startz6   = f6['best_start_pointz'][:]
startz7   = f7['best_start_pointz'][:]
startz8   = f8['best_start_pointz'][:]
startz9   = f9['best_start_pointz'][:]
#startz10   = f10['best_start_pointz'][:]

#data11    = f10['last_chi_min'][:]
#data12    = f10['last_chi_min'][:]
#data13    = f13['last_chi_min'][:]
#data14    = f14['last_chi_min'][:]
#data15    = f15['last_chi_min'][:]

print("data1.shape[0]: ", data1.shape[0])
fraction1 = np.sum((data1 >= 0.0) & (data1 <= 129.2), axis=0) / data1.shape[0]
fraction2 = np.sum((data2 >= 0.0) & (data2 <= 129.2), axis=0) / data2.shape[0]
fraction3 = np.sum((data3 >= 0.0) & (data3 <= 129.2), axis=0) / data3.shape[0]
fraction4 = np.sum((data4 >= 0.0) & (data4 <= 129.2), axis=0) / data4.shape[0]
#fraction5 = np.sum((data5 >= 0.0) & (data5 <= 129.2), axis=0) / data5.shape[0]
#fraction13 = np.sum((data13 >= 0.0) & (data13 <= 129.2), axis=0) / data13.shape[0]
#fraction14 = np.sum((data14 >= 0.0) & (data14 <= 129.2), axis=0) / data14.shape[0]


fraction6 = np.sum((data6 >= 0.0) & (data6 <= 129.2), axis=0) / data6.shape[0]
fraction7 = np.sum((data7 >= 0.0) & (data7 <= 129.2), axis=0) / data7.shape[0]
fraction8 = np.sum((data8 >= 0.0) & (data8 <= 129.2), axis=0) / data8.shape[0]
fraction9 = np.sum((data9 >= 0.0) & (data9 <= 129.2), axis=0) / data9.shape[0]
#fraction10 = np.sum((data10 >= 0.0) & (data10 <= 129.2), axis=0) / data10.shape[0]
#fraction6 = np.sum((data6 >= 0.0) & (data6 <= 129.1) & (starty6 <= 48.) & (startz6 <= 48.), axis=0) / np.sum((starty6 <= 48.) & (startz6 <= 48.), axis=0) #data6.shape[0]
#fraction7 = np.sum((data7 >= 0.0) & (data7 <= 129.1) & (starty7 <= 48.) & (startz7 <= 48.), axis=0) / np.sum((starty7 <= 48.) & (startz7 <= 48.), axis=0) #data7.shape[0]
#fraction8 = np.sum((data8 >= 0.0) & (data8 <= 129.1) & (starty8 <= 48.) & (startz8 <= 48.), axis=0) / np.sum((starty8 <= 48.) & (startz8 <= 48.), axis=0) #data8.shape[0]
#fraction9 = np.sum((data9 >= 0.0) & (data9 <= 129.1) & (starty9 <= 48.) & (startz9 <= 48.), axis=0) / np.sum((starty9 <= 48.) & (startz9 <= 48.), axis=0) #data9.shape[0]
#fraction10 = np.sum((data10 >= 0.0) & (data10 <= 129.1) & (starty10 <= 48.) & (startz10 <= 48.), axis=0) / np.sum((starty10 <= 48.) & (startz10 <= 48.), axis=0) #data10.shape[0]
#fraction11 = np.sum((data11 >= 0.0) & (data11 <= 129.1), axis=0) / data11.shape[0]
#fraction12 = np.sum((data12 >= 0.0) & (data12 <= 129.1), axis=0) / data12.shape[0]

#print("data shape: ", data1.shape[0], data2.shape[0], data3.shape[0], data4.shape[0], data5.shape[0], data6.shape[0], data7.shape[0], data8.shape[0], data9.shape[0], data10.shape[0], data11.shape[0], data12.shape[0])
#print("totpoints excl edges of nue and numu: ", np.sum((starty6 <= 48.) & (startz6 <= 48.), axis=0), np.sum((starty7 <= 48.) & (startz7 <= 48.), axis=0), np.sum((starty8 <= 48.) & (startz8 <= 48.), axis=0), np.sum((starty9 <= 48.) & (startz9 <= 48.), axis=0), np.sum((starty10 <= 48.) & (startz10 <= 48.), axis=0))

#fraction15 = np.sum(data15 <= (130.6), axis=0) / data15.shape[0]
#fig, ax = plt.subplots(1, 1)
#plt.hist(data5,bins=range(120, 140 + 1, 1))
#fig.savefig('histo.pdf',transparent=True)
#plt.close()

#print("fraction1, fraction2, fraction3, fraction4, fraction5: ", fraction1, fraction2, fraction3, fraction4, fraction5)
#print("fraction6, fraction7, fraction8, fraction9, fraction10: ", fraction6, fraction7, fraction8, fraction9, fraction10)
#array of fractions

fractions1 = [100.*fraction1[0], 100.*fraction2[0], 100.*fraction3[0], 100.*fraction4[0]]#, 100.*fraction5[0]]#, 100.*fraction5[0], 100.*fraction5[0]]
fractions2 = [100.*fraction6[0], 100.*fraction7[0], 100.*fraction8[0], 100.*fraction9[0]]#, 100.*fraction10[0]]#, 100.*fraction11[0], 100.*fraction12[0]]
#fractions1 = [100., 100.*fraction1[0], 100.*fraction2[0], 100.*fraction3[0], 100.*fraction4[0], 100.*fraction5[0]]#, 100.*fraction5[0], 100.*fraction5[0]]
#fractions2 = [100., 100.*fraction6[0], 100.*fraction7[0], 100.*fraction8[0], 100.*fraction9[0], 100.*fraction10[0]]#, 100.*fraction11[0], 100.*fraction12[0]]
#fractions3 = [0., 100.*fraction11[0], 100.*fraction12[0], 100.*fraction13[0], 100.*fraction14[0], 100.*fraction15[0], 0., 0.]
print("fractions 1: ", fractions1)
print("fractions 2: ", fractions2)

# Set the y-axis for the time required to build the MFA signal model 
fifty=[0., 1.18284, 2.37443, 6.57478, 21.1398, 116.436, 894.052, 3410.22]
twentyfive=[0., 0.248391, 0.537951, 1.32992, 3.72901, 25.3193, 216.914, 850.416]

# Set the y-axis for the time required to run FC 
FCb=[0., 143.549, 145.787, 166.47, 195.04, 365.251, 971.541, 1760.25, ]
FCa=[0., 0.248391, 0.537951, 1.32992, 3.72901, 25.3193, 216.914, 850.416]

# Define x-axis data
#x_data = np.array([0, 1, 2, 4, 8, 16, 32, 64])
x_data = np.array([4, 8, 16, 32])#, 64])

# Map the x_data array to the desired x_ticks array
#x_ticks = np.array([ '0', '', '', '100', '', '400', '800', '1600'])
x_ticks = np.array([ '100', '', '400', '800'])#, '1600'])

#------------------------------------------------------------------

# Line plot
fig, ax = plt.subplots(1, 1)
# Set the padding for the bottom of the plot
plt.subplots_adjust(bottom=0.15)
plt.plot(x_data, fractions1, data_style, color=data_colors[0], label='grid resolution 25x25x25 (with bug)')
plt.plot(x_data, fractions1, fit_style, color=data_colors[0], label='_nolegend_')
#plt.plot(x_data, fractions2, data_style, color=data_colors[1], label='grid resolution 100x50x50')
plt.plot(x_data, fractions2, data_style, color=data_colors[1], label='grid resolution 25x25x25 (bug fixed)')
plt.plot(x_data, fractions2, fit_style, color=data_colors[1], label='_nolegend_' )

# Set plot title and labels
plt.title(r'Fraction of Randomized Starting Point within the $R_{null}$ Best-Fit point', fontsize=12, fontweight='bold')
plt.xlabel(r'Number of $\Delta m^{2}$ input points for MFA model')
plt.ylabel(r'Fraction (%)')
plt.xticks(x_data, x_ticks)
plt.legend()

# Show plot
fig.savefig('fraction.pdf',transparent=True)
plt.close()

# Line plot
fig, ax = plt.subplots(1, 1)

## Set the padding for the bottom of the plot
#plt.subplots_adjust(bottom=0.15)
#
## Perform curve fitting
#popt, pcov = curve_fit(func, x_data, fifty)
#popt2, pcov2 = curve_fit(func, x_data, twentyfive)
#
## Print the optimal parameters
#print(popt)
#
## Set custom colors and linestyles
#
## Plot the data and the fit
##x = np.linspace(1, 64, 100)
#plt.plot(x_data, func(x_data, *popt2), fit_style, color=fit_color, label=f'fit: {popt2[0]:.2f}$x^2$ + {popt2[1]:.2f}$x$ + {popt2[2]:.2f}')
#plt.plot(x_data, twentyfive, data_style, color=data_colors[0], label='25x25')
#plt.plot(x_data, func(x_data, *popt), '--', color=fit_color, label=f'fit: {popt[0]:.2f}$x^2$ + {popt[1]:.2f}$x$ + {popt[2]:.2f}')
#plt.plot(x_data, fifty, data_style, color=data_colors[1], label='50x50')
#
## Set custom title font size and weight
#plt.title('Building MFA Signal (using 4 ranks)', fontsize=16, fontweight='bold')
#
## Set plot title and labels
#plt.xlabel(r'Number of $\Delta m^{2}$ input points for MFA model')
#plt.ylabel(r'Time (sec)')
#plt.xticks(x_data, x_ticks)
#plt.tick_params(axis='x', which='major', labelsize=12)
#
#plt.legend()
#
## Show plot
#fig.savefig('SignalModel.pdf',transparent=True)
