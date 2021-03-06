# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 10:35:13 2020

Last Edited: 26/06/2020

@author: sophia

Easier way to plot data from SpinRatesAreHere script. 
Allows input of cluster name and produces the following outputs:
    -Period-Mass plot(linear + log space)
    -RA-Dec map
    -Mass histogram (normal + cumulative)
    -Period histogram (linear + log space)
    -Mass model comparison (requires editing of SpinRatesAreHere.py) [18/06/2020]
    -Calculation of cluster width (physicial and angular) [18/06/2020]
    -All plots separate disked and diskless stars [23/06/2020]
    -KS test on disked and diskless, as well as low and high mass samples within a cluster [23/06/2020]
    -Plot of disk fraction as a function of period [24/06/2020]
    -Rolling quantiles at 10%, 50%, and 90% for period-mass [25/06/2020]
    -Rolling disk fraction for disk fraction vs period plot [26/06/2020]
    
Note [23/06/2020]: with updates to include separation between disked and non disked stars, hPer can no
longer be used without commenting out/changing variables that use cluster.Disk. Easiest to
just remove parameters using cluster.Disk in the DISK section rather than editing any of
the functions later on.
e.g. comment out disked=np.where....

Edit [25/06/2020]: has_mass_prot has been added for ease, and PeriodMassPlot() has a line that can be used
to plot hPer, however the rest of the functions do not have accomodations for hPer that can
just be uncommented yet.
    
"""
"""
============================	PACKAGES	=================================
"""
import matplotlib.pyplot as plt
import numpy as np

import hciplot
import sys
import os

import math
import pandas as pd
from math import cos
from math import sin
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pylab import figure
import warnings
warnings.simplefilter("ignore")

import SpinRatesAreHere as spin 


"""
End of package imports 
"""

#asks user to input name of cluster they wish to use
name_input = 'NGC2362'
#name_input = input('Name of cluster:  ')

if name_input =='NGC2264':
    cluster = spin.NGC2264()
    
elif name_input == 'hPer':
    cluster = spin.hPer()
    
elif name_input == 'USco':
    cluster = spin.USco()
    
elif name_input == 'NGC2362':
    cluster = spin.NGC2362()
    
else:
    print('Invalid name')
    exit()
    
    
"""
============================	DISK 	=================================
"""
#DO NOT USE FOR HPER AS IT IS TOO OLD TO HAVE DISK INFORMATION

#want to separate stars with and without disks in the cluster so that they can be compared
#assigning variables for disked and diskless stars
disked=np.where(cluster.Disk>0)[0]
diskless=np.where(cluster.Disk<0)[0]
unknown=np.where(cluster.Disk==0)[0]


#only want stars that have mass, disk, and period data
mass=np.where(np.isfinite(cluster.Mass))[0]
prot=np.where(np.isfinite(cluster.Prot))[0]
disk=np.where(np.isfinite(cluster.Disk))[0]

#separating samples with masses above and below 0.4 solar masses
high_mass=np.where(cluster.Mass>0.4)[0]
low_mass=np.where(cluster.Mass<=0.4)[0]

#has_mass_prot=np.array(list(set(mass)&set(prot))) #for hPer
has_disk_mass_prot=np.array(list(set(mass)&set(disk)&set(prot)))
has_unknown_mass_prot=np.array(list(set(mass)&set(unknown)&set(prot)))
has_disked_mass_prot=np.array(list(set(mass)&set(disked)&set(prot)))
has_diskless_mass_prot=np.array(list(set(mass)&set(diskless)&set(prot)))

has_disk_high_mass_prot=np.array(list(set(high_mass)&set(disk)&set(prot)))
has_disk_low_mass_prot=np.array(list(set(low_mass)&set(disk)&set(prot)))

#calculation of disk fraction for eligable stars in the cluster
#print('Stars with mass, period, and disk information= ',len(has_disk_mass_prot))
#print('Stars with disks= ',len(has_disked_mass_prot))
#print('Disk fraction(%)= ',(len(has_disked_mass_prot)/len(has_disk_mass_prot))*100)


"""
============================	PLOTS	=================================
"""

#plots scatter plot of period vs mass
def PeriodMassPlot():
    plt.scatter(cluster.Mass[disked], cluster.Prot[disked], color='r', s=4, label='Disk')
    plt.scatter(cluster.Mass[diskless], cluster.Prot[diskless], color='b', s=4, label='Diskless')
    plt.legend(loc='upper right')
    plt.ylabel('Period [days]')
    plt.xlabel('Mass [$M_{sun}$]')
    plt.title('Period-Mass plot of {name}'.format(name=name_input))
    plt.xlim(xmin=0)
    plt.savefig('{name}_disk_period_mass.png'.format(name=name_input))
    plt.show() 

#same as above but with period on a log scale
def PeriodMassLogPlot():
    plt.scatter(cluster.Mass[disked], cluster.Prot[disked], alpha=0.5, color='r', s=4, label='Disk')
    plt.scatter(cluster.Mass[diskless], cluster.Prot[diskless], alpha=0.5, color='b', s=4, label='Diskless')
    #plt.scatter(cluster.Mass[unknown], cluster.Prot[unknown], alpha=0.5, color='y', s=4, label='Unknown')
    #plt.scatter(cluster.Mass, cluster.Prot, alpha=0.5, color='b', s=4)
    
    RollingQuant()
    #plt.legend(loc='upper right')
    plt.yscale('log')
    plt.ylabel('Period [days]')
    plt.xlabel('Mass [$M_{sun}$]')
    plt.title('Period-Mass plot of {name} (period log scale)'.format(name=name_input))
    #plt.xlim(xmin=0)
    plt.savefig('{name}_rolling_quantile.png'.format(name=name_input))
    plt.show() 
       
#Compares the different mass models available in the database. To use, go into
#SpinRatesAreHere.py and comment out the if/else statements and rename self.Mass
#to self.Mass_0,1,2 respectively for the mass_type part of the function  
def MassCompare():
    plt.scatter(cluster.Mass_0, cluster.Prot, s=4, alpha=0.6, color='c', label='Irwin et al. 2008')
    plt.scatter(cluster.Mass_1, cluster.Prot, s=4, alpha=0.6, color='b', label='MESA model')
    plt.scatter(cluster.Mass_2, cluster.Prot, s=4, alpha=0.6, color='m', label='Baraffe 98 model')
    plt.legend(loc='lower right')
    plt.yscale('log')
    plt.ylabel('Period [days]')
    plt.xlabel('Mass [$M_{sun}$]')
    plt.title('Period-Mass plot of {name}'.format(name=name_input))
    plt.xlim(xmin=0)
    plt.savefig('{name}_mass_comparison_log.png'.format(name=name_input))
    plt.show()

def MassHistCompare():
    bin_list=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4]
    plt.hist(cluster.Mass_0, bins=bin_list, rwidth=0.9, alpha=0.6, edgecolor='black', color='c', cumulative=False, label='Irwin et al. 2008')
    plt.hist(cluster.Mass_1, bins=bin_list, rwidth=0.9, alpha=0.6, edgecolor='black', color='b', cumulative=False, label='MESA model')
    #plt.hist(cluster.Mass_2, bins=bin_list, rwidth=0.9, alpha=0.6, edgecolor='black', color='m', cumulative=False, label='Baraffe 98 model')
    plt.ylabel('Number of stars')
    plt.xlabel('Mass [$M_{sun}$]')
    plt.legend(loc='upper right')
    plt.xlim(xmin=0)
    plt.title('Mass histogram for {name}'.format(name=name_input))
    #plt.savefig('{name}_mass_hist_compare.png'.format(name=name_input))
    plt.show()
    
#plots RA vs Dec map 
def CoordMap():
    if name_input == 'NGC2362': #or name_input == 'USco' or
        import astropy.coordinates as ac
        cat_disked = ac.SkyCoord(cluster.RA[disked], cluster.Dec[disked], unit="deg")
        cat_diskless = ac.SkyCoord(cluster.RA[diskless], cluster.Dec[diskless], unit="deg")
        plt.scatter(cat_disked.ra.deg*15, cat_disked.dec.deg, color='r', s=4, label='Disk') #RA is *15 to convert from time to degrees
        plt.scatter(cat_diskless.ra.deg*15, cat_diskless.dec.deg, color='b', s=4, label='Diskless') 
    else:
        plt.scatter(cluster.RA[disked], cluster.Dec[disked], color='r', s=4, label='Disk')
        plt.scatter(cluster.RA[diskless], cluster.Dec[diskless], color='b', s=4, label='Diskless')
    plt.legend(loc='upper right')    
    plt.xlabel('Right Ascension [deg]')
    plt.ylabel('Declination [deg]')
    plt.title('RA-Dec Map of {name}'.format(name=name_input))
    plt.savefig('{name}_disk_ra_dec_map.png'.format(name=name_input))
    plt.show() 
    

#plots normal histogram of stellar masses
def MassHistNorm():
    #bin_list=
    fig, axes = plt.subplots(1,2, figsize=(8,4))
    
    plt.subplot(axes[0])
    axes[0].hist(cluster.Mass[disked], bins=20, rwidth=0.9, cumulative=False, color='r',label='Disk')
    plt.ylabel('Number of stars')
    plt.xlabel('Mass [$M_{sun}$]')
    plt.legend(loc='upper right')
    plt.xlim(xmin=0)
    plt.title('Mass histogram for {name}'.format(name=name_input))
    
    plt.subplot(axes[1])
    axes[1].hist(cluster.Mass[diskless], bins=20, rwidth=0.9, cumulative=False, color='b',label='Diskless')
    plt.xlabel('Mass [$M_{sun}$]')
    plt.legend(loc='upper right')
    plt.xlim(xmin=0)
    
    plt.savefig('{name}_disk_mass_hist.png'.format(name=name_input))
    plt.show()

#plots cumulative histogram of stellar masses
def MassHistCumul():
    fig, axes = plt.subplots(1,2, figsize=(8,4))
    
    plt.subplot(axes[0])
    axes[0].hist(cluster.Mass[disked], bins=20, rwidth=0.9, cumulative=True, color='r',label='Disk')
    plt.ylabel('Number of stars')
    plt.xlabel('Mass [$M_{sun}$]')
    plt.legend(loc='upper right')
    plt.xlim(xmin=0)
    plt.title('Cumulative mass histogram for {name}'.format(name=name_input))
    
    plt.subplot(axes[1])
    axes[1].hist(cluster.Mass[diskless], bins=20, rwidth=0.9, cumulative=True, color='b',label='Diskless')
    plt.xlabel('Mass [$M_{sun}$]')
    plt.legend(loc='upper right')
    plt.xlim(xmin=0)
    
    plt.savefig('{name}_disk_cumul_mass_hist.png'.format(name=name_input))
    plt.show()

#plots normal histogram of stellar periods
def PeriodHist():
    
    bin_list=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28]
    
    plt.hist(cluster.Prot[has_diskless_mass_prot], bins=bin_list, rwidth=0.9, alpha=0.5, edgecolor='blue', cumulative=False, color='c',label='Diskless')
    plt.hist(cluster.Prot[has_disked_mass_prot], bins=bin_list, rwidth=0.9, alpha=0.5, edgecolor='purple', cumulative=False, color='m',label='Disk')
    plt.ylabel('Number of stars')
    plt.xlabel('Period [days]')
    plt.xlim(xmin=0)
    plt.legend(loc='upper right')
    plt.title('Period histogram for {name}'.format(name=name_input))
    plt.savefig('{name}_disk_period_hist.png'.format(name=name_input))
    plt.show()
    
    
def DiskFractPlot():
    
    window=int(20)
    
    prot_sort=cluster.Prot[has_disk_mass_prot]
    sort=np.argsort(prot_sort)
    prot_sort=prot_sort[sort]
    length=int(len(cluster.Prot[has_disk_mass_prot])/window)
    
    bin_list=np.empty(length)
    bin_list[0]=0
    for i in range(0,(length-1)):
        bin_list[i+1]=prot_sort[(i+1)*window]+bin_list[i]
        

    a=plt.hist(cluster.Prot[has_diskless_mass_prot], bins=bin_list, rwidth=0.9, alpha=0.5, edgecolor='blue', cumulative=False, color='c',label='Diskless')
    b=plt.hist(cluster.Prot[has_disked_mass_prot], bins=bin_list, rwidth=0.9, alpha=0.5, edgecolor='purple', cumulative=False, color='m',label='Disk')
    
    plt.show()
   
    #creates an array for the midpont between each bin of the histogram
    x=np.empty(len(bin_list)-1)
    for i in range(0,len(bin_list)-1):
        x[i]=0.5*(bin_list[i]+bin_list[i+1])
    
    y=100*b[0]/(a[0]+b[0]) #calculating disk fraction
    
    #error bars in x axis are the size of the bins used in the histogram
    xerr=np.empty(len(x))
    for i in range(0,len(x)):
        xerr[i]=0.5*(bin_list[i+1]-bin_list[i])
        
    plt.errorbar(x,y,None,xerr,alpha=0.6, color='g', fmt='o')

    plt.ylabel('Disk Fraction [%]')
    plt.xlabel('Period [days]')
    plt.xlim(xmin=0)
    #plt.legend(loc='upper right')
    plt.title('Disk fraction-Period plot for {name}'.format(name=name_input))
    
    #plt.savefig('{name}_disk_fract_period_scatter.png'.format(name=name_input))
    plt.show()
    
#way of plotting the data similar to the rolling quantiles -- helps small sample size problem
def DiskFractAlt():
    
    prot_sort=cluster.Prot[has_disk_mass_prot]
    disk_sort=cluster.Disk[has_disk_mass_prot]
    sort=np.argsort(prot_sort)
    prot_sort=prot_sort[sort]
    disk_sort=disk_sort[sort]
    data_length=int(len(cluster.Prot[has_disk_mass_prot]))
    window=int(data_length*0.5)
    length=data_length-(window+1)
    
    y=np.empty(length)
    x=np.empty(length)
    for i in range(0,(length)):
        n_disk=0
        n_diskless=0
        n_unknown=0
        for n in range(i,(i+window)):
            if disk_sort[n]>0:
                n_disk=n_disk+1
            elif disk_sort[n]<0:
                n_diskless=n_diskless+1
            else:
                n_unknown=n_unknown+1
        x[i]=prot_sort[i]
        y[i]=100*n_disk/(n_disk+n_diskless+n_unknown)
   
    #a lot of the slower rotators are cut off due to the window size, added an extra section to try and keep more data in the plot
    window_slow=int(data_length*0.15) #must be lower than window
    s=np.empty(window-(window_slow+1))
    t=np.empty(window-(window_slow+1))
    

    for p in range(0,window-(window_slow+1)): 
        q_disk=0
        q_diskless=0
        q_unknown=0
        for q in range(p,(p+window_slow)):
            if disk_sort[length+q]>0:
                q_disk=q_disk+1
            elif disk_sort[length+q]<0:
                q_diskless=q_diskless+1
            else:
                q_unknown=q_unknown+1
        s[p]=prot_sort[length+p]
        t[p]=100*q_disk/(q_disk+q_diskless+q_unknown)
        
    w=[*x,*s] #joins the two sets of data together so that the curves can be overplotted to fix discontinuity between curves
    z=[*y,*t]
        
    plt.ylabel('Disk Fraction [%]')
    plt.xlabel('Period [days]')
    plt.title('Rolling Disk fraction-Period plot for {name}'.format(name=name_input))
    plt.plot(w,z,color='b',label='window size = 50%')
    plt.plot(s,t,color='r', label= 'window size = 15%')
    plt.legend(loc='upper right')
    
    plt.savefig('{name}_disk_fract_period_rolling.png'.format(name=name_input))
    plt.show()

def RollingQuant():
    
    prot_sort=cluster.Prot[has_disk_mass_prot]
    mass_sort=cluster.Mass[has_disk_mass_prot]
    
    window_size=int(0.25*len(has_disk_mass_prot))
    min_size=int(0.15*len(has_disk_mass_prot))
    
    print(min_size)

    sort=np.argsort(mass_sort) #creates an array for the positions of cluster.mass if it were in order of smallest to largest

    prot_sort=prot_sort[sort]
    mass_sort=mass_sort[sort]
    
    series=pd.Series(prot_sort)
    
    quant=series.rolling(window=window_size,min_periods=min_size,center=True).quantile(0.5,interpolation='linear')
    plt.plot(mass_sort,quant,color='k',alpha=0.9)
    
    quant=series.rolling(window=window_size,min_periods=min_size,center=True).quantile(0.9,interpolation='linear')
    plt.plot(mass_sort,quant,color='k',alpha=0.9)
    
    quant=series.rolling(window=window_size,min_periods=min_size,center=True).quantile(0.1,interpolation='linear')
    plt.plot(mass_sort,quant,color='k',alpha=0.9)

def PeriodMap():
    
    if name_input == 'NGC2362': #or name_input == 'USco'
        import astropy.coordinates as ac
        cat = ac.SkyCoord(cluster.RA[has_disk_mass_prot], cluster.Dec[has_disk_mass_prot], unit="deg")
        a=plt.scatter(cat.ra.deg*15, cat.dec.deg, c=cluster.Prot[has_disk_mass_prot],alpha=0.6,cmap = plt.cm.Spectral_r) #RA is *15 to convert from time to degrees
       
    else:
        a=plt.scatter(cluster.RA[has_disk_mass_prot], cluster.Dec[has_disk_mass_prot], c=cluster.Prot[has_disk_mass_prot],alpha=0.6,cmap = plt.cm.Spectral_r)
    
    cb=plt.colorbar(a, orientation="vertical", pad=0.01,aspect=15)
    cb.ax.set_ylabel('Period [days]', rotation=270,linespacing=5,fontsize=10,labelpad=20)
    #plt.legend(loc='upper right')    
    plt.xlabel('Right Ascension [deg]')
    plt.ylabel('Declination [deg]')
    plt.title('RA-Dec Map of {name}'.format(name=name_input))
    plt.savefig('{name}_period_ra_dec_map.png'.format(name=name_input))
    
    
    plt.show()
    
"""
============================	CALCULATIONS	=================================
"""

def WidthCalc():
    cluster.ClusterInfo()
    if name_input == 'NGC2362': #or name_input == 'USco'
        import astropy.coordinates as ac
        cat = ac.SkyCoord(cluster.RA, cluster.Dec, unit="deg")
        RA_width = 15*(max(cat.ra.deg)-min(cat.ra.deg))
        DEC_width = max(cat.dec.deg)-min(cat.dec.deg)

    else:
        RA_width = max(cluster.RA)-min(cluster.RA)
        DEC_width = max(cluster.Dec)-min(cluster.Dec)  
        
    print('Distance = ', cluster.dist)
    
    phys_RA_width = math.tan(math.radians(RA_width/2))*cluster.dist*2
    phys_DEC_width = math.tan(math.radians(DEC_width/2))*cluster.dist*2

    #prints width of cluster in deg and in pc
    print('RA width = ', RA_width, phys_RA_width)
    print('Dec width = ', DEC_width, phys_DEC_width) 
    
def KSTest():
    from scipy import stats
    ks_disk=stats.ks_2samp(cluster.Prot[has_disked_mass_prot],cluster.Prot[has_diskless_mass_prot])
    ks_mass=stats.ks_2samp(cluster.Prot[has_disk_high_mass_prot],cluster.Prot[has_disk_low_mass_prot])
    
    print(ks_disk,ks_mass)
  

"""
============================	FUNCTIONS	=================================
"""

#calls functions: comment/uncomment as desired

#PeriodMassPlot()
#PeriodMassLogPlot()
#CoordMap()
#MassHistNorm()
#MassHistCumul()
#PeriodHist()
#DiskFractPlot()
DiskFractAlt()
#MassCompare()
#MassHistCompare()
#WidthCalc()
#KSTest()
#PeriodMap()

