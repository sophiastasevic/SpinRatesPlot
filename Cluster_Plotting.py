# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 10:35:13 2020

@author: sophia

Easier way to plot data from SpinRatesAreHere script. 
Allows input of cluster name and produces the following plots:
    -Period-Mass (linear + log space)
    -RA-Dec map
    -Mass histogram (normal + cumulative)
    -Period histogram (linear + log space)
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
#name_input = 'NGC2362'
name_input = input('Name of cluster:  ')

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
============================	PLOTS	=================================
"""

#plots scatter plot of period vs mass
def PeriodMassPlot():
    plt.scatter(cluster.Mass, cluster.Prot, s=4)
    plt.ylabel('Period [days]')
    plt.xlabel('Mass [$M_{sun}$]')
    plt.title('Period-Mass plot of {name}'.format(name=name_input))
    plt.xlim(xmin=0)
    plt.savefig('{name}_period_mass.png'.format(name=name_input))
    plt.show() 

#same as above but with period on a log scale
def PeriodMassLogPlot():
    plt.scatter(cluster.Mass, cluster.Prot, s=4)
    plt.yscale('log')
    plt.ylabel('Period [days]')
    plt.xlabel('Mass [$M_{sun}$]')
    plt.title('Period-Mass plot of {name} (period log scale)'.format(name=name_input))
    plt.xlim(xmin=0)
    plt.savefig('{name}_period_mass_log.png'.format(name=name_input))
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
    if name_input == 'USco' or name_input == 'NGC2362':
        import astropy.coordinates as ac
        cat = ac.SkyCoord(cluster.RA, cluster.Dec, unit="deg")
        plt.scatter(cat.ra.deg*15, cat.dec.deg, s=4) #RA is *15 to convert from time to degrees
    else:
        plt.scatter(cluster.RA, cluster.Dec, s=4)
        
    plt.xlabel('Right Ascension [deg]')
    plt.ylabel('Declination [deg]')
    plt.title('RA-Dec Map of {name}'.format(name=name_input))
    plt.savefig('{name}_ra_dec_map.png'.format(name=name_input))
    plt.show() 

#plots normal histogram of stellar masses
def MassHistNorm():
    plt.hist(cluster.Mass, bins=20, rwidth=0.9, cumulative=False)
    plt.ylabel('Number of stars')
    plt.xlabel('Mass [$M_{sun}$]')
    plt.xlim(xmin=0)
    plt.title('Mass histogram for {name}'.format(name=name_input))
    plt.savefig('{name}_mass_hist.png'.format(name=name_input))
    plt.show()

#plots cumulative histogram of stellar masses
def MassHistCumul():
    plt.hist(cluster.Mass, bins=20, rwidth=0.9, cumulative=True)
    plt.ylabel('Number of stars')
    plt.xlabel('Mass [$M_{sun}$]')
    plt.xlim(xmin=0)
    plt.title('Cumulative mass histogram for {name}'.format(name=name_input))
    plt.savefig('{name}_cumul_mass_hist.png'.format(name=name_input))
    plt.show()

#plots normal histogram of stellar periods
def PeriodHist():
    plt.hist(cluster.Prot, bins=20, rwidth=0.9, cumulative=False)
    plt.ylabel('Number of stars')
    plt.xlabel('Period [days]')
    #plt.xscale('scalar')
    plt.title('Period histogram for {name}'.format(name=name_input))
    plt.savefig('{name}_period_hist.png'.format(name=name_input))
    plt.show()
    
"""
============================	CALCULATIONS	=================================
"""

def WidthCalc():
    cluster.ClusterInfo()
    if name_input == 'USco' or name_input == 'NGC2362':
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


#calls functions: comment/uncomment as desired
#PeriodMassPlot()
#PeriodMassLogPlot()
#CoordMap()
#MassHistNorm()
#MassHistCumul()
#PeriodHist()
#MassCompare()
#MassHistCompare()
#WidthCalc()