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
#hPer = spin.hPer()
#NGC2264 = spin.NGC2264()
#USco = spin.USco()

#asks user to input name of cluster they wish to use, not actually working yet
#for now have to input name manually, will get it working later
#name_input = input('Name of cluster:  ')
name_input = 'NGC2264'
cluster = spin.NGC2264()
mass_type=0 #set mass type from SpinRatesAreHere

#plots scatter plot of period vs mass
def PeriodMassPlot():
    plt.scatter(cluster.Prot, cluster.Mass, s=4)
    plt.xlabel('Period [days]')
    plt.ylabel('Mass [$M_{sun}$]')
    plt.title('Period-Mass plot of {name}'.format(name=name_input))
    plt.savefig('{name}_period_mass.png'.format(name=name_input))
    plt.show() 

#same as above but with period on a log scale
def PeriodMassLogPlot():
    plt.scatter(cluster.Prot, cluster.Mass, s=4)
    plt.xscale('log')
    plt.xlabel('Period [days]')
    plt.ylabel('Mass [$M_{sun}$]')
    plt.title('(Period-Mass plot of {name}'.format(name=name_input))
    plt.savefig('{name}_period_mass_log.png'.format(name=name_input))
    plt.show() 

#plots RA vs Dec map 
#RA and Dec of USco is in xx xx xx form rather than just degrees so don't use for that cluster yet 
def CoordMap():
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
    plt.title('Mass histogram for {name}'.format(name=name_input))
    plt.savefig('{name}_mass_hist.png'.format(name=name_input))
    plt.show()

#plots cumulative histogram of stellar masses
def MassHistCumul():
    plt.hist(cluster.Mass, bins=20, rwidth=0.9, cumulative=True)
    plt.ylabel('Number of stars')
    plt.xlabel('Mass [$M_{sun}$]')
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

#calls functions: comment/uncomment as desired
PeriodMassPlot()
PeriodMassLogPlot()
CoordMap()
MassHistNorm()
MassHistCumul()
PeriodHist()

