#!/usr/bin/env python
# -*- coding:Utf-8 -*-

""" Fits a power law function given a maximum, minimum value, and the number of points 

The code takes the maximum and minimum value and creates a power law fit between the maximum and minimum. 

Examples
----------

To use function, type in the following on the command line:

    python powerlaw.py Emin Emax Number



Arguments
----------

min: numerical 
    minimum value of range

max: numerical 
    maximum value of range    
    
nb: numerical
    number of points to include in the range for the power law fit

Returns
-------

tab: A list or array of values in the range, with the interpolated values calculated using the power law

Notes
--------
    
Modules needed: os, sys, time, scipy, threading, datetime

Notes created by B Hegyi, 10/22/20

"""

import os
import sys
import time
from scipy import *
from threading import Thread
from datetime import datetime

def powerlaw(min,max,nb):

    """Plot wavelength vs. cross-section line plot for species with standard cross section values contained in the xml file
    
    
    Arguments
    ----------
    min: numerical 
        minimum value of range
    
    max: numerical 
        maximum value of range    
        
    nb: numerical
        number of points to include in the range for the power law fit
    
    Returns
    -------
    
    tab: A list or array of values in the range, with the interpolated values calculated using the power law
    
    Notes
    -------
    
    
    """
    eps=1.E-4
    dzo=(max-min)/(2.*nb)
    if(dzo>1.):
        dzo=1.
    flntabm1=(nb-1.)
    hsave=0.
    itr=0

    while(True):
        itr+=1
        hnu=flntabm1*dzo/log( (max+hsave)/(min+hsave) ) - min
        resid=abs(hsave-hnu)
        if(abs(hsave)<abs(hnu)):
            scal = abs(hsave)
        else:
            scal=abs(hnu)
            hsave=hnu
        if(scal!=0):
            resid=resid/scal
        if(resid>eps):
            break


    ratio=(max+hsave)/(min+hsave)
    alrat=log(ratio)

    nbr=int(nb)
    tab=[]
    for i in range(0,nbr):
        tab.append(0)
    for i in range(0,nbr):
        tab[nbr-1-i]=(min+hsave)*pow(ratio,(i*1.)/flntabm1)-hsave
    tab[0]=max
    tab[nbr-1]=min
    return tab
		








if __name__=="__main__":
	print("Hello")
	if(len(sys.argv)<4):
		print("Need more values. Use the following to call the powerlaw: ")
		print("powerlaw Emin Emax Number")
		sys.exit()
	print(sys.argv[1])
	print(sys.argv[2])
	print(sys.argv[3])
	print(powerlaw(float(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3])))


