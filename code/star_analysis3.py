#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  7 16:54:39 2018

@author: liuqiang
"""

import sys,os,glob ,math,scipy,itertools
import numpy as np
from astropy.io import fits
import bw
from extract_code1 import Filepath,File
from star_analysis import Starimage



##################compute value of snr and remove snrerror###########################
filepath2 = input("input the fits file path:")
#######fitsfound1 represent the exit fit file path we find. 
fitsfound2 = Filepath(filepath2)
######  fitneeded1 represent the exit fits file we need.
fitsneeded2=fitsfound2.findfile()
ellipticity_list = []
for w in fitsneeded2:
    u = File(t)
    gray,h,w = u.readfile()
    target1 = Starimage(gray)
    #compute ellipticity
    ellipticity = target1.compute_ellipticity()
    #add all the ellipticity to an empty list
    ellipticity_list.append(ellipticity)
    


#################### Targets with very low SNR will be filtered out ####################
ellipticityerror = int(input("input an ellipticity abnormal critical value:"))
############## fits name that I test is 'path+num+1 to XXX.fits'
counts22 = 0           ##########  warning:  must keep fits file name in this format  
for scan1 in ellipticity_list:
    target1.ellipticity_error(scan1,filepath1,ellipticityerror,counts22)
    counts22 += 1
    
    
    