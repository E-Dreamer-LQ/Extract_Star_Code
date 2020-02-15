#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  7 16:00:46 2018

@author: liuqiang
"""

import sys,os,glob ,math,scipy,itertools
import numpy as np
from astropy.io import fits
import bw
from extract_code1 import Filepath,File
from star_analysis import Starimage


##################compute value of snr and remove snrerror###########################
filepath1 = input("input the fits file path:")
#######fitsfound1 represent the exit fit file path we find. 
fitsfound1 = Filepath(filepath1)
######  fitneeded1 represent the exit fits file we need.
fitsneeded1=fitsfound1.findfile()
counts0 = 1
snr_list = []
for t in fitsneeded1:
    s = File(t)
    gray,h,w = s.readfile()
    target = Starimage(gray)
    #compute snr.
    compute_snr = target.SNR()
    snr_list.append(compute_snr)
    


#################### Targets with very low SNR will be filtered out ####################
snrerror = int(input("input a SNR abnormal critical value:"))
############## fits name that I test is 'path+num+1 to XXX.fits'
counts11 = 0           ##########  warning:  must keep fits file name in this format  
for scan in snr_list:
    target.SNRError(scan,filepath1,snrerror,counts11)
    counts11 += 1
    
    
    
    
    
    
    
    