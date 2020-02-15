#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  7 15:04:17 2018

@author: liuqiang
"""



import sys,os,glob ,math,scipy,itertools
import numpy as np
from astropy.io import fits
import bw
from extract_code1 import Filepath,File


############################### extract star's position and saved as a .txt file  #############################
#filepath = input("please input the fit file path:")
filepath = '/home/liuqiang/test00/test'
#######a represent the exit fit file path we find. 
fitsfound = Filepath(filepath)
######  b represent the exit fits file we need.
fitsneeded=fitsfound.findfile()
counts = 1
####Go through fits file in this folder
for c in fitsneeded:
    d = File(c)
    gray,h,w = d.readfile()
    print(counts,c)
    counts += 1
    #got the array after connected component analysis
    img,copy = d.extractstar(gray,h,w)
    #save these coordinate in a .txt file
    d.saveposition(img,copy)
    

    
