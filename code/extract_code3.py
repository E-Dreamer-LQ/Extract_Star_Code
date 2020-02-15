#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  7 16:07:08 2018

@author: liuqiang
"""

import sys,os,glob ,math,scipy,itertools
import numpy as np
from astropy.io import fits
import bw
from extract_code1 import Filepath,File


############################## Slice into optional fits matrix and save ##########################
positionpath = input("please input the exiting position file path:")
positionfound = Filepath(positionpath)
positionneeded = positionfound.position()
#####Traverse the.txt file throughout the folder
for l in positionneeded:
    #change the ...fit.txt file into ...fit file ,Remove .txt at the end of each row
    l = re.sub('.txt', '', l)   
    p = File(l)
    #got the fit array :gray
    gray,h,w = p.readfile()   #warning:must promise that the fits file in this folder
    savepath = p.psfextract(gray,positionpath)
    print("\nslice the fit array completed!")
    
    
    
    
    