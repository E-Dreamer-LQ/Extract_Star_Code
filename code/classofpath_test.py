#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 20:36:59 2018

@author: liuqiang
"""

import sys,os,glob ,math,scipy,itertools
import numpy as np
from astropy.io import fits




class Filepath():
    def __init__(self,path):
        self.path = path
    
    
    def findfile(self):
        fitfilepath = self.path+'/*.fit*'
        fitfile = glob.glob(fitfilepath)
        return fitfile   
        
filepath = input("please input the fit file path:")
#######a represent the exit fit file path we find. 
a = Filepath(filepath)
######  b represent the exit fits file we need.
b=a.findfile()






        