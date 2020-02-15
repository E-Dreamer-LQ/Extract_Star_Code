#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 21:23:19 2018

@author: liuqiang
"""

import sys,os,glob ,math,scipy,itertools
import numpy as np
from astropy.io import fits
import bw
from extract_code1 import Filepath,File
import BasicCal




#create a starimage-class to compute the feature of star image 
class Starimage():
    def __init__(self,img):
        #initialize property,img represent a array
        self.img = img
        
    #Function definition SNR is used to calculate SNR of images
    def SNR(self,sigmacut=1):
        #Ref to lucky imaging paper by Staley
        #https://www.ast.cam.ac.uk/sites/default/files/SPIE_7735-211_010710.pdf
        #Equation 5
        #Data reduction strategies for lucky imaging
        #img is the observational data 2D ndarray
        #sigmacut is the value used to cut background 
        totalflux=np.sum(self.img)
        #Calculate the mean flux of all
        meanflux=np.mean(self.img)
        #Calculate the background with one sigma clipped
        imgsigma=np.var(self.img)**0.5*sigmacut
        #Calculate the background contribution
        backtotal=np.mean(self.img[self.img<=imgsigma+meanflux])*np.shape(self.img)[0]*np.shape(self.img)[1]
        #Calculate the background sigma contribution (background possion noise contribution)
        backsigma=np.var(self.img[self.img<=imgsigma+meanflux])*np.shape(self.img)[0]*np.shape(self.img)[1]
        #Calculate the snr,backdigma is varriance
        return (totalflux-backtotal)/(backsigma)**0.5
        
    def SNRError(self,snr,Path2,snrerror,counts11):
        """
        function:
            remove some abnormal fits file
        args:
            snr(float):the value of snr we compute.
            Path2:fits file exit path.
            snrerror(float):Set the exception threshold.
            counts11(int):count in the fits file name.
        return:
            none.
        """
        if (snr < snrerror):
            os.remove(Path2+'/num'+str(counts11)+'.fits')
        
    def compute_ellipticity(self):
        """
        function:
        compute the fit file's ellipticity.
        args:
            self.img
        returns:
            ELLIPTICITY(int):the value we computed.
        """
        
        ellipse1,ellipse2 = BasicCal.KSBcal(self.img)
        #ellipticity = -+arctan(b/a)  (-+ represent Left-handed right-handed)
        ELLIPTICITY = np.arctan(ellipse2/ellipse1)
        return ELLIPTICITY
    
    def ellipticity_error(self,ellipticity,Path22,ellipticity_error,counts22):
        """
        function:
            remove some abnormal fits file
        args:
            ellipticity(float):the value of snr we compute.
            Path22:fits file exit path.
            ellipticity_error(float):Set the exception threshold.
            counts22(int):count in the fits file name.
        return:
            none.
        """
        if (ellipticity > ellipticity_error):
            os.remove(Path22+'/num'+str(counts22)+'.fits')










        
    

    
        
        
        
        
        
        
        
        
        
        






    




    
    
    
    
            
        
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        






