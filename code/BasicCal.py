# -*- coding: utf-8 -*-
"""
Used for some basic definition calculation
Created on Thu May  7 11:38:43 2018

@author: Peng
"""

import numpy as np

#Function defintion the following function is used to calculate the Ellipicity under KSB Model
#Clip PSF is used to cut part of psf from the orginal psf
def clippsf(orgpsf,psfsize=None):
    #clip psf according to its maximal value to a predefined size
    #psfsize should be odd number
    #orgpsf is the original psf 2D array
    temcoor=np.where(orgpsf==np.max(orgpsf))
    #print temcoor
    cooy=temcoor[0][0]
    coox=temcoor[1][0]
    if psfsize is not None:
        halfsize=psfsize/2
        newpsf=orgpsf[cooy-halfsize-1:cooy+halfsize,coox-halfsize-1:coox+halfsize]
    else:
        halfsize=np.shape(orgpsf)[0]
        maxsize=np.min([cooy,coox,halfsize])
        newpsf=orgpsf[cooy-maxsize:cooy+maxsize+1,coox-maxsize:coox+maxsize+1]
    return newpsf

#Storexy is used to store the x and y distribution matrix nx and ny should be the size of matrix in x and y
def storeXY(nx,ny):
    XX=np.zeros((nx,ny))
    XY=np.zeros((nx,ny))
    YY=np.zeros((nx,ny))
    for i in range(0,nx):
        x=0.5+i-(nx)/2.0
        for j in range(0,ny):
            y=0.5+j-(ny)/2.0
            r=9#sqrt(x*x+y*y)
            if (r<10):
                XX[i,j]=x*x
                XY[i,j]=x*y
                YY[i,j]=y*y
    return XX,XY,YY

#get Quad is used to calculate the quadic of the matrix
def getQuad(img,XX,XY,YY): # returns the (unweighted) quadrupole matrix of an (nx x ny) img 
    #img is the orginal image with size of N M
    #XX XY and YY are the 2D matrix calculated by storeXY
    quad=np.array([[0.0,0.0],[0.0,0.0]])
    mod=0.
    quad[1][0]=(img*XY).sum()
    quad[0][0]=(img*YY).sum()-mod ##### WARNING: X AND Y HAVE BEEN SWAPPED TO
    quad[1][1]=(img*XX).sum()-mod ##### ACCOUNT FOR NUMPY BEING (Y,X)
    quad[0][1]=quad[1][0] 
    return quad     

#Used to calculate the quad parameter in KSB
def polE(quad): # returns the KSB "polarization parameters" defined in KSB Eq 3.2
    e=np.array([[0.0],[0.0]])
    q1=quad[0][0]-quad[1][1]
    q2=2.0*quad[1][0]
    T=(quad[0][0]+quad[1][1]) /2.
    T+=2.*np.sqrt(np.linalg.det(quad))/2.
    e[0]=q1/T
    e[1]=q2/T
    return e

#Main function used to calculate the Ellipisicity
def KSBcal(orgpsf):
    #need to install package of skimage 
    from skimage.measure import regionprops
#    orgpsf=clippsf(orgpsf)
    sizey,sizex=np.shape(orgpsf)
#    if(sizey != sizex):
#        return 0,0
    XX,XY,YY=storeXY(sizex,sizey)
    quad=getQuad(orgpsf,XX,XY,YY)
    ellipse1,ellipse2=polE(quad)
    return ellipse1[0],ellipse2[0]





