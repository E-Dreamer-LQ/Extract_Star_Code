#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 17:24:35 2018

@author: liuqiang
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 15:05:48 2018

@author: liuqiang
"""
"""
This program is written on the same class as the previous one.
"""

import sys,os,glob ,math,scipy,itertools
import numpy as np
from astropy.io import fits
import background
import bw
import re


###########create a filefolder-class
class Filepath():
    def __init__(self,path):
        self.path = path
    
    
    def findfile(self):
        #find all fit files in a folder
        fitfilepath = self.path+'/*.fit*'
        fitfile = glob.glob(fitfilepath)
        return fitfile 
    def position(self):
        #find the position file in a folder
        positionfilepath = self.path+'/*.txt*'
        positionfile = glob.glob(positionfilepath)
        return positionfile
        
#######create a file-class   
class File(Filepath):
    def __init__(self,filename):
        #initialize the son class's property
        super(Filepath,self).__init__()
        #initialize the father class's property
        self.filename = filename
           
    #the method read fit file
    def readfile(self): 
        #open the fits file
        hdul = fits.open(self.filename)
        #put the gray level imfornation in to matrix
        gray = hdul[0].data
        #get the row number and column number
        h = hdul[0].header['NAXIS1']
        w = hdul[0].header['NAXIS2']
        #transfer gray
        gray = gray.T
        return gray,h,w

    def extractstar(self,gray1,h1,w1):
        global mesh,mode,value
        mesh = int(input("please input the mesh size:"))
        mode = int(input("please input the bw's mode:"))
        value = int(input("please input threshold's influence coefficient:"))
        gray2,copy = background.Mesh(mesh,gray1,h1,w1,value)
        img = bw.bwlabel(gray2)
        return img,copy
        
    def saveposition(self,matrix,orimat):
        #change the name of the saved file 
        savefilepath = self.filename.split("/")[-1]
        Save = input("please input a save file path:")
        if os.path.isdir(Save) == True:
            #remove a exit folder
            os.popen('rm -rf '+Save)
            #make a save dirs and give it the execute permission
            os.makedirs(Save,0o777)
        #while we input a empty folder,it make a path automatically.
        elif Save=='':
            Save = sys.path[0]
        else:
            os.makedirs(Save,0o777)
        
        #write a file
        file = open(Save+'/'+savefilepath+'.txt','w+')
        #compute the max number of lable,the type is np.int64
        num = np.amax(matrix)
        count = 0
        tag=1#init the tag number
        while tag <= (num):
            obj = matrix.copy()
            #scan the matrix's row and column
            for i in range(len(matrix[0])):
                for j in range(len(matrix)):
                    #judge if find the connected domain,taat is to say the axis we index == tag 
                    if obj[i,j] == tag:
                        #Assign the original pixel position
                        obj[i,j] = orimat[i,j]
                        count += 1
                    else:
                        obj[i,j] = 0
        #delete regions too small to extract PSFs
        #remove the noise
            if count >= 2:
                element=np.where(obj==np.max(obj))
                #make sure the position of matrix's max value,xval represent the row value,yval represent the column value.
                xval=element[0][0]+1
                yval=element[1][0]+1
                #write the noise's position into file.
                file.write("%s %s\n"%(xval,yval))
        #   count = 0
            tag += 1
        file.close
        #the flag of ending
        print("position's extract completed!\n")
    
    
    def psfextract(self,gray11,savepath):
        self.filename += '.txt'
        #find the .txt file 
        new = np.loadtxt(self.filename)
        e = []
        f = []
        g = []
        #translate the type of ndarray into a list 
        for m in range(len(new)):
            e.append(new[m])
        for z in e:
            f.append(z)
        #translate the str array into a int array
        new_array = np.array(f,dtype=float)
        f = np.array(new_array,dtype=int)
        #the number of row
        number = new.shape[0]
        #input the slice's size randomly
        slice_size = int(input("input the size of slice: "))
        for j in range(number):
            axis = f[j]
            #axis of array
            x_axis = int(axis[0])
            y_axis = int(axis[1])   
            #slice the entire fits file array
            domain = gray11[(x_axis-slice_size):(x_axis+slice_size),(y_axis-slice_size):(y_axis+slice_size)]
            #save the fits files what we slice 
            #saveposition represent the path of these fits file which we sliced and saved
            if os.path.exists(savepath+'/num'+str(j)+'.fits'):
                os.remove(savepath+'/num'+str(j)+'.fits')
            grey = fits.PrimaryHDU(domain)
            greyHDU = fits.HDUList([grey])
            greyHDU.writeto(savepath+'/num'+str(j)+'.fits')
        #return a fits file path 
        return savepath
    
        
############################### extract star's position and saved as a .txt file  #############################
#filepath = input("please input the fit file path:")
#########a represent the exit fit file path we find. 
#fitsfound = Filepath(filepath)
########  b represent the exit fits file we need.
#fitsneeded=fitsfound.findfile()
##counts = 1
#for c in fitsneeded:
#    d = File(c)
#    gray,h,w = d.readfile()
#    print(counts,c)
#    counts += 1
#    img,copy = d.extractstar(gray,h,w)
#    d.saveposition(img,copy)
#    
    
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



      


















    
    
    
    
    


    
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        