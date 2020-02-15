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
import bw


###########create a filefolder-class
class Filepath():
    def __init__(self,path):
        self.path = path
        self.fitfile = self.findfile()
        self.positionfile = self.position()
    #define some methods
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
        """
        function:
            read a fit file.
        args:
            self.filename from started class definition.
        returns:
            gray(array):an array from a fit file.
            h(int):array's row number.
            w:array's column number.
        """
        
        #open the fits file
        hdul = fits.open(self.filename)
        #put the gray level imfornation in to matrix
        gray = hdul[0].data
        #get the row number and column number
        h = hdul[0].header['NAXIS1']      # 6274
        w = hdul[0].header['NAXIS2']       #6258  
        
        #transfer gray
        gray = gray.T
        return gray,h,w

    def extractstar(self,gray1,h1,w1):
        """
        function:
            extract star image.
        args:
            gray1(array):an array from a fit file.
            h1(int):array's row number.
            w1:array's column number.
        returns:
            img(array):an array from the connected domain.
            copy(array):copyed an array from a fit file.
        """
        
        
        global mesh,mode,value
#        mesh = int(input("please input the mesh size:"))
#        mode = int(input("please input the bw's mode:"))
#        value = int(input("please input threshold's influence coefficient:"))
        mesh = 64
        mode = 8
        value = 3
#        gray2,copy = background.Mesh(mesh,gray1,h1,w1,value)
        #parameters's initialization
        mesh_h=0
        mesh_w=0
        i=0
        j=0
        #make a copy in gray 
        copy = gray1.copy()
        #get the value of backsize
        mesh_h = mesh+j
        mesh_w = mesh+i
        #create a local area
        while mesh_h <= h1:
            while mesh_w <= w1:
                local_area = gray1[i:mesh_w,j:mesh_h]
                #compute the local area's all the pixels's mean,median,standard diviation
                mean = np.mean(local_area)
                median = np.median(local_area)
                std_d = np.std(local_area)
                #set a threshold to remove some pixels
                threshold = median+std_d*value
                #scan the local area's pixels to find useful pixels
                for x in range(i,mesh_w):
                    for y in range(j,mesh_h):
                        #judge 
                        x_axis=x%mesh
                        y_axis=y%mesh
                        #remove the pixels which fall on  median+3*std_d outside,create a 0-1 array
                        if local_area[x_axis,y_axis] > threshold:
                            gray1[x,y]=1
                        else:
                            gray1[x,y]=0
                i += mesh
                mesh_w += mesh
            i = 0
            mesh_w = i+mesh
            j += mesh
            mesh_h += mesh
        img = bw.bwlabel(gray1)
        return img,copy
        
    def saveposition(self,matrix,orimat):
        """
        function:
            save the position file.
        args:
            matrix(array):an array after dealed with.
            orimat:original fit file's array.
        returns:
            none.
        """
        
        
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
        """
        function:
            Open an coordinate file and Slice the entire fits file around a coordinate point.
        args:
            gray11(array):fits file's array we transport into.
            savepath:The preservation path of fits generated by slices.
        returns:
            savepath:The preservation path of fits generated by slices.
        """
        
        
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
    
        






      


















    
    
    
    
    


    
        

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        