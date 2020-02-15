#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 25 11:26:41 2018

@author: liuqiang
"""

import sys,os,glob ,math,scipy,itertools
import numpy as np
from astropy.io import fits    #astropy is a package for astro in python



#-------------------------------------------------------
def size(IN):
    """
    function:
        compute the IN's list's length and the list's first element's length.
        
    args:
        IN:the list of IN.  PS:IN is a nested list.
        
    returns:
        M(int):the length of the first element.if IN is a nested list, M represent the number of columns.
        N(int):the length of the list.if IN is a nested list,N represent the number of row.
    """
    
    M=len(IN[0])
    N=len(IN)
    return (M, N)
 
    
def NumberOfRuns(IN):
    """
    function:
        compute the number of Non-zero pixels which next to each other in each column.Run means the Non-zero pixels.
        
    args:
        IN:the nested list we input.
        
    returns:
        result(int):the number of runs.
    """
    
    #M represent the number of columns,N represent the number of rows.
    M,N = size(IN)
    #set the init value of Runs.
    result = 0
    if M != 0 and N != 0:
        #Perform column scanning on the image
        for col in IN:
            #first element is non-zero in a column?
            if  col[0] != 0:
                #if ture,the number of Runs + 1
                result += 1
                #scan all the columns
            for idx in range(1,M):
                #the current element is non-zero and the last element is zero?
                if col[idx] != 0 and col[idx-1] == 0:
                    result += 1
    return result


def FillRunVectors(IN):
    """
    function:
        compute every run's start row,end row and current column.
        
    args:
        IN:the nested list we input.
        
    returns:
        sr, er, c ,three lists
        sr[k] represent the first K Run's start row.
        er[k] represent the first k Run's end row.
        c[k]  represent the first k current column.
    """
    
    M,N = size(IN)
    sr = []
    er = []
    c = []
    #scan the whole image,cidx represent the iteration's time,col represent the element
    for cidx,col in enumerate(IN): 
        #Stand at the front of each column
        k = 0
        #search for the first 1 appear in a column
        while k < M:
            #try-except deal with the abnormity
            try:
                #move to the position of 1
                k += col[k:].index(1)
                c.append(cidx+1)
                sr.append(k+1)
                #next to search for the first 0 appear in a column
                try:
                    #move to the position of 0
                    k += col[k:].index(0)
                except ValueError:
                    #break the recurrent
                    k = M
                er.append(k)
            except ValueError:
                break
    #return every run's start row,end row and current column.
    return sr,er,c
        

def FirstPass(numRuns,mode,sr,er,c):
    """
    function:
        Deal with all of Runs in a image.Define a label to every exit Run.
        
    args:
        numRuns:the number of Runs.
        mode:the type of connected domain.
        sr:Run's start row.
        er:Run's end row.
        c: Run's current column.
        
    returns:
        labels, rowEquivalences, colEquivalences
        labels(list) represent the label which we make in a Run.
        rowEquivalences(list) used for storaging different labels in the same cluster. 
        colEquivalences(list) used for storaging different labels in the same cluster. 
    """
    
    #initialize
    currentColumn = 0
    nextLabel = 1
    firstRunOnPreviousColumn = -1
    lastRunOnPreviousColumn = -1
    firstRunOnThisColumn = -1
    #create a empty list for storaging
    equivList = []
    #first create the number of Runs labels 
    labels=[0] * numRuns 
    #judge the mode of connected domain
    if mode == 8:
        offset = 1
    else:
        offset = 0
    #travelsing all of Runs,we deal with the K Runs in every circulation
    for k in range(numRuns): 
        #if The K Run and The K-1 Run are in adjacent columns
        if c[k] == currentColumn + 1:
            
            firstRunOnPreviousColumn = firstRunOnThisColumn
            #Let firstRunOnThisColumn point to The k Run ,The k Run became the first Run in the column
            firstRunOnThisColumn = k
            #Let lastRunOnPreviousColumn point to the k-1 Run,The k-1 Run became the last Run in the last column
            lastRunOnPreviousColumn = k-1
            #Let currentColumn point to the k Run's column
            currentColumn = c[k]
        #if The K Run and The k-1 Run are not in adjacent columns 
        elif c[k] > currentColumn+1:
            #represent not in adjacent column
            firstRunOnPreviousColumn = -1
            lastRunOnPreviousColumn = -1
            firstRunOnThisColumn = k
            currentColumn = c[k]
        else:
            pass
        #in the situation of  adjacent columns
        if firstRunOnPreviousColumn >= 0:
            p=firstRunOnPreviousColumn
            #deal with all the Runs which next to the current the K Run
            while p <= lastRunOnPreviousColumn and sr[p] <= er[k]+offset:
                #judge if these Runs have the p Run and the current Run overlapped on the row
                if er[k] >= sr[p]-offset and sr[k] <= er[p]+offset:
                    #if it is true,please define the same labels,if it is false,storage in the list of equivList
                    if labels[k] == 0:
                        labels[k] = labels[p]
                    else:
                        if labels[k] != labels[p]:
                            equivList.insert(0, (labels[k],labels[p]))
                        else:
                            pass
                p += 1
        #Label the labels,start from 1
        if labels[k] == 0:
            labels[k] = nextLabel
            nextLabel += 1
    #create two empty list
    rowEquivalences = []
    colEquivalences = []
    if len(equivList) > 0:  
        #the list of equivList's element appear as one pair
        for item0, item1 in equivList:
            rowEquivalences.append(item0)
            colEquivalences.append(item1)
    return labels, rowEquivalences, colEquivalences


def bwlabel1(BW,mode):
    """
    function:
        get the information of a Run.the information include its start row,end row,current column 
        and its labels,its error labels.
        
    args:
        BW:the data's type we input is np.ndarray.To change it into a nested list.
        mode=8:the connected domain's mode is 8. 
        
    returns:
        sr, er, c, labels, rowEquivalences, colEquivalences
    """
    #the iteration in python need to a list,a nested list.
    BW = list(BW)
    
    for k in range(len(BW)):
        BW[k] = list(BW[k])
    numRuns = NumberOfRuns(BW)
    sr, er, c = FillRunVectors(BW)
    labels, rowEquivalences, colEquivalences = FirstPass(numRuns, mode, sr, er, c)
    return sr, er, c, labels, rowEquivalences, colEquivalences


def bwlabel1p5(labels, rowEq, colEq):
    """
    function:
        solve the problem that The Runs belong to the same culster but have different labels.
        
    args:
        labels:the label defined.
        rowEq:storage the error labels.
        colEq:storage the error labels.
        
    returns:
        labels
    """
    
    #create a empty list 
    lblist = []
    #take a node in a for circular
    for k, r in enumerate(rowEq):
        if r == -1:
            continue
        queue = list([rowEq[k], colEq[k]])
        cur_lblist = list()
        #start with a node which can reach at any node to find all the node 
        while queue:
            head = queue.pop(0)
            #package in the cur_lblist
            cur_lblist.append(head) 
            for n in range(k+1, len(rowEq)):
                if rowEq[n] == head:
                    queue.append(colEq[n])
                    #The visited node set -1
                    rowEq[n] = -1
                    colEq[n] = -1
                elif colEq[n] == head:
                    queue.append(rowEq[n])
                    rowEq[n] = -1
                    colEq[n] = -1
        #storage in the lblist
        lblist.append(cur_lblist)
    #transfer into a array
    labels = scipy.array(labels)
    for oldlabels in lblist:
        for ol in oldlabels:
            labels[labels==ol] = oldlabels[0]
    #zip()'s output type is tuple,so changed it into list to sort
    sort_labels = list(zip(labels, range(len(labels))))
    #After find out these Runs which belong to the same culster need to be sorted.
    sort_labels.sort()
    sort_idx = [k[1] for k in sort_labels]
    sort_labels = scipy.array([k[0] for k in sort_labels])
    if sort_labels[0] != 1:
        sort_labels -= (sort_labels[0]-1)
    for k in range(1, len(sort_labels)):
        cur_label = sort_labels[k]
        pre_label = sort_labels[k-1]
        if cur_label > pre_label+1:
            sort_labels[sort_labels==cur_label] = pre_label+1
    for k, l in zip(sort_idx, sort_labels):
        labels[k] = l
    return labels
    

def bwlabel2(sr, er, sc, labels, M, N):
    """
    function:
        show the result at last.
        
    args:
        sr(list):Runs's start row.
        er(list):Run's end row.
        sc(list):Run's current column.
        labels(list):the labels we defined.
        M(int):the number of column.
        N(int):the number of row.
        
    returns:
        L:the last connected domain which have right labels.
    """
    L = scipy.zeros((M, N))
    for s, e, c, l in zip(sr, er, sc, labels):
        L[c-1, (s-1):e] = l
    return L


def bwlabel(BW, mode=8):
    """
    function:
        get the right result by return bwlable2().
        
    args:
        BW:the data's type we input is np.ndarray.To change it into a nested list.
        
    returns:
        bwlabel2(sr, er, sc, labels, N, M)
    """
    sr, er, sc, labels, req, ceq = bwlabel1(BW, mode)
    labels = bwlabel1p5(labels, req, ceq)
    M, N = size(BW)
    return bwlabel2(sr, er, sc, labels, N, M)





