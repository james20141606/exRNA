#! /usr/bin/env python
from tqdm import tqdm
import argparse
import numpy as np
import time
import sys,os,errno,gc
import numba
import itertools
import h5py
import pandas as pd
parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='ind')  
parser.add_argument('-o', dest='output_file')
args = parser.parse_args()


originalnames = np.array(['length','lncRNA', 'miRNA', 'mRNA', 'piRNA', 'snoRNA', 'snRNA', 'srpRNA',
       'tRNA', 'vaultRNA', 'Y_RNA', 'tucpRNA'])
rnanames = np.array(['length', 'miRNA', 
            'Y_RNA', 'piRNA', 'srpRNA', 
            'snoRNA', 'snRNA', 'tRNA','vaultRNA',
            'tucpRNA', 'lncRNA','mRNA' ])
xsorted = np.argsort(originalnames)
ypos = np.searchsorted(originalnames[xsorted], rnanames)
transind  = xsorted[ypos]
wholenum = np.math.factorial(1)*np.math.factorial(3)*np.math.factorial(4)*np.math.factorial(3)

rankin1 = np.array(list(itertools.permutations(np.arange(0,3))))
rankin2 = np.array(list(itertools.permutations(np.arange(0,4))))
rankin3 = np.array(list(itertools.permutations(np.arange(0,3))))
rankind = np.zeros([wholenum,10]).astype('int')
for i in range(6):
    for j in range(24):
        for t in range(6):
            rankind[i*144+j*6+t] =\
            np.concatenate((rankin1[i],rankin2[j],rankin3[t]))

lengthind = np.arange(16,51).astype('int')

filename = np.loadtxt('lengthfilename.txt'
                      ,dtype='str')[int(args.ind)]  
samplearrmx = np.array(pd.read_table(filename)
                       .iloc[:,:])[:,transind]


print ('sample matrix loaded')
print (samplearrmx.shape)


mirnacountafter = np.where(samplearrmx[:,1]==0) 
mirnacountneed = np.where(samplearrmx[:,1]!=0) 
countind1,counts1 =np.unique(samplearrmx[mirnacountneed,0],return_counts=True)
countind1 = countind1.astype('int')
positionind = np.where(np.in1d(lengthind,countind1))
mirnalength= np.zeros([35])
print (counts1.shape)
print (positionind[0].shape)
mirnalength[positionind] = counts1


group1needarr = samplearrmx[mirnacountafter[0]][:,2:5]#np.concatenate((np.array([0]),np.arange(2,5)))]
group1countafter = np.where(np.sum(group1needarr,axis=1)==0)
group1needind = np.where(np.sum(group1needarr,axis=1)!=0)
group1needarr = group1needarr[group1needind]
lengthgroup1 = samplearrmx[:,0][mirnacountafter[0]][group1needind]

group2needarr = samplearrmx[mirnacountafter[0]][group1countafter][:,5:9]#np.concatenate((np.array([0]),np.arange(5,9)))]
group2countafter = np.where(np.sum(group2needarr,axis=1)==0)
group2needind = np.where(np.sum(group2needarr,axis=1)!=0)
group2needarr = group2needarr[group2needind]
lengthgroup2 = samplearrmx[:,0][mirnacountafter[0]][group1countafter[0]][group2needind]

group3needarr = samplearrmx[mirnacountafter[0]][group1countafter][group2countafter][:,9:12]#np.concatenate((np.array([0]),np.arange(9,12)))]
group3needind = np.where(np.sum(group3needarr,axis=1)!=0)
group3needarr = group3needarr[group3needind]
lengthgroup3 = samplearrmx[:,0][mirnacountafter[0]][group1countafter[0]][group2countafter[0]][group3needind]



def get_single_column_length(lengthgroup,groupneedarr,i,rank):
    groupcol = (groupneedarr[:,rank] !=0).argmax(axis=1)
    groupcountind = np.unique( lengthgroup[np.where(groupcol ==i)],return_counts=True)
    groupcountind_ = groupcountind[0].astype('int')
    returnlength = np.zeros([35])
    for i in range(35):
        returnlength[i] = groupcountind[1][np.where(groupcountind_ ==(i+16))] if np.where(groupcountind_ ==(i+16))[0].shape[0]!=0 else 0
    #positionind = np.where(np.in1d(lengthind,groupcountind_))[0]
    #returnlength = np.zeros([35])
    #returnlength[positionind] = groupcountind[1]
    return returnlength

print (group1needarr.shape,group2needarr.shape,group3needarr.shape)

counts1 = np.ndarray([6,3,35]).astype('int')
for i in range(6):
    group1countlength1 = get_single_column_length(lengthgroup1,group1needarr,0,rankin1[i])
    group1countlength2 = get_single_column_length(lengthgroup1,group1needarr,1,rankin1[i])
    group1countlength3 = get_single_column_length(lengthgroup1,group1needarr,2,rankin1[i])
    counts1[i] = np.concatenate((group1countlength1,group1countlength2
                               ,group1countlength3)).reshape(3,-1)

counts2 = np.ndarray([24,4,35]).astype('int')
for i in range(24):
    group2countlength1 = get_single_column_length(lengthgroup2,group2needarr,0,rankin2[i])
    group2countlength2 = get_single_column_length(lengthgroup2,group2needarr,1,rankin2[i])
    group2countlength3 = get_single_column_length(lengthgroup2,group2needarr,2,rankin2[i])
    group2countlength4 = get_single_column_length(lengthgroup2,group2needarr,3,rankin2[i])
    counts2[i] = np.concatenate((group2countlength1,group2countlength2
                               ,group2countlength3,group2countlength4)).reshape(4,-1)
counts3 = np.ndarray([6,3,35]).astype('int')
for i in range(6):
    group3countlength1 = get_single_column_length(lengthgroup3,group3needarr,0,rankin3[i])
    group3countlength2 = get_single_column_length(lengthgroup3,group3needarr,1,rankin3[i])
    group3countlength3 = get_single_column_length(lengthgroup3,group3needarr,2,rankin3[i])
    counts3[i] = np.concatenate((group3countlength1,group3countlength2
                               ,group3countlength3)).reshape(3,-1)

wholecounts = np.zeros([864,11,35]).astype('int')

for i in range(6):
    for j in range(24):
        for t in range(6):
            wholecounts[i*144+j*6+t] =\
            np.concatenate((mirnalength.reshape(1,-1),counts1[i],counts2[j],counts3[t]),axis=0)


print ('whole array dim: ')
print (wholecounts.shape)

def get_rnanameind():
    wholenum = 864
    rankin1 = np.array(list(itertools.permutations(np.arange(1,4))))
    rankin2 = np.array(list(itertools.permutations(np.arange(4,8))))
    rankin3 = np.array(list(itertools.permutations(np.arange(8,11))))
    rankind = np.zeros([wholenum,11]).astype('int')
    for i in range(6):
        for j in range(24):
            for t in range(6):
                rankind[i*144+j*6+t,1:] =\
                np.concatenate((rankin1[i],rankin2[j],rankin3[t]))
    return rankind
adjustarr = np.zeros([864,11,35])

transindtest = np.ndarray([864,11]).astype('int')
for i in range(864):
    originalid = get_rnanameind()[i]
    xsorted = np.argsort(originalid)
    ypos = np.searchsorted(originalid[xsorted], np.arange(0,11))
    transindtest[i] = xsorted[ypos]

for i in range(864):
    adjustarr[i] = wholecounts[i,transindtest[i]]

with h5py.File('02.matrix/new/length/adjust2/'+filename.split('/')[-1]+'lengthbyseq') as f:
    f.create_dataset('counts',data = adjustarr)


'''

{
inds=$(seq 0 1 60)
for ind in $inds; do
echo bin/mappingmxlength.py \
-i ${ind} 
done
}> Jobs/mappingwhole.txt
qsubgen -n mappingwhole -q Z-LU -a 1-60 -j 1 --bsub --task-file Jobs/mappingwhole.txt
bsub < Jobs/mappingwhole.sh
'''