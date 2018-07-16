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
parser.add_argument('-i', dest='ind')  #0~9
parser.add_argument('-o', dest='output_file')
#parser.add_argument('-r', dest='indind')
args = parser.parse_args()


originalnames = np.array(['lncRNA', 'miRNA', 'mRNA', 'piRNA', 'snoRNA', 'snRNA', 'srpRNA',
       'tRNA', 'vaultRNA', 'Y_RNA', 'tucpRNA'])
rnanames = np.array([ 'miRNA', 
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


filename = np.loadtxt('filename.txt',dtype='str')[int(args.ind)]  #'17402567-B.all.mx'
samplearrmx = np.array(pd.read_table(filename)
                       .iloc[:,1:])[:,transind]
#/Share/home/caojingyi/exRNA/process/16.one_map/02.matrix
print ('sample matrix loaded')
print (samplearrmx.shape)
start = time.clock()

mirnacountafter = np.where(samplearrmx[:,0]==0)  #filter after mirna
mirnacount = samplearrmx.shape[0] - mirnacountafter[0].shape[0]

group1needarr = samplearrmx[mirnacountafter[0]][:,1:4]
group1countafter = np.where(np.sum(group1needarr,axis=1)==0)
group1needind = np.where(np.sum(group1needarr,axis=1)!=0)
group1needarr = group1needarr[group1needind]

group2needarr = samplearrmx[mirnacountafter[0]][group1countafter][:,4:8]
group2countafter = np.where(np.sum(group2needarr,axis=1)==0)
group2needind = np.where(np.sum(group2needarr,axis=1)!=0)
group2needarr = group2needarr[group2needind]

group3needarr = samplearrmx[mirnacountafter[0]][group1countafter][group2countafter][:,8:11]
group3needind = np.where(np.sum(group3needarr,axis=1)!=0)
group3needarr = group3needarr[group3needind]

print (group1needarr.shape,group2needarr.shape,group3needarr.shape)
counts1 = np.ndarray([6,3]).astype('int')
for i in range(6):
    uniind1 = np.unique(( group1needarr[:,rankin1[i]] !=0).argmax(axis=1),return_counts=True)
    counts1[i][uniind1[0]]= uniind1[1]
counts2 = np.ndarray([24,4]).astype('int')
for i in range(24):
    uniind2 = np.unique(( group2needarr[:,rankin2[i]] !=0).argmax(axis=1),return_counts=True)
    counts2[i][uniind2[0]]= uniind2[1]
counts3 = np.ndarray([6,3]).astype('int')
for i in range(6):
    uniind3 = np.unique(( group3needarr[:,rankin3[i]] !=0).argmax(axis=1),return_counts=True)
    counts3[i][uniind3[0]]= uniind3[1]       
    
wholecounts = np.zeros([wholenum,10]).astype('int')

for i in range(6):
    for j in range(24):
        for t in range(6):
            wholecounts[i*144+j*6+t,:] =\
            np.concatenate((counts1[i],counts2[j],counts3[t]))
wholecounts = np.concatenate((np.repeat(mirnacount,wholenum).reshape(-1,1),wholecounts),axis=1)
elapsed = (time.clock() - start)
print("Time used:",elapsed)

with h5py.File('02.matrix/new/counts/5.25/'+filename.split('/')[-1]+'countsbyseq') as f:
    f.create_dataset('counts',data = wholecounts)


'''
{
ind=$(seq 0 1 2)
indind=$(seq 0 1 3)
for ind in $ind;do
for indind in $indind;do
bin/mappingmx.py \
-i ${ind} \
-r ${indind}
done
done
} | parallel -P 10


{
inds=$(seq 0 1 60)
for ind in $inds; do
echo bin/mappingmx.py \
-i ${ind} 
done
}> Jobs/mappingwhole.txt
qsubgen -n mappingwhole -q Z-LU -a 1-60 -j 1 --bsub --task-file Jobs/mappingwhole.txt
bsub < Jobs/mappingwhole.sh

#bpeek -J "mappingwhole[2]"


{
HDF5_USE_FILE_LOCKING=FALSE bin/mappingmx.py \
    -i 12 
}
'''

'''
rnanames = np.array(['miRNA', 'piRNA', 'Y_RNA', 'snRNA','srpRNA','tRNA',
                     'lncRNA','mRNA','other_genomic_region'])
needrank = int(args.rankind)
rankind = np.zeros([np.math.factorial(needrank),needrank+1]).astype('int')
rankind[:,1:] = np.array(list(itertools.permutations(np.arange(1,needrank+1))))
rangeofrank = int(rankind.shape[0]/10)
indset = {}
for i in range(9):
    indset[i] = np.arange(i*rangeofrank,(i+1)*rangeofrank)
indset[9] = np.arange(9*rangeofrank,rankind.shape[0])
#row = 5000000
#testarr = np.hstack((np.ones(np.floor(row*needrank*0.3).astype('int') ),
                    # np.zeros(row*needrank - np.floor(row*needrank*0.3).astype('int'))))
#np.random.shuffle(testarr)
#testarr = testarr.reshape(-1,needrank)
#mirnaarr = np.hstack(( np.ones(np.floor(row*0.3).astype('int')),
                   #   np.zeros(row - np.floor(row*0.3).astype('int')))).reshape(-1,1)
#np.random.shuffle(mirnaarr)
#wholearr = np.concatenate((mirnaarr,testarr),axis=1)

def get_ratio(ind):
    wholearr[:,rankind[ind]]
    counts = np.unique((wholearr[:,rankind[ind]]!=0).argmax(axis=1),return_counts=True)[1]
    #ratio = counts/np.sum(counts)
    #rnaname = rnanames[rankind[ind]]
    return counts

with h5py.File(args.output_file+args.ind) as f:
    for i in tqdm(indset[int(args.ind)]):
        f.create_dataset(str(i),data = get_ratio(i))

rangeofrank = int(rankind.shape[0]/4)
indset = {}
for i in range(4):
    indset[i] = np.arange(i*rangeofrank,(i+1)*rangeofrank)

'''