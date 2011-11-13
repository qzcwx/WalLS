import os
import pdb
import random
import time


""" remember to generate the instances for MAXSAT problem

    instances = random.sample(range(1,1001), numOfInstance)
"""

# python demo.py [NameOfProblem] [NameOfAlgorithm] [fit/mean/std] [I] [PopSize] [N] [K] [Q]
class Struct:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

def writeScript(p,a,f,i,s,n,k,q,count,scriptNo,chunkSize, track):
    scriptNo = random.randint(0,len(track)-1)    
    fName = 'run-'+str(track[scriptNo].num)+'.sh'
    fileName = open(fName, 'a')
    fileName.write('python demo.py '+p+' '+a+' '+f+' '+str(i)+' '+str(s)+' '+str(n)+' '+str(k)+' '+str(q)+' '+'&\n')
    fileName.close()

    track[scriptNo].jobs = track[scriptNo].jobs + 1
    if track[scriptNo].jobs == chunkSize:
        # the script file is full
        del track[scriptNo]

def writeHeader(numOfjobs, chunkSize):
    numOfScript = numOfjobs/chunkSize
    for i in range(numOfScript):
        fName = 'run-'+str(i+1)+'.sh'
        fileName = open(fName, 'w')
        print >>fileName, '#!/bin/bash'
        fileName.close()
    # initialize tracking table
    track = []
    for i in range(numOfScript):
        track.append(Struct(num = i+1, jobs = 0))
    return track

def writeFooter(numOfjobs, chunkSize):
    numOfScript = numOfjobs/chunkSize
    for i in range(numOfScript):
        fName = 'run-'+str(i+1)+'.sh'
        fileName = open(fName, 'a')
        print >>fileName, 'wait'
        fileName.close()
    for i in range(numOfScript):
        os.system('qsub sat.pbs -v n='+str(i+1))

count = 0
chunkSize = 24
scriptNo = 0
fileName = ''
numOfInstance = 10
totalJobs = 11400
rseed = 0


random.seed(rseed)
track = writeHeader(totalJobs, chunkSize)
for k in [0, 2, 4, 8, 16]:
    for n in [20, 50, 100] :
        for i in range(10):
            for a in ['LS', 'GA', 'CHC']:
                for f in ['fit', 'mean', 'std']:
                    for p in ['NKQ','NK']:
                        if p == 'NKQ':
                            for q in [2, 4, 8, 16]:
                                if a == 'LS':
                                    s = 1
                                    writeScript(p,a,f,i,s,n,k,q,count,scriptNo,chunkSize, track)
                                else: # if the algorithm is population-based, for GA and CHC
                                    for s in [30, 100]:
                                        writeScript(p,a,f,i,s,n,k,q,count,scriptNo,chunkSize, track)
                        elif p == 'NK': # for NK problem
                            q = 0
                            if a == 'LS':
                                s = 1
                                writeScript(p,a,f,i,s,n,k,q,count,scriptNo,chunkSize, track)
                            else: # if the algorithm is population-based, for GA and CHC
                                for s in [30, 100]:
                                    writeScript(p,a,f,i,s,n,k,q,count,scriptNo,chunkSize, track)

print 'count', count
# for MAXSAT problem
instances = random.sample(range(1,1001), numOfInstance)
p = 'SAT'
for a in ['LS', 'GA', 'CHC']:
    for i in instances:
        if a == 'LS':
            s =1
            for f in ['fit', 'mean', 'std']:
                writeScript(p,a,f,i,s,n,k,q,count,scriptNo,chunkSize, track)

        else: # if the algorithm is population-based, for GA and CHC
            for s in [30,100]:
                for f in ['fit', 'mean', 'std']:
                    writeScript(p,a,f,i,s,n,k,q,count,scriptNo,chunkSize, track)

writeFooter(totalJobs, chunkSize)
print count, 'jobs submitted', 'separated into', scriptNo, 'scripts'
