import os
import pdb
import random

""" remember to generate the instances for MAXSAT problem

    instances = random.sample(range(1,1001), numOfInstance)
"""

# python demo.py [NameOfProblem] [NameOfAlgorithm] [fit/mean/std] [I] [PopSize] [N] [K] [Q]


def writeScript(p,a,f,i,s,n,k,q,count,scriptNo,chunkSize,fileName):

    if count % chunkSize == 0:
        scriptNo = scriptNo + 1
        fName = 'run-'+str(scriptNo)+'.sh'
        fileName = open(fName, 'w')
        print >>fileName, '#!/bin/bash\n'

#    pdb.set_trace()
    fileName.write('python demo.py '+p+' '+a+' '+f+' '+str(i)+' '+str(s)+' '+str(n)+' '+str(k)+' '+str(q)+' '+'&\n')

    if count % chunkSize == chunkSize - 1:
        print >>fileName, 'wait\n'
        fileName.close()
        os.system('qsub sat.pbs -v n='+str(scriptNo))
    count = count+1
    return count, scriptNo,fileName

count = 0
chunkSize = 24
scriptNo = 0
fileName = ''
numOfInstance = 10

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
                                    count,scriptNo,fileName = writeScript(p,a,f,i,s,n,k,q,count,scriptNo,chunkSize, fileName)
                                else: # if the algorithm is population-based, for GA and CHC
                                    for s in [30, 100]:
                                        count,scriptNo,fileName = writeScript(p,a,f,i,s,n,k,q,count,scriptNo,chunkSize, fileName)
                        elif p == 'NK': # for NK problem
                            q = 0
                            if a == 'LS':
                                s = 1
                                count,scriptNo,fileName = writeScript(p,a,f,i,s,n,k,q,count,scriptNo,chunkSize, fileName)
                            else: # if the algorithm is population-based, for GA and CHC
                                for s in [30, 100]:
                                    count,scriptNo,fileName = writeScript(p,a,f,i,s,n,k,q,count,scriptNo,chunkSize, fileName)

print 'count', count
# for MAXSAT problem
instances = random.sample(range(1,1001), numOfInstance)
p = 'SAT'
for a in ['LS', 'GA', 'CHC']:
    for i in instances:
        if a == 'LS':
            s =1
            for f in ['fit', 'mean', 'std']:
                count,scriptNo,fileName = writeScript(p,a,f,i,s,n,k,q,count,scriptNo,chunkSize, fileName)

        else: # if the algorithm is population-based, for GA and CHC
            for s in [30,100]:
                for f in ['fit', 'mean', 'std']:
                    count,scriptNo,fileName = writeScript(p,a,f,i,s,n,k,q,count,scriptNo,chunkSize, fileName)

print count, 'jobs submitted', 'separated into', scriptNo, 'scripts'
