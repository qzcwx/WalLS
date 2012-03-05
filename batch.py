import os
import pdb
import random
import math
import time


""" remember to generate the instances for MAXSAT problem
    instances = random.sample(range(1,1001), numOfInstance)
"""

# python demo.py [NameOfProblem] [NameOfAlgorithm] [fit/mean/std] [I] [PopSize] [N] [K] [Q]
class Struct:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

def resultExist(probName,algoName,fitName,inst,s,c,n,k,q):
    """ check whether the results have been produced or not """
    if probName == 'NKQ':
        nameOfF = './result/'+probName+'-'+algoName+'-F'+fitName+'-C'+c+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'-K'+str(k)+'-Q'+str(q)+'.txt'
    elif probName == 'NK':
        nameOfF = './result/'+probName+'-'+algoName+'-F'+fitName+'-C'+c+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'-K'+str(k)+'.txt'

    if os.path.isfile(nameOfF)==True:
        print nameOfF, 'exists!!!'
    return os.path.isfile(nameOfF)

def countJobs(p,a,f,i,s,c,n,k,q,count):
    if resultExist(p,a,f,i,s,c,n,k,q) == False:
        return count + 1
    else:
        return count

def writeScript(p,a,f,i,s,c,n,k,q,chunkSize, track):
    if resultExist(p,a,f,i,s,c,n,k,q) == False:
        scriptNo = random.randint(0,len(track)-1)    
        fName = 'run-'+str(track[scriptNo].num)+'.sh'
        fileName = open(fName, 'a')
#        print 'python demo.py '+p+' '+a+' '+f+' '+str(i)+' '+str(s)+' '+str(n)+' '+str(k)+' '+str(q)+' '+'&\n'
        fileName.write('python demo.py '+c+' '+p+' '+a+' '+f+' '+str(i)+' '+str(s)+' '+str(n)+' '+str(k)+' '+str(q)+' '+'&\n')
        fileName.close()

        track[scriptNo].jobs = track[scriptNo].jobs + 1
        if track[scriptNo].jobs == chunkSize:
            # the script file is full
            del track[scriptNo]

def writeHeader(numOfjobs, chunkSize):
    numOfScript = int(math.ceil(numOfjobs/float(chunkSize)))
    print 'header', numOfScript
    for i in range(numOfScript):
        fName = 'run-'+str(i+1)+'.sh'
        fileName = open(fName, 'w')
        print >>fileName, '#!/bin/sh'
        fileName.close()
    # initialize tracking table
    track = []
    for i in range(numOfScript):
        track.append(Struct(num = i+1, jobs = 0))
    return track

def writeFooter(numOfjobs, chunkSize):
    numOfScript = int(math.ceil(numOfjobs/float(chunkSize)))
    print 'footer', numOfScript
    for i in range(numOfScript):
        fName = 'run-'+str(i+1)+'.sh'
        fileName = open(fName, 'a')
        print >>fileName, 'wait'
        fileName.close()
    submitJobs(numOfScript)

def submitJobs(numOfScript):
    """ submit jobs in a scheduled manner """
    submit = False 

    for i in range(numOfScript):
        if submit == False:
            numOfRun = int(os.popen('apstat | grep wchen | wc -l').read())
            maxJobs = int(open('runJobs.conf').read())
            jobToSubmit = maxJobs - numOfRun
            print 'job to submit', jobToSubmit
            submit = True 
        if submit == True:
            if jobToSubmit != 0 :
                print 'Job', i, 'submitted'
                os.system('qsub sat.pbs -v n='+str(i+1))
                jobToSubmit = jobToSubmit - 1
            else :
                submit = False
                time.sleep(10)

if __name__== "__main__":
    count = 0
    chunkSize = 1000
    scriptNo = 0
    fileName = ''
    #numOfInstance = 1
    rseed = 0
    overWrite = 0

    kRange = [2,4,8]
    nRange = [20,50,100,200,500]
    iRange = [0]
    aRange = ['rLS']
    fRange = ['fit','mean']
    cRange = ['hyperWalSearch']
    pRange = ['NK', 'NKQ']

    random.seed(rseed)

    """ the first pass count the number of jobs needed to be run """
    for k in kRange:
        for n in nRange :
            for i in iRange :
                for a in aRange :
                    for f in fRange :
                        for c in cRange :
                            for p in pRange :
                                if p == 'NKQ':
                                    for q in [2]:
                                        if a == 'LS' or a == 'rLS':
                                            s = 1
                                            count =  countJobs(p,a,f,i,s,c,n,k,q,count)
                                        else: # if the algorithm is population-based, for GA and CHC
                                            for s in [30]:
                                                count = countJobs(p,a,f,i,s,c,n,k,q,count)
                                elif p == 'NK': # for NK problem
                                    q = 0
                                    if a == 'LS' or a == 'rLS':
                                        s = 1
                                        count = countJobs(p,a,f,i,s,c,n,k,q,count)
                                    else: # if the algorithm is population-based, for GA and CHC
                                        for s in [30]:
                                            count = countJobs(p,a,f,i,s,c,n,k,q,count)
    totalJobs = count

    # clean up
    os.system('rm run-* SumSat-*')
    track = writeHeader(totalJobs, chunkSize)
    numOfscript = len(track)
    for k in kRange:
        for n in nRange :
            for i in iRange :
                for a in aRange :
                    for f in fRange :
                        for c in cRange :
                            for p in pRange :
                                if p == 'NKQ':
                                    for q in [2]:
                                        if a == 'LS' or a == 'rLS':
                                            s = 1
                                            writeScript(p,a,f,i,s,c,n,k,q,chunkSize, track)
                                        else: # if the algorithm is population-based, for GA and CHC
                                            for s in [30]:
                                                writeScript(p,a,f,i,s,c,n,k,q,chunkSize, track)
                                elif p == 'NK': # for NK problem
                                    q = 0
                                    if a == 'LS' or a == 'rLS':
                                        s = 1
                                        writeScript(p,a,f,i,s,c,n,k,q,chunkSize, track)
                                    else: # if the algorithm is population-based, for GA and CHC
                                        for s in [30]:
                                            writeScript(p,a,f,i,s,c,n,k,q,chunkSize, track)
    writeFooter(totalJobs, chunkSize)

    # for MAXSAT problem
    #instances = random.sample(range(1,1001), numOfInstance)
    #p = 'SAT'
    #for a in ['LS', 'GA', 'CHC']:
    #    for i in instances:
    #        if a == 'LS':
    #            s =1
    #            for f in ['fit', 'mean', 'std']:
    #                count = writeScript(p,a,f,i,s,n,k,q,count,scriptNo,chunkSize, track)
    #
    #        else: # if the algorithm is population-based, for GA and CHC
    #            for s in [30,100]:
    #                for f in ['fit', 'mean', 'std']:
    #                    count = writeScript(p,a,f,i,s,n,k,q,count,scriptNo,chunkSize, track)
    #

    print count, 'jobs submitted', 'separated into', numOfscript, 'scripts'
