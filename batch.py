import os
import random
import math
import time


# python demo.py [NameOfProblem] [NameOfAlgorithm] [fit/mean/std] [I] [PopSize] [N] [K] [Q]
class Struct:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

def resultExist(probName,algoName,fitName,inst,s,c,n,k,q,w,m,t,e):
    """ check whether the results have been produced or not """
    if probName == 'NKQ':
        nameOfF = './result/'+probName+'-'+algoName+'-F'+fitName+'-M'+m+'-I'+str(inst)+'-S'+str(s)+'-W'+str(w)+'-N'+str(n)+'-K'+str(k)+'-C'+str(c)+'-Q'+str(q)+'-T'+str(t)+'-E'+str(e)+'.txt'
    elif probName == 'NK' or probName == 'NonNK':
        nameOfF = './result/'+probName+'-'+algoName+'-F'+fitName+'-C'+str(c)+'-I'+str(inst)+'-S'+str(s)+'-W'+str(w)+'-N'+str(n)+'-K'+str(k)+'-E'+str(e)+'.txt'

    if os.path.isfile(nameOfF)==True:
        print nameOfF, 'exists!!!'
    return os.path.isfile(nameOfF)

def countJobs(p,a,f,i,s,c,n,k,q,w,m,t,e,count):
    if resultExist(p,a,f,i,s,c,n,k,q,w,m,t,e) == False:
        return count + 1
    else:
        return count

def writeScript(p,a,f,i,s,c,n,k,q,w,m,t,chunkSize, track):
    if resultExist(p,a,f,i,s,c,n,k,q,w,m,t,e) == False:
        scriptNo = random.randint(0,len(track)-1)
        fName = 'run-'+str(track[scriptNo].num)+'.sh'
        fileName = open(fName, 'a')
        # fileName.write('cd ~/sched/workspace/SumSat; nice -n 19 python run.py -c '+str(c)+' -p '+p+' -a '+a+' -f '+f+' -i '+str(i)+' -s '+str(s)+' -w '+str(w)+' -n '+str(n)+' -k '+str(k)+' -q '+str(q)+' -m'+m+' -t '+str(t)+' -e '+str(e)+' '+'&\n')
        fileName.write('cd ~/sched/workspace/SumSat; nice -n 19 python run.py -c '+str(c)+' -p '+p+' -a '+a+' -f '+f+' -i '+str(i)+' -s '+str(s)+' -w '+str(w)+' -n '+str(n)+' -k '+str(k)+' -q '+str(q)+' -m '+m+' -t '+str(t)+' -e '+str(e)+' '+'&\n')
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

def initTrack(numOfjobs, chunkSize):
    numOfScript = int(math.ceil(numOfjobs/float(chunkSize)))

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
    chunkSize = 10000
    scriptNo = 0
    fileName = ''
    #numOfInstance = 1
    rseed = 0
    overWrite = 0

    nRange = [20,50,100,150,200,250,300,350,400,450,500]
    kRange = [2,4]
    
    # nRange = [200]
    # kRange = [2]

    # cRange = [1000, 3000, 4260, 6000, 7500, 9000]
    # cRange = [500, 1500, 2130, 3000, 3750, 4500]
    vRange = [1]

    # eRange = [10,20,100]
    eRange = [0]

    tRange = [7]

    iRange = [0]
    aRange = ['rLS']
    fRange = ['fit']
    mRange = ['walRestFlip','walRest','partEval','countFreq']
    
    pRange = ['NK','NonNK']

    # temp = [5*a for a in range(1,11)]
    # temp.insert(0,1)
    # wRange = temp  # [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
    #wRange = [1, 5]
    wRange = [1]

    random.seed(rseed)

    """ the first pass count the number of jobs needed to be run """
    for k in kRange:
        for n in nRange :
            for i in iRange :
                for e in eRange:
                    for a in aRange :
                        for m in mRange:
                            for f in fRange :
                                for t in tRange:
                                    # for c in cRange :
                                    for v in vRange:
                                        c = int(n*v)
                                        for p in pRange :
                                            if p == 'NKQ':
                                                for q in [2]:
                                                    if a == 'LS' or a == 'rLS':
                                                        for w in wRange :
                                                            s = 1
                                                            count =  countJobs(p,a,f,i,s,c,n,k,q,w,m,t,e,count)
                                                    else: # if the algorithm is population-based, for GA and CHC
                                                       for s in [30]:
                                                           count = countJobs(p,a,f,i,s,c,n,k,q,1,m,t,e,count)
                                            elif p == 'NK' or p == 'NonNK': # for NK problem
                                                q = 0
                                                if a == 'LS' or a == 'rLS':
                                                    for w in wRange :
                                                        s = 1
                                                        count =  countJobs(p,a,f,i,s,c,n,k,q,w,m,t,e,count)
                                                else: # if the algorithm is population-based, for GA and CHC
                                                    for s in [30]:
                                                        count = countJobs(p,a,f,i,s,c,n,k,q,1,m,t,e,count)

    totalJobs = count

    # clean up
    os.system('rm run-*')
    os.system('rm  SumSat-*')
    # writeHeader(totalJobs, chunkSize)
    track = initTrack(totalJobs, chunkSize)
    numOfscript = len(track)


    for k in kRange:
        for n in nRange :
            for i in iRange :
                for e in eRange:
                    for a in aRange :
                        for m in mRange:
                            for f in fRange :
                                for t in tRange:
                                    # for c in cRange :
                                    for v in vRange:
                                        c = int(n*v)
                                        for p in pRange :
                                            if p == 'NKQ':
                                                for q in [2]:
                                                    if a == 'LS' or a == 'rLS':
                                                        for w in wRange :
                                                            s = 1
                                                            writeScript(p,a,f,i,s,c,n,k,q,w,m,t,chunkSize, track)
                                                    else: # if the algorithm is population-based, for GA and CHC
                                                        for s in [30]:
                                                            writeScript(p,a,f,i,s,c,n,k,q,1,m,t,chunkSize, track)
                                            elif p == 'NK' or p == 'NonNK': # for NK problem
                                                q = 0
                                                if a == 'LS' or a == 'rLS':
                                                    for w in wRange :
                                                        s = 1
                                                        writeScript(p,a,f,i,s,c,n,k,q,w,m,t,chunkSize, track)
                                                else: # if the algorithm is population-based, for GA and CHC
                                                    for s in [30]:
                                                        writeScript(p,a,f,i,s,c,n,k,q,1,m,t,chunkSize, track)

    # writeFooter(totalJobs, chunkSize)

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

