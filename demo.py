import nkLandscape as nk
import nkqLandscape as nkq
import WalshAnalysis as wal
import geneticAlgorithm as ga
import AutoCorrelation as ac
#import matplotlib.pyplot as plt
import LocalOptima as lo
import CHC as chc
import MAXSAT as mx
import LocalSearch as ls
import LocalOptima as lo
import tool as tl
import os
import numpy as np
import random
import math
import time
import pdb
import sys
import batch

""" consider as a minimization problem """
tl.checkParam(sys.argv)

rseed = 0
nameOfDir = './result/'
runtimeDir = './runtime/'
waltimeDir = './walshtime/'
prefixNK = './benchmark/NK/'
prefixNKQ = './benchmark/NKQ/'

random.seed(rseed)

compMeth = tl.getArgv() # bf(brute force) / wal (walsh analysis)
probName = tl.getArgv()
algoName = tl.getArgv()
fitName = tl.getArgv() # fit/mean/std
overWrite = int(tl.getArgv())

#if compMeth == 'wal' and fitName != 'mean':
#    print 'ERROR: Walsh analysis can only be applied to compute mean'
#    sys.exit()

inst = int(tl.getArgv())
s = tl.getArgv() # get the setting for population size
n = int(tl.getArgv())
if probName != 'SAT':
    k = int(tl.getArgv())

maxFit = 1000 * n
#maxFit = 0
runs = 30
popSize = 50 # always keep popSize to even number
q = 0

#maxFit = 1000
#runs = 20
#popSize = 4

crossoverR = 0.8 # typically in (0.6, 0.9)
mutationR = 1.0/float(n) # typically between 1/popSize and 1/dim
# for CHC 
D = n/4.0
DR = 0.35
M = 1

print 'probName', probName, 'inst', inst, 'n', n, 'k', k 

if algoName.find('LS') != -1:
    popSize = 1

if probName == 'SAT':
    """ with SAT, we are forced to set n to 100 """

    """ 
    TODO : 
        need to perform multiple runs for each instance 
    """

    model = mx.MAXSAT()
    res = []

    model.setInstance(inst)
    print 'Instance', inst

    if algoName.find('LS') != -1:
        algo = ls.LocalSearch(model.compFit, maxFit, n)
    elif algoName.find('GA') != -1:
        algo = ga.GeneticAlgorithm( model.compFit, maxFit, popSize, n )
    elif algoName.find('CHC') != -1:
        algo = chc.CHC()

    tAll = np.zeros(runs)

    for i in range(runs):
        start = os.times()[0]
        if algoName.find('GA') != -1:
            res.append(algo.run(crossoverR, mutationR, fitName, minimize = False))
        elif algoName == 'LS':
            res.append(algo.run(fitName, minimize = False, restart = False))
        elif algoName == 'rLS':
            res.append(algo.run(fitName, minimize = False, restart = True))
        elif algoName.find('CHC') != -1:
            res.append(algo.run(model.compFit, maxFit, popSize, n, D, DR, M, fitName, minimize = False))
        tAll[i] = os.times()[0] - start

    if probName == 'SAT':
        nameOfF = nameOfDir+probName+'-'+algoName+'-F'+fitName+'-C'+compMeth+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'.txt'
    f = open(nameOfF, 'w')
    for i in range(len(res)):
        if fitName != 'fit':
            print >>f,"%g\t%g\t%g" % (res[i]['sol'], res[i]['fitG'], res[i]['nEvals'])
        else:
            print >>f,"%g\t%g" % (res[i]['sol'], res[i]['nEvals'])
    f.close()

    """ store runtime to files """
    if probName == 'SAT':
        nameOfF = runtimeDir+probName+'-'+algoName+'-F'+fitName+'-C'+compMeth+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'.txt'

    f = open(nameOfF, 'w')
    for i in range(runs):
        print >>f,"%g" % (tAll[i])
    f.close()

else:
    res = []

    if probName == 'NK':
        if overWrite == 1 or (overWrite == 0 and not batch.resultExist(probName,algoName,fitName,inst,popSize,compMeth,n,k,q)):
            model = nk.NKLandscape(n,k,prefixNK+'NK-N'+str(n)+'-K'+str(k)+'-I'+str(inst))
            #model = nk.NKLandscape(n,k)
        else:
            sys.exit()
    elif probName == 'NKQ':
        q = int(tl.getArgv())
        if overWrite == 1 or (overWrite == 0 and not batch.resultExist(probName,algoName,fitName,inst,popSize,compMeth,n,k,q)):
            model = nkq.NKQLandcape(n, k, q, prefixNKQ+'NKQ-N'+str(n)+'-K'+str(k)+'-I'+str(inst)+'-Q'+str(q))
            #model = nkq.NKQLandcape(n, k, q)
        else:
            sys.exit()

    if compMeth == 'walWalk' or compMeth == 'walRest' or compMeth == 'supm' or compMeth == 'bitImp' or compMeth == 'walSearch' or compMeth == 'checkOptWal' or compMeth == 'checkHyper' or compMeth == 'checkHyperRank' or compMeth == 'hyperSearch' or compMeth == 'hyperSqSearch' or 'hyperWalSearch':
        start = os.times()[0]
        w = model.WalshCofLinearLinklist()
        # print Walsh coefficients
        walTime = os.times()[0] - start

        """ store runtime to files """
        if probName == 'NKQ':
            nameOfF = waltimeDir+probName+'-'+algoName+'-F'+fitName+'-C'+compMeth+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'-K'+str(k)+'-Q'+str(q)+'.txt'
        elif probName == 'NK':
            nameOfF = waltimeDir+probName+'-'+algoName+'-F'+fitName+'-C'+compMeth+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'-K'+str(k)+'.txt'

        f = open(nameOfF, 'w')
        print >>f,"%g" % (walTime)
        f.close()

#    bit,fit = tl.compFit(model)
#    a = sorted(zip(bit,fit), key=lambda a_entry: a_entry[1]) 
#    print 'opti\n', a[0][0], a[0][1]
#    print

#    for i in a:
#        print i[0], '%.2f' %(i[1])

#    for i in zip(bit,fit):
#        print i[0],'%.3f' %(i[1])
        
#    print 'bit',bit
#    print 'fit',fit
#    print 'mean',np.mean(fit)
#    print 'w', w

#    numOpt = lo.localOpt(bit, fit)
#    print numOpt

    if algoName.find('LS') != -1:
        algo = ls.LocalSearch(model, maxFit, n)
    elif algoName.find('GA') != -1:
        algo = ga.GeneticAlgorithm( model.compFit, maxFit, popSize, n )
    elif algoName.find('CHC') != -1:
        algo = chc.CHC()

    tAll = np.zeros(runs)
    for i in range(runs):
#        print 'run', i, ':probName', probName, 'algoName', algoName, 'fitName', fitName, 'I', inst, 'n', n, 'k', k 
        start = os.times()[0]
        if algoName.find('GA') != -1:
            res.append(algo.run(crossoverR, mutationR, fitName))
        elif algoName == 'LS':
            res.append(algo.run(fitName, minimize = True, restart = False,compM = compMeth ))
        elif algoName == 'rLS':
            res.append(algo.run(fitName, minimize = True, restart = True,compM = compMeth))
        elif algoName.find('CHC') != -1:
            res.append(algo.run(model.compFit, maxFit, popSize, n, D, DR, M, fitName))
        tAll[i] = os.times()[0] - start

#    trace = res[0]['trace']
#    for i in trace:
##        print 'Eval', i.fitEval, 'fit', i.fit
#        print 'Eval', i.fitEval, 'fit', i.fit, 'fitG', i.fitG
#
#    plt.plot([i.fitEval for i in trace],[i.fit for i in trace],'.-')
#    plt.plot([i.fitEval for i in trace],[i.fitG for i in trace],'.-')
#    plt.show()

    """ store results to files """
    if probName == 'NKQ':
        nameOfF = nameOfDir+probName+'-'+algoName+'-F'+fitName+'-C'+compMeth+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'-K'+str(k)+'-Q'+str(q)+'.txt'
    elif probName == 'NK':
        nameOfF = nameOfDir+probName+'-'+algoName+'-F'+fitName+'-C'+compMeth+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'-K'+str(k)+'.txt'

#    """ print the mean over multiple runs """
#    r = np.zeros(runs)
#    for i in range(runs):
#        r[i] = res[i]['sol']
#    print np.mean(r)
    f = open(nameOfF, 'w')
    for i in range(runs):
        if fitName != 'fit':
            print >>f,"%g\t%g\t%g" % (res[i]['sol'], res[i]['fitG'], res[i]['nEvals'])
        else:
            print >>f,"%g\t%g" % (res[i]['sol'], res[i]['nEvals'])
    f.close()

    """ store runtime to files """
    if probName == 'NKQ':
        nameOfF = runtimeDir+probName+'-'+algoName+'-F'+fitName+'-C'+compMeth+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'-K'+str(k)+'-Q'+str(q)+'.txt'
    elif probName == 'NK':
        nameOfF = runtimeDir+probName+'-'+algoName+'-F'+fitName+'-C'+compMeth+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'-K'+str(k)+'.txt'

    f = open(nameOfF, 'w')
    for i in range(runs):
        if compMeth == 'hyperSearch' or compMeth == 'hyperSqSearch':
            print >>f,"%0.2e\t%0.2e\t%0.2e\t%0.2e" % (tAll[i], res[i]['init'],res[i]['update'], res[i]['hyper'])
        else:
            print >>f,"%0.2e\t%0.2e\t%0.2e" % (tAll[i], res[i]['init'],res[i]['update'])
    f.close()

    print
