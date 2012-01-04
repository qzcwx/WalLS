import nkLandscape as nk
import nkqLandscape as nkq
import WalshAnalysis as wal
import geneticAlgorithm as ga
import AutoCorrelation as ac
import matplotlib.pyplot as plt
import CHC as chc
import MAXSAT as mx
import LocalSearch as ls
import LocalOptima as lo
import tool as tl
import numpy as np
import random
import math
import time
import pdb
import sys

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

#if compMeth == 'wal' and fitName != 'mean':
#    print 'ERROR: Walsh analysis can only be applied to compute mean'
#    sys.exit()
inst = int(tl.getArgv())
s = tl.getArgv() # get the setting for population size
n = int(tl.getArgv())
if probName != 'SAT':
    k = int(tl.getArgv())


maxFit = 1000 * n
#maxFit = 10
runs = 1
popSize = 50 # always keep popSize to even number

#maxFit = 1000
#runs = 20
#popSize = 4

crossoverR = 0.8 # typically in (0.6, 0.9)
mutationR = 1.0/float(n) # typically between 1/popSize and 1/dim
# for CHC 
D = n/4.0
DR = 0.35
M = 1

print 'probName', probName, 'algoName', algoName, 'fitName', fitName,  inst, 'n', n, 'k', k 

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
        start = time.time()
        if algoName.find('GA') != -1:
            res.append(algo.run(crossoverR, mutationR, fitName, minimize = False))
        elif algoName == 'LS':
            res.append(algo.run(fitName, minimize = False, restart = False))
        elif algoName == 'rLS':
            res.append(algo.run(fitName, minimize = False, restart = True))
        elif algoName.find('CHC') != -1:
            res.append(algo.run(model.compFit, maxFit, popSize, n, D, DR, M, fitName, minimize = False))
        tAll[i] = time.time() - start

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
    for i in range(len(tAll)):
        print >>f,"%g" % (tAll[i])
    f.close()

else:
    res = []

    if probName == 'NK':
        model = nk.NKLandscape(n,k,prefixNK+'NK-N'+str(n)+'-K'+str(k)+'-I'+str(inst))
    elif probName == 'NKQ':
        q = int(tl.getArgv())
        model = nkq.NKQLandcape(n, k, q, prefixNKQ+'NKQ-N'+str(n)+'-K'+str(k)+'-I'+str(inst)+'-Q'+str(q))

    if compMeth == 'wal':
        start = time.time()
        w = model.WalshCofLinearLinklist()
        walTime = time.time() - start

        """ store runtime to files """
        if probName == 'NKQ':
            nameOfF = waltimeDir+probName+'-'+algoName+'-F'+fitName+'-C'+compMeth+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'-K'+str(k)+'-Q'+str(q)+'.txt'
        elif probName == 'NK':
            nameOfF = waltimeDir+probName+'-'+algoName+'-F'+fitName+'-C'+compMeth+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'-K'+str(k)+'.txt'

        f = open(nameOfF, 'w')
        print >>f,"%g" % (walTime)
        f.close()

#    bit,fit = tl.compFit(model)
#    for i in zip(bit,fit):
#        print i
#    print 'bit',bit
#    print 'fit',fit
#    print 'mean',np.mean(fit)
#    print 'w', w

    if algoName.find('LS') != -1:
        algo = ls.LocalSearch(model, maxFit, n)
    elif algoName.find('GA') != -1:
        algo = ga.GeneticAlgorithm( model.compFit, maxFit, popSize, n )
    elif algoName.find('CHC') != -1:
        algo = chc.CHC()

    tAll = np.zeros(runs)
    for i in range(runs):
        start = time.time()
        if algoName.find('GA') != -1:
            res.append(algo.run(crossoverR, mutationR, fitName))
        elif algoName == 'LS':
            res.append(algo.run(fitName, minimize = True, restart = False,compM = compMeth ))
        elif algoName == 'rLS':
            res.append(algo.run(fitName, minimize = True, restart = True,compM = compMeth))
        elif algoName.find('CHC') != -1:
            res.append(algo.run(model.compFit, maxFit, popSize, n, D, DR, M, fitName))
        tAll[i] = time.time() - start

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

    f = open(nameOfF, 'w')
    for i in range(len(res)):
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
    for i in range(len(tAll)):
        print >>f,"%g" % (tAll[i])
    f.close()
