import nkLandscape as nk
import nkqLandscape as nkq
import WalshAnalysis as wal
import geneticAlgorithm as ga
import AutoCorrelation as ac
#import matplotlib.pyplot as plt
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
import os
import sys

""" consider as a minimization problem """
tl.checkParam(sys.argv)

rseed = 1
nameOfDir = 'result'
prefixNK = './benchmark/NK/'
prefixNKQ = './benchmark/NKQ/'

random.seed(rseed)

""" create files for storing results """
if os.path.isdir(nameOfDir) == False:
    os.mkdir(nameOfDir)

probName = sys.argv[1]
algoName = sys.argv[2]
fitName = sys.argv[3] # fit/mean/std
inst = int(sys.argv[4])
s = sys.argv[5]
n = int(sys.argv[6])
if probName != 'SAT':
    k = int(sys.argv[7])


maxFit = 1000 * n
runs = 30
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

print 'probName', probName, 'algoName', algoName, 'fitName', fitName,  'instanceNo', inst, 'n', n, 'maxFit', maxFit

if probName == 'SAT':
    """ need to perform multiple runs for each instance """
    """ with SAT, we are forced to set n to 100 """
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

#    for i in range(runs):
#        if algoName.find('GA') != -1:
#            res.append(algo.run(crossoverR, mutationR, fitName, minimize = False))
#        elif algoName.find('LS') != -1:
#            res.append(algo.run(fitName, minimize = False))
#        elif algoName.find('CHC') != -1:
#            res.append(algo.run(model.compFit, maxFit, popSize, n, D, DR, M, fitName, minimize = False))

    if probName == 'SAT':
        nameOfF = './result/'+probName+'-'+algoName+'-F'+fitName+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'.txt'
    f = open(nameOfF, 'w')
    for i in range(len(res)):
        if fitName != 'fit':
            print >>f,"%g\t%g\t%g" % (res[i]['sol'], res[i]['fitG'], res[i]['nEvals'])
        else:
            print >>f,"%g\t%g" % (res[i]['sol'], res[i]['nEvals'])
    f.close()

else:
    res = []

    if probName == 'NK':
        model = nk.NKLandscape(n,k,prefixNK+'NK-N'+str(n)+'-K'+str(k)+'-I'+str(inst))
    elif probName == 'NKQ':
        q = int(sys.argv[8])
        model = nkq.NKQLandcape(n, k, q, prefixNKQ+'NKQ-N'+str(n)+'-K'+str(k)+'-I'+str(inst)+'-Q'+str(q))

    if algoName.find('LS') != -1:
        algo = ls.LocalSearch(model.compFit, maxFit, n)
    elif algoName.find('GA') != -1:
        algo = ga.GeneticAlgorithm( model.compFit, maxFit, popSize, n )
    elif algoName.find('CHC') != -1:
        algo = chc.CHC()

#    for i in range(runs):
#        if algoName.find('GA') != -1:
#            res.append(algo.run(crossoverR, mutationR, fitName))
#        elif algoName.find('LS') != -1:
#            res.append(algo.run(fitName))
#        elif algoName.find('CHC') != -1:
#            res.append(algo.run(model.compFit, maxFit, popSize, n, D, DR, M, fitName))

    """ store to files """
    if probName == 'NKQ':
        nameOfF = './result/'+probName+'-'+algoName+'-F'+fitName+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'-K'+str(k)+'-Q'+str(q)+'.txt'
    elif probName == 'NK':
        nameOfF = './result/'+probName+'-'+algoName+'-F'+fitName+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'-K'+str(k)+'.txt'

    f = open(nameOfF, 'w')
    for i in range(len(res)):
        if fitName != 'fit':
            print >>f,"%g\t%g\t%g" % (res[i]['sol'], res[i]['fitG'], res[i]['nEvals'])
        else:
            print >>f,"%g\t%g" % (res[i]['sol'], res[i]['nEvals'])
    f.close()
