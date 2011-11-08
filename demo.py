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

rseed = 0
nameOfDir = 'result'

random.seed(rseed)

""" create files for storing results """
if os.path.isdir(nameOfDir) == False:
    os.mkdir(nameOfDir)

probName = sys.argv[1]
algoName = sys.argv[2]
fitName = sys.argv[3] # fit/mean/std
numOfInstance = int(sys.argv[4])
s = sys.argv[5]
n = int(sys.argv[6])
if probName != 'SAT':
    k = int(sys.argv[7])


maxFit = 100 * n
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

print 'probName', probName, 'algoName', algoName, 'fitName', fitName,  'num of instances', numOfInstance, 'n', n, 'maxFit', maxFit

if probName == 'SAT':
    """ need to perform multiple runs for each instance """
    """ with SAT, we are forced to set n to 100 """
    model = mx.MAXSAT()
    res = []
    instances = random.sample(range(1,1001), numOfInstance)
    for inst in instances:
        model.setInstance(inst)
        print 'Instance', inst

        if algoName.find('LS') != -1:
            algo = ls.LocalSearch(model.compFit, maxFit, n)
        elif algoName.find('GA') != -1:
            algo = ga.GeneticAlgorithm( model.compFit, maxFit, popSize, n )
        elif algoName.find('CHC') != -1:
            algo = chc.CHC()

        for i in range(runs):
            if algoName.find('GA') != -1:
                res.append(algo.run(crossoverR, mutationR, fitName))
            elif algoName.find('LS') != -1:
                res.append(algo.run(fitName))
            elif algoName.find('CHC') != -1:
                res.append(algo.run(model.compFit, maxFit, popSize, n, D, DR, M, fitName))

else:
    res = []
    for inst in range(numOfInstance):

        if probName == 'NK':
            model = nk.NKLandscape(n,k)
        elif probName == 'NKQ':
            q = int(sys.argv[8])
            model = nkq.NKQLandcape(n, k, q)

        if algoName.find('LS') != -1:
            algo = ls.LocalSearch(model.compFit, maxFit, n)
        elif algoName.find('GA') != -1:
            algo = ga.GeneticAlgorithm( model.compFit, maxFit, popSize, n )
        elif algoName.find('CHC') != -1:
            algo = chc.CHC()

        for i in range(runs):
            if algoName.find('GA') != -1:
                res.append(algo.run(crossoverR, mutationR, fitName))
            elif algoName.find('LS') != -1:
                res.append(algo.run(fitName))
            elif algoName.find('CHC') != -1:
                res.append(algo.run(model.compFit, maxFit, popSize, n, D, DR, M, fitName))

""" store to files """
if probName == 'NKQ':
    nameOfF = './result/'+probName+'-'+algoName+'-F'+fitName+'-N'+str(n)+'-K'+str(k)+'-Q'+str(q)+'.txt'
elif probName == 'SAT':
    nameOfF = './result/'+probName+'-'+algoName+'-F'+fitName+'-N'+str(n)+'.txt'
elif probName == 'NK':
    nameOfF = './result/'+probName+'-'+algoName+'-F'+fitName+'-N'+str(n)+'-K'+str(k)+'.txt'

f = open(nameOfF, 'w')
for i in range(len(res)):
    if algoName.find('SAT') != -1:
        print >>f,"%g\t%g\t%g" % (res[i]['sol'], res[i]['fitG'], res[i]['nEvals'])
    else:
        print >>f,"%g\t%g" % (res[i]['sol'], res[i]['nEvals'])
f.close()

