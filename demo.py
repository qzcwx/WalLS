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
numOfInstance = int(sys.argv[3])
n = int(sys.argv[4])
k = int(sys.argv[5])


maxFit = 1000 * n
runs = 30
popSize = 50 # always keep popSize to even number

#maxFit = 20
#runs = 1
#popSize = 4

crossoverR = 0.8 # typically in (0.6, 0.9)
mutationR = 1.0/float(n) # typically between 1/popSize and 1/dim
# for CHC 
D = n/4.0
DR = 0.35
M = 1


print 'probName', probName, 'algoName', algoName, 'num of instances', numOfInstance, 'n', n, 'k', k, 'maxFit', maxFit

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
            if algoName.find('SATGA') != -1:
                res.append(algo.runNeigh(crossoverR, mutationR))
            elif algoName.find('GA') != -1:
                res.append(algo.run(crossoverR, mutationR))
            elif algoName.find('SATLS') != -1:
                res.append(algo.runNeigh())
            elif algoName.find('LS') != -1:
                res.append(algo.run())
            elif algoName.find('SATCHC') != -1:
                res.append(algo.runNeigh(model.compFit, maxFit, popSize, n, D, DR, M))
            elif algoName.find('CHC') != -1:
                res.append(algo.run(model.compFit, maxFit, popSize, n, D, DR, M))

else:
    res = []
    for inst in range(numOfInstance):

        if probName == 'NK':
            model = nk.NKLandscape(n,k)
        elif probName == 'NKQ':
            q = int(sys.argv[6])
            model = nkq.NKQLandcape(n, k, q)

        if algoName.find('LS') != -1:
            algo = ls.LocalSearch(model.compFit, maxFit, n)
        elif algoName.find('GA') != -1:
            algo = ga.GeneticAlgorithm( model.compFit, maxFit, popSize, n )
        elif algoName.find('CHC') != -1:
            algo = chc.CHC()

        for i in range(runs):
            if algoName.find('SATGA') != -1:
                res.append(algo.runNeigh(crossoverR, mutationR))
            elif algoName.find('GA') != -1:
                res.append(algo.run(crossoverR, mutationR))
            elif algoName.find('SATLS') != -1:
                res.append(algo.runNeigh())
            elif algoName.find('LS') != -1:
                res.append(algo.run())
            elif algoName.find('SATCHC') != -1:
                res.append(algo.runNeigh(model.compFit, maxFit, popSize, n, D, DR, M))
            elif algoName.find('CHC') != -1:
                res.append(algo.run(model.compFit, maxFit, popSize, n, D, DR, M))

""" store to files """
if probName == 'NKQ':
    nameOfF = './result/'+probName+'-'+algoName+'-N'+str(n)+'-K'+str(k)+'-Q'+str(q)+'.txt'
if probName == 'SAT':
    nameOfF = './result/'+probName+'-'+algoName+'-N'+str(n)+'.txt'
else:
    nameOfF = './result/'+probName+'-'+algoName+'-N'+str(n)+'-K'+str(k)+'.txt'

f = open(nameOfF, 'w')
for i in range(len(res)):
    if algoName.find('SAT') != -1:
        print >>f,"%g\t%g\t%g" % (res[i]['sol'], res[i]['fitG'], res[i]['nEvals'])
    else:
        print >>f,"%g\t%g" % (res[i]['sol'], res[i]['nEvals'])
f.close()

