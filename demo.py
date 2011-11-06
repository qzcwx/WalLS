import nkLandscape as nk
import nkqLandscape as nkq
import WalshAnalysis as wal
import geneticAlgorithm as ga
import AutoCorrelation as ac
import matplotlib.pyplot as plt
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

rseed = 3
nameOfDir = 'result'

random.seed(rseed)

""" create files for storing results """
if os.path.isdir(nameOfDir) == False:
    os.mkdir(nameOfDir)

probName = sys.argv[1]
algoName = sys.argv[2]
n = int(sys.argv[3])
k = int(sys.argv[4])

if probName == 'NK':
    model = nk.NKLandscape(n,k)
elif probName == 'SAT':
    model = mx.MAXSAT()
elif probName == 'NKQ':
    q = int(sys.argv[5])
    model = nkq.nkqlandcape(n, k, q)

maxFit = 1000 * n
runs = 1
popSize = 3 # always keep popSize to even number
crossoverR = 0.8 # typically in (0.6, 0.9)
mutationR = 1.0/float(n) # typically between 1/popSize and 1/dim
#maxFit = 100000

print 'probName', probName, 'algoName', algoName,  'n', n, 'k', k, 'maxFit', maxFit

if probName == 'SAT':
    """ need to perform multiple runs for each instance """
    res = []
    noInstance = 1
    model.setInstance(noInstance)
    print 'No. Instance', noInstance

    if algoName.find('LS') != -1:
        algo = ls.LocalSearch(model.compFit, maxFit, n)
    elif algoName.find('GA') != -1:
        algo = ga.GeneticAlgorithm( model.compFit, maxFit, popSize, n )

    for i in range(runs):
        if algoName.find('SATGA') != -1:
            res.append(algo.runNeigh(crossoverR, mutationR))
        elif algoName.find('GA') != -1:
            res.append(algo.run(crossoverR, mutationR))
        elif algoName.find('SATLS') != -1:
            res.append(algo.runNeigh())
        elif algoName.find('LS') != -1:
            res.append(algo.run())

else:
    if algoName.find('LS') != -1:
        algo = ls.LocalSearch(model.compFit, maxFit, n)
    elif algoName.find('GA') != -1:
        algo = ga.GeneticAlgorithm( model.compFit, maxFit, popSize, n )

    res = []
    for i in range(runs):
        if algoName.find('SATGA') != -1:
            res.append(algo.runNeigh(crossoverR, mutationR))
        elif algoName.find('GA') != -1:
            res.append(algo.run(crossoverR, mutationR))
        elif algoName.find('SATLS') != -1:
            res.append(algo.runNeigh())
        elif algoName.find('LS') != -1:
            res.append(algo.run())
        pdb.set_trace()

""" store to files """
if probName == 'NKQ':
    nameOfF = './result/'+probName+'-'+algoName+'-N'+str(n)+'-K'+str(k)+'-Q'+str(q)+'.txt'
else:
    nameOfF = './result/'+probName+'-'+algoName+'-N'+str(n)+'-K'+str(k)+'.txt'

f = open(nameOfF, 'w')
for i in range(len(res)):
    if algoName.find('SAT') != -1:
        print >>f,"%g\t%g\t%g" % (res[i]['sol'], res[i]['fitG'], res[i]['nEvals'])
    else:
        print >>f,"%g\t%g" % (res[i]['sol'], res[i]['nEvals'])
f.close()

