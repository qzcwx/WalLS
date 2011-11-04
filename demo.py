import nkLandscape as nk
import nkqLandscape as nkq
import WalshAnalysis as wal
import geneticAlgorithm as ga
import AutoCorrelation as ac
import MAXSAT as mx
import LocalSearch as ls
import LocalOptima as lo
import numpy as np
import random
import math
import time
import pdb
import os
import sys

def globalOpt(model):
    """ find the global optimum on the fly """
    n = model.getN()
    for i in range(int(math.pow(2,n))):
        bit = bin(i)
        bit = bit[2:]
        if len(bit) < n:
            bit = (n - len(bit))*'0' + bit
        if i == 0:
            bestFit = model.compFit(bit)
        else:
            curFit = model.compFit(bit)
            if curFit < bestFit:
                bestFit = curFit
    return bestFit

def compFit(model):
    n = model.getN()
    fit = np.zeros(math.pow(2,n))
    bitStr = nk.genSeqBits(n)
    for i in range(int(math.pow(2,n))):
       fit[i] = model.compFit(bitStr[i])
    return bitStr, fit

""" consider as a minimization problem """
n = 50
k = 20
runs = 25
maxFit = 10000
popSize = 50 # always keep popSize to even number
crossoverR = 0.8 # typically in (0.6, 0.9)
mutationR = 1.0/float(n) # typically between 1/popSize and 1/dim
rseed = 3
s = 0
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
    model = nkq.NKQLandcape(n, k, q)

print 'probName', probName, 'algoName', algoName,  'n', n, 'k', k

if probName == 'SAT':
    """ need to perform multiple runs for each instance """
    res = []
    noInstance = 1
    #for noInstance in range(1,1001):
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
            print algo.run()
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
