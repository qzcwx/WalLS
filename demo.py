import nkLandscape as nk
import WalshAnalysis as wal
import geneticAlgorithm as ga
import AutoCorrelation as ac
import LocalOptima as lo
import numpy as np
import matplotlib.pyplot as plt
import random
import math
import time
import pdb

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
n = 100
k = 1
runs = 25
maxFit = 1000
popSize = 50 # always keep popSize to even number
crossoverR = 0.8 # typically in (0.6, 0.9)
mutationR = 1.0/n # typically between 1/popSize and 1/dim
rseed = 0
s = 1

random.seed(rseed)

sols1 = np.zeros(runs)
sols2 = np.zeros(runs)

print 'n', n, 'k', k, 'runs', runs, 'maxFit', maxFit, 'popSize', popSize, 'crossoverR', crossoverR, 'mutationR', mutationR, 'seed', rseed
start = time.time()
model = nk.NKLandscape(n,k)
#print 'Global Optimum', globalOpt(model)
#algo = ga.GeneticAlgorithm( model.compFit, maxFit, popSize, n )
print 'model building t0\t\t', time.time() - start

start = time.time()
w = model.WalCofLinear()
print 'Walsh t1\t\t\t', time.time() - start

start = time.time()
autoCo =  ac.autoCorr(s, w) 
#autoCo = [ ac.autoCorr(i,w) for i in range(0,20) ]
print "Auto Correlation t2\t\t", time.time() - start

#start = time.time()
#bitStr, fit = compFit(model)
#print "Evaluate all candidate t3\t", time.time() - start
#
#start = time.time()
#numOpt = lo.localOpt(bitStr, fit)
#print 'Local Optimum t4\t\t', time.time() - start

print 'AutoCorrelation', autoCo
print 'w', w
# print 'Num of LocalOptimum', numOpt

#for k in range(0,20):
#    for n in [i for i in range(2, 20) if i >= k+1 ]:
#        mutationR = 1.0/n # typically between 1/popSize and 1/dim
#        print 
#        print 'n', n, 'k', k, 'runs', runs, 'maxFit', maxFit, 'popSize', popSize, 'crossoverR', crossoverR, 'mutationR', mutationR, 'seed', rseed
#        start = time.time()
#        model = nk.NKLandscape(n,k)
#        #print 'Global Optimum', globalOpt(model)
#        algo = ga.GeneticAlgorithm( model.compFit, maxFit, popSize, n )
#        print 'model building t0\t\t', time.time() - start
#
#
#        start = time.time()
#        w = model.WalCofLinear()
#        print 'Walsh t1\t\t\t', time.time() - start
#
#        start = time.time()
#        autoCo =  ac.autoCorr(s, w) 
#        print "Auto Correlation t2\t\t", time.time() - start
#
#        start = time.time()
#        bitStr, fit = compFit(model)
#        print "Evaluate all candidate t3\t", time.time() - start
#
#        start = time.time()
#        numOpt = lo.localOpt(bitStr, fit)
#        print 'Local Optimum t4\t\t', time.time() - start
#
#        start = time.time()
#        for i in range(runs):
#            sols1[i] = algo.run(crossoverR, mutationR)
#            sols2[i] = algo.runNeigh(crossoverR, mutationR)
#        print 'Run GAs t5\t\t\t', time.time() - start
#        print 'AutoCorrelation', autoCo, 'Num of LocalOptimum', numOpt, 'Global Optimum', min(fit)
#        print '1. median', np.median(sols1), 'mean', np.mean(sols1), 'worst', max(sols1), 'best', min(sols1), 'std', np.std(sols1)
#        print '2. median', np.median(sols2), 'mean', np.mean(sols2), 'worst', max(sols2), 'best', min(sols2), 'std', np.std(sols2)
#        print 
