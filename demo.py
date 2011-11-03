import nkLandscape as nk
import WalshAnalysis as wal
import geneticAlgorithm as ga
import AutoCorrelation as ac
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


algoName = sys.argv[1]
n = int(sys.argv[2])
k = int(sys.argv[3])

print 'n', n, 'k', k, 'runs', runs, 'seed', rseed

start = time.time()
model = nk.NKLandscape(n,k)
print 'model building t0\t\t', time.time() - start

if algoName.find('LS') != -1:
    algo = ls.LocalSearch(model.compFit, maxFit, n)
elif algoName.find('GA') != -1:
    algo = ga.GeneticAlgorithm( model.compFit, maxFit, popSize, n )

start = time.time()
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
nameOfF = './result/'+algoName+'-N'+str(n)+'-K'+str(k)+'.txt'
f = open(nameOfF, 'w')

for i in range(len(res)):
    if algoName.find('SAT') != -1:
        print >>f,"%g\t%g\t%g" % (res[i]['sol'], res[i]['fitG'], res[i]['nEvals'])
    else:
        print >>f,"%g\t%g" % (res[i]['sol'], res[i]['nEvals'])

f.close()

#start = time.time()
#w = model.WalCofLinear()
#print 'Walsh t1\t\t\t', time.time() - start

#start = time.time()
#w0 = model.WalshCofLinearLinklist()
##print 'w\n', w0
#print 'Walsh t10\t\t\t', time.time() - start

#start = time.time()
#autoCo =  ac.autoCorr(s, w) 
##autoCo = [ ac.autoCorr(i,w) for i in range(0,20) ]
#print 'autoCo\n', autoCo
#print "Auto Correlation t2\t\t", time.time() - start

#start = time.time()
##autoCo0 = ac.autoCorrLinerSpace(s, w0)
#autoCo0 = [ ac.autoCorrLinerSpace(i,w0) for i in range(0,20) ]
#print 'autoCoLinear\n', autoCo0
#print "Auto Correlation t20\t\t", time.time() - start

#start = time.time()
#bitStr, fit = compFit(model)
#print "Evaluate all candidate t3\t", time.time() - start
#print 'Global Optimum', min(fit)
#print 'bitStr', bitStr
#print 'fit', fit
#
#start = time.time()
#numOpt = lo.localOpt(bitStr, fit)
#print 'Local Optimum t4\t\t', time.time() - start
#
#print 'Num of LocalOptimum', numOpt

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
