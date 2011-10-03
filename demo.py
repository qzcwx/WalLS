import nkLandscape as nk
import geneticAlgorithm as ga
import numpy as np
import matplotlib.pyplot as plt
import random
import math
import time
import pdb

def compFit(model):
    n = model.getN()
    fit = np.zeros(math.pow(2,n))
    bitStr = nk.genSeqBits(n)
    for i in range(int(math.pow(2,n))):
       fit[i] = model.compFit(bitStr[i])
    return bitStr, fit

""" consider as a minimization problem """
n = 20
k = 15
runs = 25
maxFit = 3000
popSize = 30
rseed = 0

random.seed(rseed)
model = nk.NKLandscape(n,k)
#bitStr, fit = compFit(model)
#print 'Global Optimum', min(fit)
algo = ga.GeneticAlgorithm( model.compFit, maxFit, popSize, n )

sols1 = np.zeros(runs)
sols2 = np.zeros(runs)

for i in range(runs):
    print '*********************************************'
    print 'Run', i
    sols1[i] = algo.run()
    print 'solution 1\n', sols1[i] 
    sols2[i] = algo.runNeigh()
    print 'solution 2\n', sols2[i]
    print '*********************************************'
    print 

print '1. mean', np.mean(sols1), 'worst', max(sols1), 'best', min(sols1)
print '2. mean', np.mean(sols2), 'worst', max(sols2), 'best', min(sols2)
