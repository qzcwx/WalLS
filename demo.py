import nkLandscape as nk
import geneticAlgorithm as ga
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

""" consider as a minimization problem """
n = 15
k = 10
runs = 25
maxFit = 5000
popSize = 50
rseed = 0

random.seed(rseed)
model = nk.NKLandscape(n,k)
print 'Global Optimum', globalOpt(model)
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
