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
n = 15
k = 10
maxFit = 1000
popSize = 20

rseed = 0

random.seed(rseed)
model = nk.NKLandscape(n,k)
bitStr, fit = compFit(model)
print 'Global Optimum', min(fit)
algo = ga.GeneticAlgorithm( model.compFit, maxFit, popSize, n )
algo.runNeigh()
algo.run()
