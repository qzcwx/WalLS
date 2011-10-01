import nkLandscape as nk
import geneticAlgorithm as ga
import numpy as np
import matplotlib.pyplot as plt
import random
import math
import time
import pdb

""" consider as a minimization problem """
n = 30
k = 0
maxFit = 10000
popSize = 30

rseed = 0

random.seed(rseed)
model = nk.NKLandscape(n,k)
algo = ga.GeneticAlgorithm( model.compFit, maxFit, popSize, n )
algo.run()
