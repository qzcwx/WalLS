import nkLandscape as nk
import nkqLandscape as nkq
#import matplotlib.pyplot as plt
import random
import numpy as np
import math
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

def plotDist():
    """ 
        Plot the distribution of fitness of when K and Q vary, 1000 samples,
        generate one instance of NK-Q landscapes, for:
            * N = 20, 50, 100
            * K = 0, 2, 4, 8, 16
            * q = 2, 4, 8, 16
    """
    nSamples = 1000
    for n in [20, 50, 100]:
        for k in [0, 2, 4, 8, 16]:
            for q in [2, 4, 8, 16]:
                print 'N=',n,'K=',k,'Q=',q
                model = nkq.NKQLandcape(n, k, q)
                res = np.zeros(nSamples)
                for r in range(nSamples):
                    randBitStr = []
                    for j in range(n):
                        if random.random()<0.5:
                            randBitStr.append('0')
                        else:
                            randBitStr.append('1')
                    res[r] = model.compFit(randBitStr)
                plt.figure()
                plt.hist(res)
                plt.title('N='+str(n)+', K='+str(k)+',Q='+str(q))
                plt.savefig('N='+str(n)+', K='+str(k)+',Q='+str(q))

def checkParam(argv):
    if len(argv) == 1:
        print 'Usage: python demo.py [NameOfProblem] [NameOfAlgorithm] [fit/mean/std] [I] [N] [K] [Q]'
        #print 'Usage: python demo.py [NameOfProblem] [NameOfAlgorithm] [I] [N] [K] [Q]'
        sys.exit()
 
def plotDistMaxK():
    """ 
        Plot the distribution of fitness of when K and Q vary, 1000 samples,
        generate one instance of NK-Q landscapes, for:
            * N = 20, 50, 100
            * K = 0, 2, 4, 8, 16
            * q = 2, 4, 8, 16
    """
    nSamples = 1000
    for n in [5, 10, 15]:
            k = n-1
            for q in [2, 4, 8, 16]:
                print 'N=',n,'K=',k,'Q=',q
                model = nkq.NKQLandcape(n, n-1, q)
                res = np.zeros(nSamples)
                for r in range(nSamples):
                    randBitStr = []
                    for j in range(n):
                        if random.random()<0.5:
                            randBitStr.append('0')
                        else:
                            randBitStr.append('1')
                    res[r] = model.compFit(randBitStr)
                plt.figure()
                plt.hist(res)
                plt.title('N='+str(n)+' K='+str(k)+' Q='+str(q))
                plt.savefig('N='+str(n)+' K='+str(k)+' Q='+str(q))
#plotDist()
#plotDistMaxK()
