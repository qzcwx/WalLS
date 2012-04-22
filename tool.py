import LocalSearch as ls
#import matplotlib.pyplot as plt
import random
import numpy as np
import math
import sys

argvCount = 1

def getArgv():
    global argvCount
    arg = sys.argv[argvCount]
    argvCount = argvCount + 1
    return arg

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
                bestBit = bit
    return bestFit, bestBit

def compFit(model):
    n = model.getN()
    fit = np.zeros(math.pow(2,n))
    bitStr = genSeqBits(n)
    for i in range(int(math.pow(2,n))):
       fit[i] = model.compFit(bitStr[i])
    return bitStr, fit

def plotDistNKQ():
    """ 
        Plot the distribution of fitness of when K and Q vary, 1000 samples,
        generate one instance of NK-Q landscapes, for:
            * N = 20, 50, 100
            * K = 0, 2, 4, 8, 16
            * q = 2, 4, 8, 16
    """
    nSamples = 1000

    inst = 0
    prefixNKQ = './benchmark/NKQ/'

    for n in [20, 50, 100]:
        for k in [0, 2, 4, 8, 16]:
            for q in [2, 4, 8, 16]:
                print 'N=',n,'K=',k,'Q=',q
                model = nkq.NKQLandcape(n, k, q, prefixNKQ+'NKQ-N'+str(n)+'-K'+str(k)+'-I'+str(inst)+'-Q'+str(q))
                res = np.zeros((nSamples, 3))
                for r in range(nSamples):
                    randBitStr = []
                    for j in range(n):
                        if random.random()<0.5:
                            randBitStr.append('0')
                        else:
                            randBitStr.append('1')
                    res[r][0] = evalSol(randBitStr,model,'fit',True)
                    res[r][1] = evalSol(randBitStr,model,'mean',True)
                    res[r][2] = evalSol(randBitStr,model,'std',True)
                    #res[r] = model.compFit(randBitStr)

                        
                plt.figure()
                plt.hist(res, histtype='bar',
                            label=['fit', 'mean', 'std'])
                plt.title('N='+str(n)+',K='+str(k)+',Q='+str(q)+',I='+str(inst))
                plt.legend()
                plt.savefig('N='+str(n)+'K='+str(k)+'Q='+str(q)+'I='+str(inst))

#def plotDistNK():
#    """ 
#        Plot the distribution of fitness of when K and Q vary, 1000 samples,
#        generate one instance of NK landscapes, for:
#            * N = 20, 50, 100
#            * K = 0, 2, 4, 8, 16
#            * q = 2, 4, 8, 16
#    """
#    nSamples = 1000
#
#    inst = 0
#    prefixNK = './benchmark/NK/'
#
#    for n in [20, 50, 100]:
#        for k in [0, 2, 4, 8, 16]:
#            print 'N=',n,'K=',k
#            model = nk.NKLandscape(n,k,prefixNK+'NK-N'+str(n)+'-K'+str(k)+'-I'+str(inst))
#            res = np.zeros((nSamples, 3))
#            for r in range(nSamples):
#                randBitStr = []
#                for j in range(n):
#                    if random.random()<0.5:
#                        randBitStr.append('0')
#                    else:
#                        randBitStr.append('1')
#                res[r][0] = evalSol(randBitStr,model,'fit',True)
#                res[r][1] = evalSol(randBitStr,model,'mean',True)
#                res[r][2] = evalSol(randBitStr,model,'std',True)
#                #res[r] = model.compFit(randBitStr)
#
#                    
#            plt.figure()
#            plt.hist(res, histtype='bar',
#                        label=['fit', 'mean', 'std'])
#            plt.title('N='+str(n)+',K='+str(k)+',I='+str(inst))
#            plt.legend()
#            plt.savefig('N='+str(n)+'K='+str(k)+'I='+str(inst))


def checkParam(argv):
    if len(argv) == 1:
        print 'Usage: python demo.py [ComputeMethod] [NameOfProblem] [NameOfAlgorithm] [fit/mean/std] [overwrite] [I] [PopSize] [N] [K] [Q]'
        print """
Constrains: 
1) for SAT, I = [1,100]
2) for SAT, N = 100, K/Q = None
3) for NK/NKQ, I = [0,9]
4) for NK, Q = None
5) for NKQ, Q = {2, 4, 8, 16}
6) for LS,  PopSize = 1
7) [ComputeMethod] = {bf (brute force), walWalk (walsh coefficient with random walk), walRest (walsh coefficient with random walk), supm (super move), bitImp (bit impact), walSearch (Search completely relying on Walsh terms), null (run nothing)}. 
    If 'walSearch', rLS, fit. 
        """
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
                #model = nkq.NKQLandcape(n, n-1, q)

                model = nkq.NKQLandcape(n, k, q, prefixNKQ+'NKQ-N'+str(n)+'-K'+str(k)+'-I'+str(inst)+'-Q'+str(q))
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

def evalSol(s,model,fitName,minimize):
    if fitName == 'fit':
        return model.compFit(s)
    else :
        fitN = np.zeros(len(s))
        for j in range(len(s)):
            # flip the jth bit in bit-string
            neighStr = np.copy(s)
            if neighStr[j] == '1':
                neighStr[j] = '0'
            else:
                neighStr[j] = '1'
            fitN[j] = model.compFit(neighStr)
        if fitName == 'mean':
            return np.mean(fitN)
        elif minimize == True :
            return np.mean(fitN) - np.std(fitN)
        elif minimize == False :
            return np.mean(fitN) + np.std(fitN)

def plotFitRank():
    nSamples = 1000

    inst = 0
    prefixNKQ = './benchmark/NKQ/'

#    for n in [20, 50, 100]:
    for n in [100]:
        for k in [16]:
            for q in [2, 4, 8, 16]:
                print 'N=',n,'K=',k,'Q=',q
                maxFit = 1000 * n
                model = nkq.NKQLandcape(n, k, q, prefixNKQ+'NKQ-N'+str(n)+'-K'+str(k)+'-I'+str(inst)+'-Q'+str(q))
                algo = ls.LocalSearch(model.compFit, maxFit, n)
                res = algo.run('fit')
                fSort = fitRank(res['bit'], model.compFit)

                plt.figure()
                plt.stem(range(1,n+1), fSort)
                plt.title('N='+str(n)+' K='+str(k)+' Q='+str(q))
                plt.savefig('N='+str(n)+' K='+str(k)+' Q='+str(q))

def fitRank(s, fit):

    fitN = np.zeros(len(s))
    for j in range(len(s)):
        # flip the jth bit in bit-string
        neighStr = np.copy(s)
        if neighStr[j] == '1':
            neighStr[j] = '0'
        else:
            neighStr[j] = '1'
        fitN[j] = fit(neighStr)
    fitNSort = sorted(fitN)

    return fitNSort


def listMerge(l1, l2):
    """  
    merge two sorted list
    """

    if not l1:
        return l2
    elif not l2:
        return l1
    else:
        merge = []
        i = 0 
        j = 0 
        while True:
            if l1[i] < l2[j]:
                if not merge or merge[-1] != l1[i]:
                    merge.append(l1[i])
                i = i+1
            else:
                if not merge or merge[-1] != l2[j]:
                    merge.append(l2[j])
                j = j+1
            if i>=len(l1) :
                while j<len(l2):
                    if merge[-1] != l2[j]:
                        merge.append(l2[j])
                    j = j+1
                break
            elif j>=len(l2):
                while i<len(l1):
                    if merge[-1] != l1[i]:    
                        merge.append(l1[i])
                    i = i+1
                break
        return merge
    
#plotDistNKQ()
#plotDistNK()
#plotDistMaxK()
#plotFitRank()

def genSeqBits(n):
    bitStr = []
    for i in range(int(math.pow(2,n))):
       bit = bin(i)
       bit = bit[2:]
       if len(bit) < n:
           bit = (n - len(bit))*'0' + bit
       bitStr.append(bit)
    return bitStr

def printSep():
    print 80*'*'
