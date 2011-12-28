""" 
Local Search module:
    This local search module gets the problem (in terms of objective function),
    and return a solution found by hill climber
"""

import numpy as np
import WalshAnalysis as wal
import random
import math
import copy
import sys
import time
import pdb

class Struct:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

class LocalSearch:
    def __init__(self, model, MaxFit, dim):
        self.func = model.compFit
        self.model = model
        self.MaxFit = MaxFit
        self.dim = dim

    def initIndiv(self, dim):
        """ initial the search inidividual with random bit string """
        indiv = Struct(fit = 0, bit = '0') 
        randBitStr = []
        for j in range(dim):
            if random.random()<0.5:
                randBitStr.append('0')
            else:
                randBitStr.append('1')
        indiv = Struct( fit = 0, bit = randBitStr)
        return indiv

    def initIndivNeigh(self, dim):
        """ initial the search inidividual with random bit string """
        indiv = Struct(fit = 0, fitG = 0, bit = '0') 
        randBitStr = []
        for j in range(dim):
            if random.random()<0.5:
                randBitStr.append('0')
            else:
                randBitStr.append('1')
        indiv = Struct( fit = 0, fitG = 0, bit = randBitStr)
        return indiv

    def run(self, fitName, minimize, restart, compM = 'wal'):
            if compM != 'wal':
                if fitName =='fit': 
                    return self.runFit(minimize,restart)
                else:
                    return self.runNeigh(fitName, minimize,restart)
            else :
                if fitName == 'fit':
                    return self.runFitWal(fitName, minimize, restart)
                elif fitName == 'mean':
                    return self.runMeanWal(fitName, minimize, restart)

    def runFitWal(self,fitName, minimize, restart):
        """ 
        Walsh Analysis for speeding up steepest descent local search
        """
        self.oldindiv = self.initIndiv(self.dim)
        self.fitEval = 0
        
        self.transWal()

        self.indiv = copy.deepcopy(self.oldindiv)

        self.initFitWal()

        self.oldindiv = self.evalPop(self.oldindiv)
        self.bsf = copy.deepcopy(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)

#        self.initC()
        self.WA = []
#        print 'C', self.C
#        self.trace = [Struct(fitEval= self.fitEval,fit = self.oldindiv.fit, fitG = self.oldindiv.fitG)]
        compPSumT = 0
        updateT = 0
        while self.fitEval < self.MaxFit:
            # generate neighborhood and compute their fitness
            neighPop = self.neighWal()
            # compute the fitG (mean) of each neighborhood individuals
            improveN = False
            nCount = 0
            #print 
            #print 'current', self.oldindiv.bit, 'fit', self.oldindiv.fit, 'fitG', self.oldindiv.fitG
            oldFit = self.oldindiv.fit
            for n in neighPop:
                self.indiv = copy.deepcopy(n)
                self.indiv.fit = oldFit - 2*self.sumArr[nCount]
                self.fitEval = self.fitEval + 1
                #print 'neigh: ', self.indiv.bit, 'fit', self.indiv.fit, 'fitG', self.indiv.fitG
                if self.selectionFit(minimize) == True:
                    #print 'better neigh!'
                    improveN = True
                    changeBit = nCount

                nCount = nCount + 1

#            print 'improveN', improveN
            #print self.fitEval
            #pdb.set_trace()
#            self.trace.append(Struct(fitEval= self.fitEval,fit = self.oldindiv.fit, fitG = self.oldindiv.fitG))
            if improveN == False:
                if restart == True:
                    self.restart(fitName, minimize)
                else:
#                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit,'trace':self.trace}
                    #print 'compPSum', compPSumT
                    #print 'update', updateT
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                self.updateFit(changeBit)

        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit}


    def runMeanWal(self,fitName, minimize, restart):
        """ 
        steepest descent local search with respect to mean of neighs by Walsh Analysis
        """
        self.oldindiv = self.initIndivNeigh(self.dim)
        self.fitEval = 0
        
        start = time.time()
        self.transWal()
        print 'transWal', time.time() - start

        self.indiv = copy.deepcopy(self.oldindiv)

        start = time.time()
        self.initWal()
        print 'initWal', time.time() - start

        self.oldindiv = self.evalPop(self.oldindiv)
        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (sum(self.sumArr))
        self.bsf = copy.deepcopy(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)

#        self.initC()
        self.WA = []
#        print 'C', self.C
#        self.trace = [Struct(fitEval= self.fitEval,fit = self.oldindiv.fit, fitG = self.oldindiv.fitG)]
        compPSumT = 0
        updateT = 0
        while self.fitEval < self.MaxFit:
            # generate neighborhood and compute their fitness
            neighPop = self.neighWal()
            # compute the fitG (mean) of each neighborhood individuals
            improveN = False
            nCount = 0
            #print 
            #print 'current', self.oldindiv.bit, 'fit', self.oldindiv.fit, 'fitG', self.oldindiv.fitG
            oldFit = self.oldindiv.fit
            oldFitG = self.oldindiv.fitG
            for n in neighPop:
                self.indiv = copy.deepcopy(n)

                self.indiv.fit = oldFit - 2*self.sumArr[nCount]
                self.fitEval = self.fitEval + 1

                start = time.time()
                self.indiv.fitG = oldFitG - 2*self.sumArr[nCount] + 4/float(self.dim) * self.compCsum(nCount)
                compPSumT = compPSumT + time.time() - start
                #print 'neigh: ', self.indiv.bit, 'fit', self.indiv.fit, 'fitG', self.indiv.fitG
                if self.selectionFitNeigh(minimize) == True:
                    #print 'better neigh!'
                    improveN = True
                    changeBit = nCount

                nCount = nCount + 1

#            print 'improveN', improveN
            #print self.fitEval
            #pdb.set_trace()
#            self.trace.append(Struct(fitEval= self.fitEval,fit = self.oldindiv.fit, fitG = self.oldindiv.fitG))
            if improveN == False:
                if restart == True:
                    #print 'restart'
                    oldbit = self.oldindiv.bit
                    self.restart(fitName, minimize)
                    newbit = self.oldindiv.bit
                    #print oldbit, newbit
                    diff = self.diffBits(oldbit, newbit)
                    start = time.time()
                    for i in diff:
                        self.update(i)
                    updateT = updateT + time.time() - start
                else:
#                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit,'trace':self.trace}
                    #print 'compPSum', compPSumT
                    #print 'update', updateT
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                start = time.time()
                self.update(changeBit)
                #print 'update'
                #print 'C', self.C
                self.printC()
                updateT = updateT + time.time() - start


#        return { 'nEvals': self.fitEval, 'sol': self.bsf.fit, 'fitG': self.bsf.fitG, 'bit':self.bsf.bit,'trace':self.trace}
        print 'compPSum', compPSumT
        print 'update', updateT
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'fitG': self.bsf.fitG, 'bit':self.bsf.bit}

    def runFit(self, minimize,restart):
        self.oldindiv = self.initIndiv(self.dim)
        self.fitEval = 0
        self.oldindiv = self.evalPop(self.oldindiv)
        self.bsf = copy.deepcopy(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)
#        self.trace = [Struct(fitEval= self.fitEval,fit = self.oldindiv.fit)]
        while self.fitEval < self.MaxFit:
            neighs = self.neighbors()
            improveN = False
            for i in neighs:
                self.indiv.bit = np.copy(i)
                self.indiv = self.evalPop(self.indiv)
                if  self.selectionFit(minimize) == True:
                    improveN = True

#            self.trace.append(Struct(fitEval= self.fitEval,fit = self.oldindiv.fit))
            if improveN == False:
                if restart == True:
                    self.restart('fit', minimize)
                else:
                    #return {'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit,'trace':self.trace}
                    return {'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit}
        #return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit,'trace':self.trace}
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit}


    def runNeigh(self, fitName, minimize,restart):
        self.oldindiv = self.initIndivNeigh(self.dim)
        self.fitEval = 0
        self.oldindiv = self.evalPopNeigh(self.oldindiv, fitName, minimize)
        self.bsf = copy.deepcopy(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)
        while self.fitEval < self.MaxFit:
            neighs = self.neighbors()
            improveN = False
            #print 
            #print 'current', self.oldindiv.bit, 'fit', self.oldindiv.fit, 'fitG', self.oldindiv.fitG
            for i in neighs:
                self.indiv.bit = np.copy(i)
                self.indiv = self.evalPopNeigh(self.indiv, fitName, minimize)
                #print 'neigh: ', self.indiv.bit, 'fit', self.indiv.fit, 'fitG', self.indiv.fitG
                if self.selectionFitNeigh(minimize) == True:
                    improveN = True
                    #print 'better neigh!'
#            print 'improveN', improveN
#            print self.fitEval
            #pdb.set_trace()
            if improveN == False:
                if restart == True:
                    #print 'restart'
                    self.restart(fitName, minimize)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
        return { 'nEvals': self.fitEval, 'sol': self.bsf.fit, 'fitG': self.bsf.fitG, 'bit':self.bsf.bit}

    def diffBits(self, a, b):
        diff = []
        for i in range(self.dim):
            if a[i] != b[i]:
                diff.append(i)
        return diff

    def restart(self, fitName, minimize):

        if fitName == 'fit' and minimize == True :
            if self.bsf.fit > self.oldindiv.fit:
                self.bsf = copy.deepcopy(self.oldindiv)
        elif fitName == 'fit' and minimize == False :
            if self.bsf.fit < self.oldindiv.fit:
                self.bsf = copy.deepcopy(self.oldindiv)
        elif minimize == True :
            if self.bsf.fitG > self.oldindiv.fitG:
                self.bsf = copy.deepcopy(self.oldindiv)
        elif minimize == False :
            if self.bsf.fitG < self.oldindiv.fitG:
                self.bsf = copy.deepcopy(self.oldindiv)

        #print 'before restart'
#        print self.oldindiv.bit, self.oldindiv.fit, self.oldindiv.fitG
#        print self.bsf.bit, self.bsf.fit, self.bsf.fitG

        if fitName == 'fit':
            self.oldindiv = self.initIndiv(self.dim)
            self.oldindiv = self.evalPop(self.oldindiv)
        else :
            self.oldindiv = self.initIndivNeigh(self.dim)
            self.oldindiv = self.evalPopNeigh(self.oldindiv, fitName, minimize)

#        print 'after restart'
#        print self.oldindiv.bit, self.oldindiv.fit, self.oldindiv.fitG
#        print self.bsf.bit, self.bsf.fit, self.bsf.fitG

    def neighbors(self):
        neighs = []
        for j in range(self.dim):
            # flip the jth bit in bit-string
            neighStr = np.copy(self.oldindiv.bit)
            if neighStr[j] == '1':
                neighStr[j] = '0'
            else:
                neighStr[j] = '1'
            neighs.append( neighStr )
        return np.array(neighs)

    def evalPop(self, indiv):
        indiv.fit = self.func(indiv.bit)
        self.fitEval = self.fitEval + 1
        return copy.deepcopy(indiv)

    def evalPopNeigh(self, indiv, fitName, minimize):
        """ evaluate the individual itself """
        indiv.fit = self.func(indiv.bit)
        self.fitEval = self.fitEval + 1
        """ evaluate all neighborhood """
        fitN = np.zeros(self.dim)
        for j in range(self.dim):
            # flip the jth bit in bit-string
            neighStr = np.copy(indiv.bit)
            if neighStr[j] == '1':
                neighStr[j] = '0'
            else:
                neighStr[j] = '1'
            fitN[j] = self.func(neighStr)
        if fitName == 'mean':
            indiv.fitG = np.mean(fitN)
        elif minimize == True :
            indiv.fitG = np.mean(fitN) - np.std(fitN)
        elif minimize == False :
            indiv.fitG = np.mean(fitN) + np.std(fitN)

        return copy.deepcopy(indiv)

    def selectionFit(self, minimize):
        if minimize == True:
            if self.oldindiv.fit > self.indiv.fit:
                self.oldindiv = copy.deepcopy(self.indiv)
                return True
            else:
                return False
        else: # for maximization
            if self.oldindiv.fit < self.indiv.fit:
                self.oldindiv = copy.deepcopy(self.indiv)
                return True
            else:
                return False

    def selectionFitNeigh(self, minimize):
        if minimize == True :
            if self.oldindiv.fitG > self.indiv.fitG:
                self.oldindiv = copy.deepcopy(self.indiv)
                return True
            else:
                return False
        else: # maximization
            if self.oldindiv.fitG < self.indiv.fitG:
                self.oldindiv = copy.deepcopy(self.indiv)
                return True
            else:
                return False

    """ 
    for Walsh Local Search
    """
    def transWal(self):
        """
        translate bitstring represented Walsh terms into arrays of bits that they touches
        """
#        self.printW()
        self.WA = np.tile(Struct(arr = [], w = 0), len(self.model.w.keys())) # array representing Walsh terms
        c = 0
        for k in self.model.w.keys(): 
            a = []
            for i in range(self.dim):
                if k[i] == '1':
                    a.append(i)
            self.WA[c] = Struct(arr = a, w = self.model.w[k]) 
            c = c + 1
#            self.printWA()

    def binCount(self, arr, bit):
        """
        count the one bit of union self.WA[i].arr and bit
        """
        s = 0
        for i in arr:
            if bit[i] == '1':
                s = s + 1
        return s

    def initFitWal(self):
        """ 
        compute the sum array for the first time, according to the initial solution
        """
        self.sumArr = np.zeros(self.dim)
        self.WAS = np.tile(Struct(arr = [], w = 0), len(self.model.w.keys()))# Walsh coefficients with sign, represented in Array
        for i in range(len(self.WA)):
            W = int(math.pow(-1, self.binCount(self.WA[i].arr, self.indiv.bit))) * self.WA[i].w
            self.WAS[i] = Struct(arr = self.WA[i].arr, w = W)
            for j in self.WA[i].arr:
                self.sumArr[j] = self.sumArr[j] + W


    def initWal(self):
        """ 
        1. 
        compute the sum array for the first time, according to the initial solution
        2. 
        initialize a full matrix, Coincidence, C_ij, which is a symmetric one, and diagonal is empty (as we don't want to recompute SumArr), empty whenever i >= j.
        Two options for implementing C matrix
            a*. precompute C matrix 
            b. construct it on the fly using dict()
        """
        self.sumArr = np.zeros(self.dim)
        self.WAS = np.tile(Struct(arr = [], w = 0), len(self.model.w.keys()))# Walsh coefficients with sign, represented in Array
        self.lookup = dict()
        self.infectBit = dict()
        self.C = np.zeros((self.dim,self.dim)) # coincidence matrix
        #self.C = dict() # coincidence matrix
#        self.W = dict() # Walsh coefficient where sign is included, self.W should be updated as well 
        #        print "bit str", self.indiv.bit
        #        print 'bit', self.indiv.bit
        for i in range(len(self.WA)):
            W = int(math.pow(-1, self.binCount(self.WA[i].arr, self.indiv.bit))) * self.WA[i].w
            self.WAS[i] = Struct(arr = self.WA[i].arr, w = W)
            comb = self.genComb(len(self.WA[i].arr))

            for j in self.WA[i].arr:
                self.sumArr[j] = self.sumArr[j] + W

            for j in comb: # for each list in comb
                j0 = self.WA[i].arr[int(j[0])]
                j1 = self.WA[i].arr[int(j[1])]
#                if (j0, j1) in self.C.keys():
#                    self.C[j0,j1] = self.C[j0,j1] + W
#                elif W != 0:
#                    self.C[j0,j1] = W
                self.C[j0,j1] = self.C[j0,j1] + W

            # add list of order >= 3 Walsh terms for the purpose of updating C matrix
            if len(self.WA[i].arr) >= 3:
                for j in self.WA[i].arr:
                    if j not in self.infectBit: 
                        self.infectBit[j] = [Struct(arr=self.WA[i].arr, WI=i)]
                    else :
                        self.infectBit[j].append(Struct(arr=self.WA[i].arr, WI=i))

#        print 'C', self.C

#        print 'sum array', self.sumArr
        
#        self.sumArr = np.zeros(self.dim)
#        for k in self.model.w.keys(): # N * 2^K
#            self.W[k] = self.model.w[k] * math.pow(-1,wal.bc(k,self.indiv.bit))
#            for i in range(self.dim):
#                if k[i] == '1':
#                    self.sumArr[i] = self.sumArr[i] + self.W[k]

#        for i in zip(self.model.w.keys(), self.model.w.values(), self.W.values()):
#            print i
#
#        print 'sum array', self.sumArr

    def compPSum(self,bitStr):
        """
        use Elementary Landscape Analysis to obtain the average of neighs of given
        individual
        """
        p = np.zeros(self.model.k+1)
        for k in self.model.w.keys():
            oneC = k.count('1')
            if  oneC !=0 :
                p[oneC-1] = p[oneC-1] + self.model.w[k] * math.pow(-1,wal.bc(k,bitStr))
    #            else :
    #                p[0] = p[0] + self.model.w[k]

        s = 0
        for i in range(self.model.k+1):
            s = s + (i+1)*p[i]
        return s

    def genComb(self,N):
        """ 
        Generate C_k^0 sequence, index are stored, because they are more general, Implemented in an *incremental* fashion.
        """
        if N in self.lookup.keys(): # the key exists before
            return self.lookup[N]
        else : 
            comb =  []
            c = 0
            for i in range(N):
                for j in [ k for k in range(N) if k > i]:
                   arr = np.zeros(2)
                   arr[0] = i
                   arr[1] = j
                   comb.append(arr)
                   c = c + 1    
            self.lookup[N] = comb
            return comb

    def compCsum(self,p):
        """
        \sigma_{i=1}^{N} C_{ip}: be careful with C_{ii}, i \in N
        """
        s = 0

#        for i in self.C.keys():
#            if p in i: 
#                s = s + self.C[i]

        for i in range(p):
            s = s + self.C[i,p]

        for i in range(p+1, self.dim):
            s = s + self.C[p,i]

        s = s + self.sumArr[p]

        return s

    def updateFit(self,p):
        """
        super fast and simple update
        """
        self.sumArr[p] = - self.sumArr[p]

    def update(self, p):
        """
        By keeping track of coincidence matrix, 
        Cij stands for S_i(y_j) = S_i(x) - C_ij
        partially update the Sum Array and self.WAS, given the bit which is changed
        """
        self.sumArr[p] = - self.sumArr[p]
        
        for i in range(p): # i < p
#            if (i,p) in self.C.keys() and self.C[i,p]!=0:
#                self.sumArr[i] = self.sumArr[i] - 2*self.C[i,p]
            self.sumArr[i] = self.sumArr[i] - 2*self.C[i,p]

        for i in range(p+1, self.dim): # i > p
#            if (p,i) in self.C.keys() and self.C[p,i]!=0:
#                self.sumArr[i] = self.sumArr[i] - 2*self.C[p,i]
            self.sumArr[i] = self.sumArr[i] - 2*self.C[p,i]

        for i in range(p): # update C matrix, i < p
#            if (i,p) in self.C.keys() and self.C[i,p]!=0:
#                self.C[i,p] = - self.C[i,p]
            self.C[i,p] = - self.C[i,p]

        for i in range(p+1, self.dim): # i >p
#            if (p,i) in self.C.keys() and self.C[p,i]!=0:
#                self.C[p,i] = - self.C[p,i]
            self.C[p,i] = - self.C[p,i]

        # update the rest of elements in C matrix
        # pdb.set_trace()
        if p in self.infectBit.keys():
            for i in self.infectBit[p]:
                arr = copy.deepcopy(i.arr)
                arr.remove(p)
                comb = self.genComb(len(arr))
                for k in range(len(comb)):
                    k0 = arr[int(comb[k][0])]
                    k1 = arr[int(comb[k][1])]
                    self.C[k0,k1] = self.C[k0,k1] - 2 * self.WAS[i.WI].w

        for i in range(len(self.WAS)): # update WAS
            if p in self.WAS[i].arr :
                self.WAS[i].w = - self.WAS[i].w


    def neighWal(self):
        """ 
        generate neighborhoods and compute their fitnesses (mean of neighs) by Walsh Analysis
        """
        neighs = self.neighbors()
        neighIndiv = np.tile(Struct(fit = 0, fitG=0, bit = '0' ), (self.dim))
        for i in range(self.dim):
            neighIndiv[i] = Struct(fit = 0, bit = neighs[i]) 
        return neighIndiv

    def neighFitWal(self):
        """ 
        generate neighborhoods and compute their fitnesses (real one) by Walsh Analysis
        """
        neighs = self.neighbors()
        neighIndiv = np.tile(Struct(fit = 0,  bit = '0' ), (self.dim))
        for i in range(self.dim):
            neighIndiv[i] = Struct(fit = 0, bit = neighs[i]) 
        return neighIndiv

    def printW(self):
        """
        print all the original walsh terms
        """
        for k in self.model.w.keys():
            print k, self.model.w[k]

    def printWA(self):
        """
        print all walsh terms with array representation
        """
        for i in range(len(self.WA)):
           print self.WA[i].arr, self.WA[i].w

    def printC(self):
        a = 0
        c = 0
        for i in range(self.dim):
            for j in range(i+1,self.dim):
                a = a + 1
                if self.C[i,j] != 0:
                    print '1',
                    c = c + 1
                else:
                    print '0',
            print 
        print c,a, c/float(a)
