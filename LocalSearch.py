""" Local Search module:
    This local search module gets the problem (in terms of objective function),
    and return a solution found by hill climber """

import numpy as np
import WalshAnalysis as wal
import random
import math
import copy
import sys
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
        if fitName =='fit': 
            return self.runFit(minimize,restart)
        else :
            if compM != 'wal':
                return self.runNeigh(fitName, minimize,restart)
            else :
                return self.runWal(fitName, minimize,restart)

    def runWal(self,fitName,  minimize, restart):
        """ 
        steepest descent local search with respect to mean of neighs by Walsh Analysis
        """
        self.oldindiv = self.initIndivNeigh(self.dim)
        self.fitEval = 0
        self.oldindiv = self.evalPop(self.oldindiv)
        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * self.compPSum(self.oldindiv.bit)
#        print 'initial', self.oldindiv.bit, 'fitG', self.oldindiv.fitG, 'fit', self.oldindiv.fit
        self.bsf = copy.deepcopy(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)

        self.compSumArr()
#        self.trace = [Struct(fitEval= self.fitEval,fit = self.oldindiv.fit, fitG = self.oldindiv.fitG)]
        while self.fitEval < self.MaxFit:
            # generate neighborhood and compute their fitness
            neighPop = self.neighWal()
            # compute the fitG (mean) of each neighborhood individuals
            improveN = False
            nCount = 0
#            print 'current', self.oldindiv.bit, 'fitG', self.oldindiv.fitG, 'fit', self.oldindiv.fit
            for n in neighPop:
                n = self.evalPop(n)
                n.fitG = n.fit - 2/float(self.dim) * self.compPSum(n.bit) 
                self.indiv = copy.deepcopy(n)
                if self.selectionFitNeigh(minimize) == True:
                    improveN = True
                    changeBit = nCount

#            print 'improveN', improveN
#            self.trace.append(Struct(fitEval= self.fitEval,fit = self.oldindiv.fit, fitG = self.oldindiv.fitG))
            if improveN == False:
                if restart == True:
                    self.restart(fitName, minimize)
                else:
#                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit,'trace':self.trace}
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                self.updateSumArr(changeBit)

            nCount = nCount + 1

#        return { 'nEvals': self.fitEval, 'sol': self.bsf.fit, 'fitG': self.bsf.fitG, 'bit':self.bsf.bit,'trace':self.trace}
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
            for i in neighs:
                self.indiv.bit = np.copy(i)
                self.indiv = self.evalPopNeigh(self.indiv, fitName, minimize)
                if self.selectionFitNeigh(minimize) == True:
                    improveN = True
            if improveN == False:
                if restart == True:
                    self.restart(fitName, minimize)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
        return { 'nEvals': self.fitEval, 'sol': self.bsf.fit, 'fitG': self.bsf.fitG, 'bit':self.bsf.bit}

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

        if fitName == 'fit':
            self.oldindiv = self.initIndiv(self.dim)
            self.oldindiv = self.evalPop(self.oldindiv)
        else :
            self.oldindiv = self.initIndivNeigh(self.dim)
            self.oldindiv = self.evalPopNeigh(self.oldindiv, fitName, minimize)

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

    def updateSumArr(self, changeBit):
        """
        partially update the Sum Array and self.W, given the bit which is changed
        """
        for k in self.model.w.keys():
            if k[changeBit] == '1':
                self.sumArr[changeBit] = self.sumArr[changeBit] - 2*(self.W[k])
                for i in [i for i in range(self.dim) if i != changeBit]:
                    if k[i] == '1':
                        self.sumArr[i] = self.sumArr[i] - 2*(self.W[k])

        for k in self.W.keys():
            if k[changeBit] == '1':
                self.W[k] = - self.W[k]

    def compSumArr(self):
        """ 
        compute the sum array for the first time, according to the initial 
        solution
        """
        self.sumArr = np.zeros(self.dim)

        self.W = dict() # Walsh coefficient where sign is included, self.W should be updated as well 
#        print "bit str", self.indiv.bit
        for k in self.model.w.keys():
            self.W[k] = self.model.w[k] * math.pow(-1,wal.bc(k,self.indiv.bit))
            for i in range(self.dim):
                if k[i] == '1':
                    self.sumArr[i] = self.sumArr[i] + self.W[k]

#        for i in zip(self.model.w.keys(), self.model.w.values(), self.W.values()):
#            print i

#        print 'sum array', self.sumArr

    def neighWal(self):
        """ 
        generate neighborhoods and compute their fitnesses by Walsh Analysis
        """
        neighs = self.neighbors()
        neighIndiv = np.tile(Struct(fit = 0, fitG=0, bit = '0' ), (self.dim))
        for i in range(self.dim):
            neighIndiv[i] = Struct(fit = self.indiv.fit-2*self.sumArr[i], bit = neighs[i]) 
        return neighIndiv
