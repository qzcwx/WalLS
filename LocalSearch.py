""" Local Search module:
    This local search module gets the problem (in terms of objective function),
    and return a solution found by hill climber """

import numpy as np
import random
import math
import pdb

class Struct:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

class LocalSearch:
    def __init__(self, func, MaxFit, dim):
        self.func = func
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

    def run(self, fitName, minimize=True):
        if fitName =='fit': 
            return self.runFit(minimize)
        else :
            return self.runNeigh(fitName, minimize)

    def runFit(self, minimize):
        self.indiv = self.initIndiv(self.dim)
        self.fitEval = 0
        self.evalPop()
        self.oldindiv = Struct( fit = self.indiv.fit, bit = self.indiv.bit)
        while self.fitEval < self.MaxFit:
            neighs = self.neighbors()
            improveN = False
            for i in neighs:
                self.indiv.bit = np.copy(i)
                self.evalPop()
                if  self.selectionFit(minimize) == True:
                    improveN = True
            if improveN == False:
                return {'nEvals': self.fitEval, 'sol': self.oldindiv.fit}
        return {'nEvals': self.fitEval, 'sol': self.oldindiv.fit}

    def runNeigh(self,fitName, minimize):
        self.indiv = self.initIndivNeigh(self.dim)
        self.fitEval = 0
        self.evalPopNeigh(fitName, minimize)
        self.oldindiv = Struct( fit = self.indiv.fit, fitG = self.indiv.fitG, bit = self.indiv.bit )
        while self.fitEval < self.MaxFit:
            neighs = self.neighbors()
            improveN = False
            for i in neighs:
                self.indiv.bit = np.copy(i)
                self.evalPopNeigh(fitName, minimize)
                if self.selectionFitNeigh(minimize) == True:
                    improveN = True
            if improveN == False:
                return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG }
        return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG }

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

    def evalPop(self):
        self.indiv.fit = self.func(self.indiv.bit)
        self.fitEval = self.fitEval + 1

    def evalPopNeigh(self, fitName, minimize):
        """ evaluate the individual itself """
        self.indiv.fit = self.func(self.indiv.bit)
        self.fitEval = self.fitEval + 1
        """ evaluate all neighborhood """
        fitN = np.zeros(self.dim)
        for j in range(self.dim):
            # flip the jth bit in bit-string
            neighStr = np.copy(self.indiv.bit)
            if neighStr[j] == '1':
                neighStr[j] = '0'
            else:
                neighStr[j] = '1'
            fitN[j] = self.func(neighStr)
        if fitName == 'mean':
            self.indiv.fitG = np.mean(fitN)
        elif minimize == True :
            self.indiv.fitG = np.mean(fitN) - np.std(fitN)
        elif minimize == False :
            self.indiv.fitG = np.mean(fitN) + np.std(fitN)

    def selectionFit(self, minimize):
        if minimize == True:
            if self.oldindiv.fit > self.indiv.fit:
                self.oldindiv = Struct( fit = self.indiv.fit, bit = self.indiv.bit)
                return True
            else:
                return False
        else: # for maximization
            if self.oldindiv.fit < self.indiv.fit:
                self.oldindiv = Struct( fit = self.indiv.fit, bit = self.indiv.bit)
                return True
            else:
                return False

    def selectionFitNeigh(self, minimize):
        if minimize == True :
            if self.oldindiv.fitG > self.indiv.fitG:
                self.oldindiv = Struct( fit = self.indiv.fit, fitG = self.indiv.fitG, bit = self.indiv.bit )
                return True
            else:
                return False
        else: # maximization
            if self.oldindiv.fitG < self.indiv.fitG:
                self.oldindiv = Struct( fit = self.indiv.fit, fitG = self.indiv.fitG, bit = self.indiv.bit )
                return True
            else:
                return False
