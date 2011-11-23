""" Local Search module:
    This local search module gets the problem (in terms of objective function),
    and return a solution found by hill climber """

import numpy as np
import random
import math
import copy
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
        self.oldindiv = self.initIndiv(self.dim)
        self.fitEval = 0
        self.oldindiv = self.evalPop(self.oldindiv)
        self.bsf = copy.deepcopy(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)
        while self.fitEval < self.MaxFit:
            neighs = self.neighbors()
            improveN = False
            for i in neighs:
                self.indiv.bit = np.copy(i)
                self.indiv = self.evalPop(self.indiv)
                if  self.selectionFit(minimize) == True:
                    improveN = True
            if improveN == False:
                self.restart('fit', minimize)
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit}


    def runNeigh(self,fitName, minimize):
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
                self.restart(fitName, minimize)
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

