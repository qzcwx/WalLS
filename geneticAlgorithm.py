# Genetic algorithm with:

import numpy as np
import random
import pdb
from operator import itemgetter

class GeneticAlgorithm:
    def __init__(self, func, MaxFit, popSize, dim):
        """ constructor """
        self.func = func
        self.MaxFit = MaxFit
        self.popSize = popSize
        self.dim = dim
        self.pop = self.popInit(popSize, dim)
        self.fit = np.zeros(popSize)
    def popInit(self, popSize, dim):
        """ initial population with random bit string """
        pop = []
        for i in range(popSize):
            randBitStr = []
            for j in range(dim):
                if random.random()<0.5:
                    randBitStr.append('0')
                else:
                    randBitStr.append('1')
            pop.append(randBitStr)
        return np.array(pop)
    def run(self):
        """ run the experiment """
        gen = 0
        self.fitEval = 0
        self.evalPop()
        print 'initial', self.fit
        while self.fitEval < self.MaxFit:
            #pdb.set_trace()
            self.oldpop = np.copy(self.pop)
            self.oldfit = np.copy(self.fit)
            gen = gen + 1
            self.mutation()
            self.evalPop()
            self.selectionFit()
            print 'bestVal', min(self.fit)
            print 
    def evalPop(self):
        """ evaluate the population """
        for i in range(self.popSize):
            self.fit[i] = self.func(self.pop[i])
        self.fitEval = self.fitEval + self.popSize 
    def mutation(self):
        """ one-bit flip mutation """
        #print self.pop
        for i in range(self.popSize):
            flipBit = random.sample(range(self.dim),1)
            flipBit = flipBit[0]
            #print 'flip bit', flipBit
            if self.pop[i][flipBit] == '1':
                self.pop[i][flipBit] = '0'
            else:
                self.pop[i][flipBit] = '1'
    def selectionFit(self):
        """ truncation selection based on fitness """
#        print 'fit', self.fit
        # concatenate genotype of population
#        print 'old\n', self.oldpop
#        print 'pop\n', self.pop
        popAll = np.concatenate( (self.oldpop, self.pop), axis = 0 )
        #print 'all\n', popAll

        # concatenate phenotype of population
        oldfitIndex = np.vstack( (self.oldfit, range(self.popSize))  )
        #print 'oldfitIndex', oldfitIndex
        fitIndex = np.vstack( (self.fit, range(self.popSize, 2*self.popSize)) )
        #print 'fitIndex', fitIndex
        fitAll = np.hstack( (oldfitIndex, fitIndex) ) 
        #print 'fitAll', fitAll

        # reshape array
        fitAll = [ [fitAll[0][i], fitAll[1][i]] for i in range(2*self.popSize) ]

        # sort on 
#        tempFit = np.copy(fitAll)
#        tempFit.sort()
        fitAll.sort()
#        print 'after tempFit', tempFit
#        print 'after fitAll', fitAll

        self.pop = [ list(popAll[fitAll[i][1]]) for i in range(self.popSize) ]
        self.fit = [ fitAll[i][0] for i in range(self.popSize) ]

#        print 'self.pop', self.pop
#        print 'self.oldfit', self.oldfit
#        print 'self.fit', self.fit
