# Genetic algorithm 
# g(x): the estimate best value of neighborhood, mean(N(x)) - std(N(x))
# crossover: one-point crossover with crossover rate
# mutation: one-bit flip mutation with mutation rate
# selection: truncation

import numpy as np
import random
import math
import pdb
from operator import itemgetter

class GeneticAlgorithm:
    def __init__(self, func, MaxFit, popSize, dim):
        """ constructor """
        self.func = func
        self.MaxFit = MaxFit
        self.popSize = popSize
        self.dim = dim
        self.fit = np.zeros(popSize)
        if (MaxFit % popSize) == 0:
            self.numOfGen = math.floor(MaxFit/popSize) - 1
        else:
            self.numOfGen = math.floor(MaxFit/popSize) 
        self.crossoverR = 0.5
        self.mutationR = 0.5
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
        """ run the experiment with respect to real fitness """
        self.pop = self.popInit(self.popSize, self.dim)
        gen = 0
        self.fitEval = 0
        allFit = np.zeros(self.numOfGen)
        self.evalPop()
        while self.fitEval < self.MaxFit:
            #pdb.set_trace()
            self.oldpop = np.copy(self.pop)
            self.oldfit = np.copy(self.fit)
            gen = gen + 1
            self.crossover()
            if random.random()<self.mutationR:
                self.mutation()
            self.evalPop()
            self.selectionFit()
            bestFit = min(self.fit)
            allFit[gen-1] = bestFit
        print 'allFit\n', allFit
        return min(allFit)
    def runNeigh(self):
        """ run the experiment with g(x) """
        self.pop = self.popInit(self.popSize, self.dim)
        gen = 0
        self.fitEval = 0
        allFit = np.zeros(self.numOfGen)
        allFitG = np.zeros(self.numOfGen)
        self.fitG = np.zeros(self.popSize)
        self.evalPopNeigh()
        while self.fitEval < self.MaxFit:
            #pdb.set_trace()
            self.oldpop = np.copy(self.pop)
            self.oldfit = np.copy(self.fit)
            self.oldfitG = np.copy(self.fitG)
            gen = gen + 1
            self.crossover()
            if random.random()<self.mutationR:
                self.mutation()
            self.evalPopNeigh()
            self.selectionNeigh()
            bestFit = min(self.fit)
            bestFitG = min(self.fitG)
            allFit[gen-1] = bestFit
            allFitG[gen-1] = bestFitG
        print 'all fit\n', allFit
        print 'all fit G\n', allFitG
        return min(allFit)
    def evalPop(self):
        """ evaluate the population """
        for i in range(self.popSize):
            self.fit[i] = self.func(self.pop[i])
        self.fitEval = self.fitEval + self.popSize 
    def evalPopNeigh(self):
        """ evaluate the population with g(x) """
        # consider the real fitness
        for i in range(self.popSize):
            self.fit[i] = self.func(self.pop[i])
        self.fitEval = self.fitEval + self.popSize 
        # consider the fitness g(x)
        for i in range(self.popSize):
            fitN = np.zeros(self.dim)
#            print 'fit of indiv', self.fit[i]
            for j in range(self.dim):
                # flip the jth bit in bit-string
                neighStr = np.copy(self.pop[i])
#                print 'original indiv ',i 
#                print neighStr
                if neighStr[j] == '1':
                    neighStr[j] = '0'
                else:
                    neighStr[j] = '1'
#                print 'neigh', neighStr
                fitN[j] = self.func(neighStr)
#            print 'fitness of neigh\n', fitN
           #self.fitG[i] = np.mean(fitN) + np.std(fitN)
            self.fitG[i] = np.mean(fitN) - np.std(fitN)
#            print 'g(x)', self.fitG[i]
#            print 
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
    def crossover(self):
        """ one point cross """
        for i in range( int(self.popSize * self.crossoverR /2.0) ):
            selectedIndex = random.sample(range(self.popSize), 2)
#            print 'selectedIndex', selectedIndex
            pop0 = np.copy(self.pop[selectedIndex[0]])
            pop1 = np.copy(self.pop[selectedIndex[1]])
#            print 'pop0', pop0
#            print 'pop1', pop1
            crossPos = int(random.sample(range(self.dim - 1), 1)[0])
#            print 'crossPos', crossPos
            for j in range(self.dim - crossPos - 1):
                pop0[j+crossPos+1] = np.copy( self.pop[selectedIndex[1]][j+crossPos+1] )
                pop1[j+crossPos+1] = np.copy( self.pop[selectedIndex[0]][j+crossPos+1] )
            self.pop[selectedIndex[0]] = np.copy(pop0)
            self.pop[selectedIndex[1]] = np.copy(pop1)
#            print 'self.pop0', pop0
#            print 'self.pop1', pop1
    def selectionFit(self):
        """ truncation selection based on fitness """
        # concatenate genotype of population
        popAll = np.concatenate( (self.oldpop, self.pop), axis = 0 )

        # concatenate phenotype of population
        oldfitIndex = np.vstack( (self.oldfit, range(self.popSize))  )
        fitIndex = np.vstack( (self.fit, range(self.popSize, 2*self.popSize)) )
        fitAll = np.hstack( (oldfitIndex, fitIndex) ) 

        # reshape array
        fitAll = [ [fitAll[0][i], fitAll[1][i]] for i in range(2*self.popSize) ]

        # sort on all fitnesses
        fitAll.sort()

        self.pop = [ list(popAll[fitAll[i][1]]) for i in range(self.popSize) ]
        self.fit = [ fitAll[i][0] for i in range(self.popSize) ]
    def selectionNeigh(self):
        """ selection base on the statistics of neighborhoods """
        # concatenate genotype of population
        popAll = np.concatenate( (self.oldpop, self.pop), axis = 0 )

        # concatenate fitness of population 
        fitnessAll = np.hstack( (self.oldfit, self.fit) ) 
#        print 'fitness All', fitnessAll

        # concatenate g(x) of population
        oldfitIndex = np.vstack( (self.oldfitG, range(self.popSize))  )
        fitIndex = np.vstack( (self.fitG, range(self.popSize, 2*self.popSize)) )
        fitAll = np.hstack( (oldfitIndex, fitIndex) ) 

        # reshape array
        fitAll = [ [fitAll[0][i], fitAll[1][i]] for i in range(2*self.popSize) ]

        # sort on all fitnesses
        fitAll.sort()
#        print 'fitAll index', fitAll
#        print 'sorted index', [  fitAll[i][1]  for i in range(self.popSize) ]

        self.pop = [ list( popAll[fitAll[i][1]] ) for i in range(self.popSize) ]
        for i in range(self.popSize):
            self.fit[i] = fitnessAll[fitAll[i][1]]
        self.fitG = [ fitAll[i][0] for i in range(self.popSize) ]
#        print 'old pop\n', self.oldpop
#        print 'pop\n', self.pop
#        print
#        print 'old fit\n', self.oldfit
#        print 'fit\n', self.fit
#        print
#        print 'old fitG\n', self.oldfitG
#        print 'fitG\n', self.fitG
