# Genetic algorithm
# g(x): the estimate best value of neighborhood, mean(N(x)) - std(N(x))
# crossover: one-point crossover with crossover rate
# mutation: one-bit flip mutation with mutation rate
# selection: deterministic tournament selection with size 2 / fitness-proportional selection
# Reference: Introduction to Evolutionary Computing - Genetic Algorithms by A.E. Eiben and J.E. Smith

import numpy as np
import random
import math
import pdb
import copy

class Struct:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

class GeneticAlgorithm:
    def __init__(self, func, MaxFit, popSize, dim):
        """ constructor """
        self.func = func
        self.MaxFit = MaxFit
        self.popSize = popSize
        self.dim = dim
        if (MaxFit % popSize) == 0:
            self.numOfGen = math.floor(MaxFit/popSize) - 1
        else:
            self.numOfGen = math.floor(MaxFit/popSize) 
    def popInit(self, popSize, dim):
        """ initial population with random bit string """
        pop = np.tile(Struct(fit = 0, bit = '0'), (self.popSize))
        for i in range(popSize):
            randBitStr = []
            for j in range(dim):
                if random.random()<0.5:
                    randBitStr.append('0')
                else:
                    randBitStr.append('1')
            pop[i] = Struct( fit = 0, bit = randBitStr)
        return pop
    def popInitNeigh(self, popSize, dim):
        """ initial population with random bit string """
        pop = np.tile(Struct(fit = 0, fitG = 0, bit = '0'), (self.popSize))
        for i in range(popSize):
            randBitStr = []
            for j in range(dim):
                if random.random()<0.5:
                    randBitStr.append('0')
                else:
                    randBitStr.append('1')
            pop[i] = Struct( fit = 0, fitG = 0, bit = randBitStr)
        return pop
    def run(self, crossoverR, mutationR,  fitName, tourSize = 2, minimize=True):
        if fitName == 'fit':
            return self.runFit(crossoverR, mutationR, tourSize, minimize)
        else :
            return self.runNeigh(crossoverR, mutationR, tourSize, fitName, minimize)

    def runFit(self, crossoverR, mutationR, tourSize, minimize):
        """ run the experiment with respect to real fitness """
        self.pop = self.popInit(self.popSize, self.dim)
        self.crossoverR = crossoverR
        self.mutationR = mutationR
        self.tourSize = tourSize
        gen = 0
        self.fitEval = 0
        allFit = np.zeros(self.numOfGen)
        self.evalPop()
        while self.fitEval < self.MaxFit:
            self.oldpop = np.copy(self.pop)
            gen = gen + 1
            self.crossover()
            self.mutation()
            self.evalPop()
            self.selectionFit(minimize)
            if minimize == True:
                bestFit = min( [ self.pop[i].fit for i in range(len(self.pop)) ] )
            else:
                bestFit = max( [ self.pop[i].fit for i in range(len(self.pop)) ] )
            allFit[gen-1] = bestFit
        return {'nEvals': self.fitEval, 'sol': min(allFit)}

    def runNeigh(self, crossoverR, mutationR, tourSize, fitName, minimize):
        """ run the experiment with g(x) """
        self.pop = self.popInit(self.popSize, self.dim)
        self.crossoverR = crossoverR
        self.mutationR = mutationR
        self.tourSize = tourSize
        gen = 0
        self.fitEval = 0
        allFit = np.zeros(self.numOfGen)
        allFitG = np.zeros(self.numOfGen)
        self.fitG = np.zeros(self.popSize)
        self.evalPopNeigh(fitName, minimize)
        while self.fitEval < self.MaxFit:
            self.oldpop = np.copy(self.pop)
            gen = gen + 1
            self.crossover()
            self.mutation()
            self.evalPopNeigh(fitName, minimize)
            self.selectionNeigh(minimize)
            if minimize == True:
                bestFit = min( [ self.pop[i].fit for i in range(len(self.pop)) ] )
                bestFitG = min( [ self.pop[i].fitG for i in range(len(self.pop)) ] )
            else:
                bestFit = max( [ self.pop[i].fit for i in range(len(self.pop)) ] )
                bestFitG = max( [ self.pop[i].fitG for i in range(len(self.pop)) ] )
            allFit[gen-1] = bestFit
            allFitG[gen-1] = bestFitG
        return {'nEvals': self.fitEval, 'sol': min(allFit),  'fitG': min(allFitG)}
    def evalPop(self):
        """ evaluate the population """
        for i in range(self.popSize):
            self.pop[i].fit = self.func(self.pop[i].bit)
        self.fitEval = self.fitEval + self.popSize 
    def evalPopNeigh(self, fitName, minimize):
        """ evaluate the population with g(x) """
        # consider the real fitness
        for i in range(self.popSize):
            self.pop[i].fit = self.func(self.pop[i].bit)
        self.fitEval = self.fitEval + self.popSize 
        # consider the fitness g(x)
        for i in range(self.popSize):
            fitN = np.zeros(self.dim)
            for j in range(self.dim):
                # flip the jth bit in bit-string
                neighStr = np.copy(self.pop[i].bit)
                if neighStr[j] == '1':
                    neighStr[j] = '0'
                else:
                    neighStr[j] = '1'
                fitN[j] = self.func(neighStr)
            if fitName == 'mean':
                self.pop[i].fitG = np.mean(fitN)
            elif minimize == True : 
                self.pop[i].fitG = np.mean(fitN) - np.std(fitN)
            elif minimize == False :
                self.pop[i].fitG = np.mean(fitN) + np.std(fitN)

    def mutation(self):
        """ one-bit flip mutation """
        for i in range(self.popSize):
            if random.random() < self.mutationR:
                flipBit = random.sample(range(self.dim),1)[0]
                if self.pop[i].bit[flipBit] == '1':
                    self.pop[i].bit[flipBit] = '0'
                else:
                    self.pop[i].bit[flipBit] = '1'
    def crossover(self):
        """ one point cross """
        # shuffle the population
        random.shuffle(self.pop)
        for i in range(self.popSize/2):
            if random.random() < self.crossoverR:
                pop0 = np.copy(self.pop[2*i].bit)
                pop1 = np.copy(self.pop[2*i+1].bit)
                crossPos = int(random.sample(range(self.dim - 1), 1)[0])
                for j in range(self.dim - crossPos - 1):
                    pop0[j+crossPos+1] = np.copy( self.pop[2*i+1].bit[j+crossPos+1] )
                    pop1[j+crossPos+1] = np.copy( self.pop[2*i].bit[j+crossPos+1] )
                self.pop[2*i].bit = np.copy(pop0)
                self.pop[2*i+1].bit = np.copy(pop1)

    def selectionFit(self, inminimize):
        self.tournamentSelection(neigh= False, minimize=inminimize)
    def selectionNeigh(self, inminimize):
        self.tournamentSelection(neigh= True, minimize=inminimize)

    def tournamentSelection(self, neigh, minimize):
        popAll = np.hstack((self.oldpop, self.pop))
        if neigh == False : 
            popSelect = np.tile(Struct(fit = 0, bit = '0'), (self.popSize))
            fit = [ popAll[i].fit for i in range(len(popAll)) ]
        else:
            popSelect = np.tile(Struct(fit = 0, fitG = 0, bit = '0'), (self.popSize))
            fit = [ popAll[i].fitG for i in range(len(popAll)) ]
        for i in range(self.popSize):
            tourI = random.sample(range(2*self.popSize), self.tourSize)
            tourPop = popAll[tourI]
            tourFit = [ fit[j] for j in tourI ]
            if minimize == True :
                bestI = tourFit.index(min(tourFit))
            else : 
                bestI = tourFit.index(max(tourFit))
            popSelect[i] = copy.deepcopy( popAll[tourI[bestI]] )
#            pdb.set_trace()

        self.pop = popSelect
            
    #############################################################################
    def proportionalSelectionFit():
        popAll = np.hstack((self.oldpop, self.pop))
        length = len(popAll)
        fitAll = [popAll[i].fit for i in range(length)]
        fitSum = sum(fitAll)
        
        prob = [(fitSum - fitAll[i])/((length-1)*fitSum) for i in range(length)]

        accumProb = np.zeros(length)
        for i in range(length):
            if i == 0:
                accumProb[i] = prob[i]
            else:
                accumProb[i] = accumProb[i-1] + prob[i]

        for i in range(length/2):
            selectedI = self.rouletteWheel(accumProb, length)
            self.pop[i] = popAll[selectedI]
    def proportionalSelectionNeigh():
        popAll = np.hstack((self.oldpop, self.pop))
        length = len(popAll)
        fitAll = [popAll[i].fitG for i in range(length)]
        fitSum = sum(fitAll)
        
        prob = [(fitSum - fitAll[i])/((length-1)*fitSum) for i in range(length)]

        accumProb = np.zeros(length)
        for i in range(length):
            if i == 0:
                accumProb[i] = prob[i]
            else:
                accumProb[i] = accumProb[i-1] + prob[i]

        for i in range(length/2):
            selectedI = self.rouletteWheel(accumProb, length)
            self.pop[i] = popAll[selectedI]
    def rouletteWheel(self, p, length):
        rand = random.random()
        if rand < p[0]:
            return 0
        for i in range(1, length):
            if rand < p[i]:
                return i
