import numpy as np
import random
import math
import pdb
import copy

class Struct:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

class CHC:
    def __init__(self):
        pass

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

    def run(self, func, MaxFit, popSize, dim, D, DR, M, fitName, minimize):
        if fitName == 'fit':
            return self.runFit(func, MaxFit, popSize, dim, D, DR, M, minimize)
        else :
            return self.runNeigh( func, MaxFit, popSize, dim, D, DR, M, fitName, minimize)

    def runFit(self, func, MaxFit, popSize, dim, D, DR, M, minimize):
        """
            func:   objective function
            MaxFit: maximum number of fitness evaluation
            popSize:population size 
            dim:    length of chromosome
            D:      Difference threshold
            DR:     Divergence rate
            M:      number OF best individuals for restarting
            minimize: minimize the problem or not
        """
        self.func = func
        self.MaxFit = MaxFit
        self.popSize = popSize
        self.dim = dim
        self.D = D
        self.Dinit = D
        self.M = M
        self.DR = DR
        if (MaxFit % popSize) == 0:
            self.numOfGen = math.floor(MaxFit/popSize) - 1
        else:
            self.numOfGen = math.floor(MaxFit/popSize) 
        self.pop = self.popInit(self.popSize, self.dim)
        gen = 0
        self.fitEval = 0
        allFit = []
        self.pop = self.evalPop(self.pop)
#        print 'initial pop:', [ self.pop[i].bit for i in range(len(self.pop)) ]
#        print 'initial fit:', [ self.pop[i].fit for i in range(len(self.pop)) ]
        while self.fitEval < self.MaxFit:
            gen = gen + 1
            matePool = self.reproduce()
            if self.D < 0:

#                print 'before restart, pop:', [ self.pop[i].bit for i in range(len(self.pop)) ]
#                print 'before restart, fit:', [ self.pop[i].fit for i in range(len(self.pop)) ]
                     
                self.divergePop(neigh = False, minimize = minimize)
                self.pop = self.evalPop(self.pop)
                self.D = self.Dinit

#                print 'after restart, pop:', [ self.pop[i].bit for i in range(len(self.pop)) ]
#                print 'after restart, fit:', [ self.pop[i].fit for i in range(len(self.pop)) ]
                continue

            offspring = self.HUXcrossover(matePool, neigh=False)
#            print 'offspring bits\n', [offspring[i].bit for i in range(len(offspring))]
            offspring = self.evalPop(offspring)
#            print 'offspring fit\n', [offspring[i].fit for i in range(len(offspring))]
            self.selectionFit(offspring, neigh = False, minimize=minimize)
#            print 'pop size', len(self.pop)
#            print 'pop fit after selection\n', [self.pop[i].fit for i in range(self.popSize)]
#            print 'pop bit after selection\n', [self.pop[i].bit for i in range(self.popSize)]
            if minimize == True :
                bestFit = min( [ self.pop[i].fit for i in range(len(self.pop)) ] )
            else :
                bestFit = max( [ self.pop[i].fit for i in range(len(self.pop)) ] )
#            print 'bestFit', bestFit
            allFit.append( bestFit )

#        print allFit
#        print 
#        print 
        return {'nEvals': self.fitEval, 'sol': min(allFit)}

    def runNeigh(self, func, MaxFit, popSize, dim, D, DR, M, fitName, minimize):
        """
            func:   objective function
            MaxFit: maximum number of fitness evaluation
            popSize:population size 
            dim:    length of chromosome
            D:      Difference threshold
            DR:     Divergence rate
            M:      number OF best individuals for restarting
            minimize: minimize the problem or not
        """
        self.func = func
        self.MaxFit = MaxFit
        self.popSize = popSize
        self.dim = dim
        self.D = D
        self.Dinit = D
        self.M = M
        self.DR = DR
        if (MaxFit % popSize) == 0:
            self.numOfGen = math.floor(MaxFit/popSize) - 1
        else:
            self.numOfGen = math.floor(MaxFit/popSize) 
        self.pop = self.popInitNeigh(self.popSize, self.dim)
        gen = 0
        self.fitEval = 0
        allFit = []
        allG = []
        self.pop = self.evalPopNeigh(self.pop, fitName, minimize)
#        print 'initial pop:', [ self.pop[i].bit for i in range(len(self.pop)) ]
#        print 'initial fit:', [ self.pop[i].fit for i in range(len(self.pop)) ]
#        print 'initial fitG:', [ self.pop[i].fitG for i in range(len(self.pop)) ]
        while self.fitEval < self.MaxFit:
            gen = gen + 1
#            print '********************GEN*********************', gen
            matePool = self.reproduce()
            if self.D < 0:

#                print 'before restart, pop:', [ self.pop[i].bit for i in range(len(self.pop)) ]
#                print 'before restart, fit:', [ self.pop[i].fit for i in range(len(self.pop)) ]
#                print 'before restart, G:', [ self.pop[i].fitG for i in range(len(self.pop)) ]
                     
                self.divergePop(neigh = True, minimize=minimize)
                self.pop = self.evalPopNeigh(self.pop, fitName, minimize)
                self.D = self.Dinit

#                print 'after restart, pop:', [ self.pop[i].bit for i in range(len(self.pop)) ]
#                print 'after restart, fit:', [ self.pop[i].fit for i in range(len(self.pop)) ]
#                print 'after restart, G:', [ self.pop[i].fitG for i in range(len(self.pop)) ]
                continue

            offspring = self.HUXcrossover(matePool, neigh=True)
#            print 'offspring bits\n', [offspring[i].bit for i in range(len(offspring))]
            offspring = self.evalPopNeigh(offspring, fitName, minimize)
#            print 'offspring fit\n', [offspring[i].fit for i in range(len(offspring))]
#            print 'offspring fitG\n', [offspring[i].fitG for i in range(len(offspring))]
            self.selectionFit(offspring, neigh=True, minimize=minimize)
#            print 'pop bit after selection\n', [self.pop[i].bit for i in range(self.popSize)]
#            print 'pop fit after selection\n', [self.pop[i].fit for i in range(self.popSize)]
#            print 'pop G after selection\n', [self.pop[i].fitG for i in range(self.popSize)]
            if minimize == True:
                bestFit = min( [ self.pop[i].fit for i in range(len(self.pop)) ] )
                bestG = min( [ self.pop[i].fitG for i in range(len(self.pop)) ] )
            else:
                bestFit = max( [ self.pop[i].fit for i in range(len(self.pop)) ] )
                bestG = max( [ self.pop[i].fitG for i in range(len(self.pop)) ] )
#            print 'bestFit', bestFit
#            print 'bestG', bestG
            allFit.append( bestFit )
            allG.append( bestG ) 

#        print allFit
#        print 
#        print 
        return {'nEvals': self.fitEval, 'sol': min(allFit), 'fitG': min(allG)}

    def reproduce(self):
        """ 
            only parents with hamming distance more than D are allowed for matting
        """
        random.shuffle(self.pop)
#        print 'pop after shuffle', [ self.pop[i].bit for i in range(len(self.pop)) ]
        matePool = []
        for i in range(self.popSize/2):
            pop0 = np.copy(self.pop[2*i].bit)
            pop1 = np.copy(self.pop[2*i+1].bit)
            diffbit = self.hammingDiff(pop0, pop1)
#            print 'pop0', pop0
#            print 'pop1', pop1
#            print 'diff bit', diffbit
            if len(diffbit)/2.0 > self.D:
                matePool.append(Struct( parents =  [pop0, pop1], bits = diffbit ))
        if not matePool:
            self.D = self.D - 1 
        
#        print 'matepool parents', [matePool[i].parents for i in range(len(matePool))]
#        print 'matepool diff bits', [matePool[i].bits for i in range(len(matePool))]
        return matePool

    def hammingDiff(self, str1, str2):
        diffbit = []
        for i in range(self.dim):
            if np.logical_xor(int(str1[i]), int(str2[i])) == True:
                diffbit.append(i)
        return diffbit

    def HUXcrossover(self, matePool, neigh):
        for i in range(len(matePool)):
            exchangebits = random.sample( matePool[i].bits, int(math.ceil(len(matePool[i].bits)/2.0)) )
            for j in exchangebits:
                if matePool[i].parents[0][j] == '0':
                    matePool[i].parents[0][j] = '1'
                    matePool[i].parents[1][j] = '0'
                else:
                    matePool[i].parents[0][j] = '0'
                    matePool[i].parents[1][j] = '1'
        # compose population from matePool
        if neigh == False:
            offspring = np.tile(Struct(fit = 0, bit = '0'), 2*len(matePool))
            for i in range(len(matePool)):
                offspring[2*i] = Struct( fit=0, bit = matePool[i].parents[0])
                offspring[2*i+1] = Struct( fit=0, bit = matePool[i].parents[1])
        else:
            offspring = np.tile(Struct(fit = 0, fitG = 0, bit = '0'), 2*len(matePool))
            for i in range(len(matePool)):
                offspring[2*i] = Struct( fit=0, fitG=0, bit = matePool[i].parents[0])
                offspring[2*i+1] = Struct( fit=0, fitG=0, bit = matePool[i].parents[1])
        return offspring
   
    def evalPop(self, pop):
        """ evaluate the population """
        for i in range(len(pop)):
            pop[i].fit = self.func(pop[i].bit)
        self.fitEval = self.fitEval + len(pop)
        return pop
    
    def evalPopNeigh(self, pop, fitName, minimize):
        """ evaluate the population with g(x) """
        # consider the real fitness
        for i in range(len(pop)):
            pop[i].fit = self.func(pop[i].bit)
        self.fitEval = self.fitEval + len(pop)
        # consider the fitness g(x)
        for i in range(len(pop)):
            fitN = np.zeros(self.dim)
            for j in range(self.dim):
                # flip the jth bit in bit-string
                neighStr = np.copy(pop[i].bit)
                if neighStr[j] == '1':
                    neighStr[j] = '0'
                else:
                    neighStr[j] = '1'
                fitN[j] = self.func(neighStr)
            if fitName == 'mean':
                self.pop[i].fitG = np.mean(fitN)
            elif minimize == True: 
                self.pop[i].fitG = np.mean(fitN) - np.std(fitN)
            elif minimize == False:
                self.pop[i].fitG = np.mean(fitN) + np.std(fitN)
        return pop

    def selectionFit(self, offspring, neigh, minimize):
        popAll = np.hstack((offspring, self.pop))
        if neigh == False :
            self.pop = np.copy (sorted(popAll,key=lambda p: p.fit, reverse = not minimize)[:self.popSize])
        else:
            self.pop = np.copy (sorted(popAll,key=lambda p: p.fitG, reverse = not minimize)[:self.popSize])
    
    def divergePop(self, neigh, minimize):
        # get the best individual
        if neigh == False:
            bestMpop = np.copy (sorted(self.pop, key=lambda p: p.fit, reverse = not minimize)[:self.M])
        else:
            bestMpop = np.copy (sorted(self.pop, key=lambda p: p.fitG, reverse = not minimize)[:self.M])
        self.pop[:self.M] = bestMpop
        bestIndiv = copy.deepcopy(bestMpop[0])
        for i in range(self.M, self.popSize):
            flipBits = random.sample(range(self.dim), int(self.DR*self.dim))
            self.pop[i] = copy.deepcopy(bestIndiv)
            for j in flipBits:
                if bestIndiv.bit[j] == '0':
                    self.pop[i].bit[j] = '1'
                else:
                    self.pop[i].bit[j] = '0'
