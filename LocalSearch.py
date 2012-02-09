""" 
Local Search module:
    This local search module gets the problem (in terms of objective function),
    and return a solution found by hill climber
"""

import numpy as np
#import matplotlib.pyplot as plt
import WalshAnalysis as wal
import random
import math
import copy
import sys
import time
import pdb
from sets import Set

class Struct:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

class LocalSearch:
    def __init__(self, model, MaxFit, dim):
        self.func = model.compFit
        self.model = model
        self.MaxFit = MaxFit
        self.dim = dim
        self.threshold = 1e-15

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

    def run(self, fitName, minimize, restart, compM ):
        if compM == 'bf':
            if fitName =='fit': 
                return self.runFit(minimize,restart)
            else:
                return self.runNeigh(fitName, minimize,restart)
        elif compM == 'walWalk':
            if fitName == 'fit':
                return self.runFitS2walk(fitName, minimize, restart)
            elif fitName == 'mean':
                return self.runMeanSC2walk(fitName, minimize, restart)
                #return self.runMeanWal(fitName, minimize, restart)
        elif compM == 'walRest':
            if fitName == 'fit':
                #TODO compose the method for restart
                return self.runFitS2rest(fitName, minimize, restart)
            elif fitName == 'mean':
                return self.runMeanSCrest(fitName, minimize, restart)
        elif compM == 'supm':
            if fitName == 'fit':
                return self.runFitsm(fitName, minimize, restart)
            elif fitName == 'mean':
                return self.runMeansm(fitName, minimize, restart)
        elif compM == 'bitImp':
            return self.bitImpact()


    def runFitsm(self,fitName, minimize, restart):
        """ 
        steepest descent local search running on S, with supermove enable
        """
        self.fitEval = 0

        start = time.time()
        self.transWal()
        self.oldindiv = self.initIndiv(self.dim)
        self.oldindiv = self.evalPop(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)
        self.initWal()
        #self.printInter()
        
        self.bsf = copy.deepcopy(self.oldindiv)
        
        self.WA = []

        #self.printInter()
        
        init = False
        updateT = 0
        initT = time.time() - start
        start = time.time()
        while self.fitEval < self.MaxFit:
            if init == False:
                improveN, bestI = self.genFitBestsm(minimize)
                init = True
            else:
                improveN, bestI = self.updateFitBestsm(minimize)

#            print improveN, bestI
#            print self.sumArr
#            print self.Buffer
        
            if improveN == False:
                if restart == True:
                    updateT = updateT + time.time() - start
                    startR = time.time()
                    oldbit = self.oldindiv.bit
                    self.oldindiv = self.evalPop(self.oldindiv)
#                    print 'restart'

                    self.fitEval = self.fitEval - 1
                    self.restart(fitName, minimize, False)
                    init = False
                    diff = self.diffBits(oldbit, self.oldindiv.bit)
                    
                    for i in diff:
                        self.update(i)
                        self.updateWAS(i)
                    initT = initT + time.time() - startR
                    start = time.time()
                else:
                    self.oldindiv = self.evalPop(self.oldindiv)
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
#                print self.P
                self.fitEval = self.fitEval + len(self.P) * self.dim
                for i in self.P:
                    self.update(i)
                    self.updateWAS(i)
                    #takenBits[bestI] = True
    #                updateCT = updateCT + time.time() - start
    #                bestBitsCount[bestI] = bestBitsCount[bestI] + 1
                    if self.oldindiv.bit[i] == '1':
                        self.oldindiv.bit[i] = '0'
                    else:
                        self.oldindiv.bit[i] = '1'

        self.bsf = self.evalPop(self.bsf)
        updateT = updateT + time.time() - start
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit,'init':initT, 'update':updateT}

    def runFitSwalk(self,fitName, minimize, restart):
        """ 
        steepest descent local search running on S
        """
        self.fitEval = 0
        start = time.time()
        self.transWal()
        self.oldindiv = self.initIndiv(self.dim)
        self.oldindiv = self.evalPop(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)
        self.initWal()
        self.bsf = copy.deepcopy(self.oldindiv)
        self.WA = []

        walkLen = 10
        init = False
        updateT = 0
        initT = time.time() - start
        start = time.time()
        while self.fitEval < self.MaxFit:
            if init == False:
                improveN, bestI, evalCount = self.genFitBest(minimize)
                init = True
            else:
                improveN, bestI, evalCount = self.updateFitBest(bestI,minimize)

            self.fitEval = self.fitEval + evalCount
        
            if improveN == False:
                if restart == True:
                    updateT = updateT + time.time() - start
                    startR = time.time()
                    self.oldindiv = self.evalPop(self.oldindiv)

                    diff = self.walk(fitName, minimize,False, walkLen)
                    init = False

                    for i in diff:
                        self.update(i)
                        self.updateWAS(i)
                    initT = initT + time.time() - startR
                    start = time.time()
                else:
                    self.oldindiv = self.evalPop(self.oldindiv)
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                self.update(bestI)
                self.updateWAS(bestI)
                if self.oldindiv.bit[bestI] == '1':
                    self.oldindiv.bit[bestI] = '0'
                else:
                    self.oldindiv.bit[bestI] = '1'
        self.bsf = self.evalPop(self.bsf)
        updateT = updateT + time.time() - start
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit,'init':initT, 'update':updateT}

    def runFitS2walk(self,fitName, minimize, restart):
        """ 
        steepest descent local search running on S
        """
        self.fitEval = 0
        start = time.time()
        self.transWal()
        self.oldindiv = self.initIndiv(self.dim)
        self.oldindiv = self.evalPop(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)
        self.initWal()
        self.bsf = copy.deepcopy(self.oldindiv)
        self.WA = []

        walkLen = 10
        init = False
        updateT = 0

        initT = time.time() - start
        start = time.time()

        while self.fitEval < self.MaxFit:
            if init == False:
                improveN, bestI, evalCount = self.genFitBest2(minimize)
                init = True
            else:
                improveN, bestI, evalCount = self.updateFitBest2(bestI,minimize)
            self.fitEval = self.fitEval + evalCount

#            print 'oldindiv',self.oldindiv.bit
#            print 'improveA',self.improveA
#            print improveN, bestI, self.fitEval
#            print self.oldindiv.bit
#            print improveN, bestI, evalCount, self.fitEval, 
#            print 
#            pdb.set_trace()
        
            if improveN == False:
                if restart == True:
                    updateT = updateT + time.time() - start
                    startR = time.time()
                    self.oldindiv = self.evalPop(self.oldindiv)
                    diff = self.walk(fitName, minimize,False, walkLen)
                    init = False

#                    print self.bsf.fit
#                    pdb.set_trace()

                    for i in diff:
                        self.update(i)
                        self.updateWAS(i)
                    initT = initT + time.time() - startR
                    start = time.time()
                else:
                    self.oldindiv = self.evalPop(self.oldindiv)
                    return {'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.bsf.bit,'init':initT, 'update':updateT}
            else : # improveN is TRUE 
                for i in bestI:
                    self.update(i)
                    self.updateWAS(i)
                    if self.oldindiv.bit[i] == '1':
                        self.oldindiv.bit[i] = '0'
                    else:
                        self.oldindiv.bit[i] = '1'

        if fitName == 'fit' and minimize == True :
            if self.bsf.fit > self.oldindiv.fit:
                self.bsf = copy.deepcopy(self.oldindiv)
        elif fitName == 'fit' and minimize == False :
            if self.bsf.fit < self.oldindiv.fit:
                self.bsf = copy.deepcopy(self.oldindiv)
        updateT = updateT + time.time() - start
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit,'init':initT, 'update':updateT}

    def runFitWal(self,fitName, minimize, restart):
        """ 
        Walsh Analysis for speeding up steepest descent local search
        """
        self.oldindiv = self.initIndiv(self.dim)
        self.fitEval = 0
        
        self.transWal()
        self.indiv = copy.deepcopy(self.oldindiv)
        self.initWal()
#        self.initInter()

        self.oldindiv = self.evalPop(self.oldindiv)
        self.bsf = copy.deepcopy(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)

        self.WA = []
#        print 'C', self.C
#        self.trace = [Struct(fitEval= self.fitEval,fit = self.oldindiv.fit, fitG = self.oldindiv.fitG)]
        init = False
        while self.fitEval < self.MaxFit:
            # generate neighborhood and compute their fitness
            neighPop = self.neighWal()
            # compute the fitG (mean) of each neighborhood individuals
            improveN = False
            nCount = 0
            #print 
            #print 'current', self.oldindiv.bit, 'fit', self.oldindiv.fit
            oldFit = self.oldindiv.fit
#            if init == False: 
#                # for initialization, all neighs should be evaluated
#                self.fitArr = np.zeros(self.dim)
            for i in range(self.dim):
                self.indiv = copy.deepcopy(neighPop[i])
                self.indiv.fit = oldFit - 2*self.sumArr[nCount]
#                self.fitArrp[i] = self.indiv.fit
                self.fitEval = self.fitEval + 1
                #print 'neigh: ', self.indiv.bit, 'fit', self.indiv.fit, 'fitG', self.indiv.fitG
                if self.selectionFit(minimize) == True:
                    #print 'better neigh!'
                    improveN = True
                    changeBit = nCount

                nCount = nCount + 1
#            else:
#                for i in range(lastChangeBit): 

#            print 'improveN', improveN
            #print self.fitEval
            #pdb.set_trace()
#            self.trace.append(Struct(fitEval= self.fitEval,fit = self.oldindiv.fit, fitG = self.oldindiv.fitG))
            if improveN == False:
                if restart == True:
                    oldbit = self.oldindiv.bit
                    self.restart(fitName, minimize)
                    newbit = self.oldindiv.bit
                    #print oldbit, newbit
                    diff = self.diffBits(oldbit, newbit)
                    for i in diff:
                        self.update(i)
                else:
#                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit,'trace':self.trace}
                    #print 'compPSum', compPSumT
                    #print 'update', updateT
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                self.update(changeBit)
                lastChangeBit = changeBit

            if init == False:
                init = True

        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit}


    def runFitSrest(self,fitName, minimize, restart):
        """ 
        steepest descent local search running on S
        """
        self.fitEval = 0
        
        start = time.time()
        self.transWal()
        self.oldindiv = self.initIndiv(self.dim)
        self.oldindiv = self.evalPop(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)
        self.initWal()
        
        self.bsf = copy.deepcopy(self.oldindiv)
        #bestBitsCount = np.zeros(self.dim)
        
        self.WA = []
        #takenBits = self.dim*[False]
        
        init = False
        #lenImproveA = []
        updateT = 0
        initT = time.time() - start
#        genT = 0
#        updateT = 0
#        updateCT = 0
#        endT = 0
#        resT = 0
#        bufferSize = 5
#        withinLimit = 0
#        exceedLimit = 0
        start = time.time()
        while self.fitEval < self.MaxFit:
            if init == False:
                #start = time.time()
                improveN, bestI, evalCount = self.genFitBest(minimize)
                #genT = genT + time.time() - start
                init = True
            else:
#                start = time.time()
                improveN, bestI, evalCount = self.updateFitBest(bestI,minimize)
            self.fitEval = self.fitEval + evalCount
#                updateT = updateT + time.time() - start
            #lenImproveA.append( len(self.improveA) )
#            sBuffer = 0
#            sBufferNotChange = 0
#            for i in self.improveA:
#                if bestI != None :
#                    if (bestI in self.Inter) and (i not in self.Inter[bestI].arr):
#                        sBuffer = sBuffer + 1
#                        if takenBits[i] == False:
#                            sBufferNotChange = sBufferNotChange + 1
#                    elif bestI not in self.Inter:
#                        sBuffer = sBuffer + 1
#                        if takenBits[i] == False:
#                            sBufferNotChange = sBufferNotChange + 1
#            #print sBuffer, len(self.improveA)
#            if sBuffer <= bufferSize:
#               withinLimit = withinLimit + 1
#            else:
#               exceedLimit = exceedLimit + 1
        
            if improveN == False:
                if restart == True:
                    updateT = updateT + time.time() - start
                    startR = time.time()
                    oldbit = self.oldindiv.bit
                    self.oldindiv = self.evalPop(self.oldindiv)

#                    start = time.time()
                    self.fitEval = self.fitEval - 1
                    self.restart(fitName, minimize, False)
                    init = False
                    diff = self.diffBits(oldbit, self.oldindiv.bit)
#                    resT = resT + time.time() - start
                    
                    #takenBits = self.dim*[False]
                    #print 'restart'
                    for i in diff:
                        self.update(i)
                        self.updateWAS(i)
                    initT = initT + time.time() - startR
                    start = time.time()
                else:
                    self.oldindiv = self.evalPop(self.oldindiv)
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
#                start = time.time()
                self.update(bestI)
                self.updateWAS(bestI)
                #takenBits[bestI] = True
#                updateCT = updateCT + time.time() - start
#                bestBitsCount[bestI] = bestBitsCount[bestI] + 1
                if self.oldindiv.bit[bestI] == '1':
                    self.oldindiv.bit[bestI] = '0'
                else:
                    self.oldindiv.bit[bestI] = '1'

#        print bestBitsCount 
#        print self.InterCount
#        plt.plot(bestBitsCount,'o')
#        plt.plot(self.InterCount,'o')
#        plt.savefig('n'+str(self.model.n)+'k'+str(self.model.k)+'.png')
#        print np.corrcoef([bestBitsCount,self.InterCount])[1,0]
        
#        start = time.time()
        self.bsf = self.evalPop(self.bsf)
#        endT = endT + time.time() - start
        #print np.mean(lenImproveA), np.std(lenImproveA)
#        print 'genBest', genT
#        print 'updateBest', updateT
#        print 'restart', resT
#        print 'updateCT',updateCT
#        print 'endT', endT
        updateT = updateT + time.time() - start
        #print 'Percentage of be withing buffer size limit', bufferSize, ':', withinLimit / (withinLimit + exceedLimit + 0.0)
        #print  withinLimit / (withinLimit + exceedLimit + 0.0)
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit,'init':initT, 'update':updateT}

    def runMeanSCwalk(self,fitName, minimize, restart):
        """ 
        steepest descent local search with respect to mean of neighs by Walsh Analysis
        """
        self.fitEval = 0
        
        start = time.time()
        self.transWal()
        self.oldindiv = self.initIndivNeigh(self.dim)
        self.oldindiv = self.evalPop(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)
        self.initWal()
        
        self.initSC()
        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (sum(self.sumArr))
        self.indiv.fitG = self.oldindiv.fitG

        self.bsf = copy.deepcopy(self.oldindiv)
        self.WA = []

        init = False
        updateT = 0
        walkLen = 10
        initT = time.time() - start
        start = time.time()
        while self.fitEval < self.MaxFit:

            if init == False:
                improveN, bestI, evalCount = self.genMeanBest(minimize)
                init = True
            else :
                improveN, bestI, evalCount = self.updateMeanBest(bestI,minimize)

            self.fitEval = self.fitEval + evalCount

            if improveN == False:
                if restart == True:
                    updateT = updateT + time.time() - start
                    startR = time.time()
                    self.oldindiv = self.evalPop(self.oldindiv)
                    self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.sumArr))
                    diff = self.walk( fitName, minimize, False, walkLen )
                    init = False
                    
                    for i in diff:
                        self.update(i)
                        self.updateSC(i)
                        self.updateWAS(i)
                    initT = initT + time.time() - startR
                    start = time.time()
                else:
                    self.oldindiv = self.evalPop(self.oldindiv)
                    self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.sumArr))
                    self.fitEval = self.fitEval - 1
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                self.update(bestI)
                self.updateSC(bestI)
                self.updateWAS(bestI)
                if self.oldindiv.bit[bestI] == '1':
                    self.oldindiv.bit[bestI] = '0'
                else:
                    self.oldindiv.bit[bestI] = '1'

        self.bsf = self.evalPop(self.bsf)
        self.fitEval = self.fitEval - 1
        diff = self.diffBits(self.bsf.bit, self.oldindiv.bit)
        for i in diff:
            self.update(i)
        self.bsf.fitG = self.bsf.fit - 2/float(self.dim) * (np.sum(self.sumArr))
        updateT = updateT + time.time() - start
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'fitG': self.bsf.fitG, 'bit':self.bsf.bit,'init':initT, 'update':updateT}

    def runMeanSC2walk(self,fitName, minimize, restart):
        """ 
        steepest descent local search with respect to mean of neighs by Walsh Analysis
        """
        self.fitEval = 0
        
        start = time.time()
        self.transWal()
        self.oldindiv = self.initIndivNeigh(self.dim)
        self.oldindiv = self.evalPop(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)
        self.initWal()
        
        self.initSC()
        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (sum(self.sumArr))
        self.indiv.fitG = self.oldindiv.fitG

        self.bsf = copy.deepcopy(self.oldindiv)
        self.WA = []

        init = False
        updateT = 0
        walkLen = 10

        initT = time.time() - start
        start = time.time()

        while self.fitEval < self.MaxFit:

            if init == False:
                improveN, bestI, evalCount = self.genMeanBest2(minimize)
                init = True
            else :
                improveN, bestI, evalCount = self.updateMeanBest2(bestI,minimize)

            self.fitEval = self.fitEval + evalCount

            if improveN == False:
                if restart == True:
                    updateT = updateT + time.time() - start
                    startR = time.time()
                    self.oldindiv = self.evalPop(self.oldindiv)
                    self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.sumArr))
                    diff = self.walk( fitName, minimize, False, walkLen )
                    init = False
                    
                    for i in diff:
                        self.update(i)
                        self.updateSC(i)
                        self.updateWAS(i)
                    initT = initT + time.time() - startR
                    start = time.time()
                else:
                    self.oldindiv = self.evalPop(self.oldindiv)
                    self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.sumArr))
                    self.fitEval = self.fitEval - 1
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                for i in bestI:
                    self.update(i)
                    self.updateSC(i)
                    self.updateWAS(i)
                    if self.oldindiv.bit[i] == '1':
                        self.oldindiv.bit[i] = '0'
                    else:
                        self.oldindiv.bit[i] = '1'

        self.bsf = self.evalPop(self.bsf)
        self.fitEval = self.fitEval - 1
        diff = self.diffBits(self.bsf.bit, self.oldindiv.bit)
        for i in diff:
            self.update(i)
        self.bsf.fitG = self.bsf.fit - 2/float(self.dim) * (np.sum(self.sumArr))
        updateT = updateT + time.time() - start
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'fitG': self.bsf.fitG, 'bit':self.bsf.bit,'init':initT, 'update':updateT}


    def runMeanSCrest(self,fitName, minimize, restart):
        """ 
        steepest descent local search with respect to mean of neighs by Walsh Analysis
        """
#        startA = time.time()
        self.fitEval = 0
        
        start = time.time()
        self.transWal()
        self.oldindiv = self.initIndivNeigh(self.dim)
        self.oldindiv = self.evalPop(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)
        self.initWal()
        
#        start = time.time()
        self.initSC()
        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (sum(self.sumArr))
        self.indiv.fitG = self.oldindiv.fitG
#        print '+ initial', time.time()-start

        self.bsf = copy.deepcopy(self.oldindiv)
        takenBits = self.dim*[False]
        self.WA = []

        init = False
        updateT = 0
        initT = time.time() - start
        #lenImproveA = []
#        genT = 0
#        updateT = 0
#        updateCT = 0
#        fitGT = 0
#        updateSCT = 0
#        endT = 0
#        resT = 0
#        bufferSize = 5
#        withinLimit = 0
#        exceedLimit = 0
        start = time.time()
        while self.fitEval < self.MaxFit:
            self.fitEval = self.fitEval + self.dim
            #print 'SC', self.SC

            if init == False:
#                start = time.time()
                improveN, bestI, evalCount = self.genMeanBest(minimize)
#                genT = genT + time.time() - start
                init = True
            else :
#                start = time.time()
                improveN, bestI, evalCount = self.updateMeanBest(bestI,minimize)

            self.fitEval = self.fitEval + evalCount

#                updateT = updateT + time.time() - start
#            sBuffer = 0
#            sBufferNotChange = 0
#            for i in self.improveA:
#                if bestI != None :
#                    if (bestI in self.Inter) and (i not in self.Inter[bestI].arr):
#                        sBuffer = sBuffer + 1
#                        if takenBits[i] == False:
#                            sBufferNotChange = sBufferNotChange + 1
#                    elif bestI not in self.Inter:
#                        sBuffer = sBuffer + 1
#                        if takenBits[i] == False:
#                            sBufferNotChange = sBufferNotChange + 1
#            #print sBuffer, len(self.improveA)
#            if sBuffer <= bufferSize:
#               withinLimit = withinLimit + 1
#            else:
#               exceedLimit = exceedLimit + 1

            #print 
            #print 'Z', self.Z
            #print self.fitEval, self.SC
            #print improveN, bestI
            #lenImproveA.append( len(self.improveA) )
#            print self.oldindiv.bit
#            print 

            if improveN == False:
                if restart == True:
                    updateT = updateT + time.time() - start
                    startR = time.time()
                    oldbit = self.oldindiv.bit
                    self.oldindiv = self.evalPop(self.oldindiv)
                    self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.sumArr))
#                    fitGT = fitGT + time.time() - start  
                    
#                    start = time.time()
                    self.restart(fitName, minimize,False)
                    init = False
                    diff = self.diffBits(oldbit, self.oldindiv.bit)
                    #print oldbit, self.oldindiv.bit
#                    resT = resT + time.time() - start
                    #print 'restart'
                    takenBits = self.dim*[False]
                    
                    for i in diff:
#                        start = time.time()
                        self.update(i)
#                        updateCT = updateCT + time.time() - start

#                        start = time.time()
                        self.updateSC(i)
#                        updateSCT = updateSCT + time.time() - start 
                        self.updateWAS(i)
                        #print 'S', self.sumArr
                        #print 'Z', self.Z
                    initT = initT + time.time() - startR
                    start = time.time()
                else:
                    self.oldindiv = self.evalPop(self.oldindiv)
                    self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.sumArr))
                    self.fitEval = self.fitEval - 1
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
#                start = time.time()
                self.update(bestI)
#                updateCT = updateCT + time.time() - start
                #takenBits[bestI] = True

#                start = time.time()
                self.updateSC(bestI)
#                updateSCT = updateSCT + time.time() - start

                self.updateWAS(bestI)
                if self.oldindiv.bit[bestI] == '1':
                    self.oldindiv.bit[bestI] = '0'
                else:
                    self.oldindiv.bit[bestI] = '1'
#        print 'while loop', time.time() - startWhile

        
        #start = time.time()
        self.bsf = self.evalPop(self.bsf)
        self.fitEval = self.fitEval - 1
        diff = self.diffBits(self.bsf.bit, self.oldindiv.bit)
        for i in diff:
            self.update(i)
        self.bsf.fitG = self.bsf.fit - 2/float(self.dim) * (np.sum(self.sumArr))
        updateT = updateT + time.time() - start
        #print np.mean(lenImproveA), np.std(lenImproveA)
#        print 'genBest', genT
#        print 'updateBest', updateT
#        print 'fitGT', fitGT
#        print 'restartT', resT
#        print 'updateCT',updateCT
#        print 'updateSCT', updateSCT
#        print 'endT', endT
#        print 'allT', time.time() - startA
        #print 'Percentage of be withing buffer size limit', bufferSize, ':', withinLimit / (withinLimit + exceedLimit + 0.0)
        #print  withinLimit / (withinLimit + exceedLimit + 0.0)
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'fitG': self.bsf.fitG, 'bit':self.bsf.bit,'init':initT, 'update':updateT}

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

        self.WA = []
        compPSumT = 0
        updateT = 0
        while self.fitEval < self.MaxFit:
            # generate neighborhood and compute their fitness
            neighPop = self.neighWal()
            # compute the fitG (mean) of each neighborhood individuals
            improveN = False
            nCount = 0
            oldFit = self.oldindiv.fit
            oldFitG = self.oldindiv.fitG
            for n in neighPop:
                self.indiv = copy.deepcopy(n)

                self.indiv.fit = oldFit - 2*self.sumArr[nCount]
                self.fitEval = self.fitEval + 1

                start = time.time()
                self.indiv.fitG = oldFitG - 2*self.sumArr[nCount] + 4/float(self.dim) * self.compPhisum(nCount)
                compPSumT = compPSumT + time.time() - start
                #print 'neigh: ', self.indiv.bit, 'fit', self.indiv.fit, 'fitG', self.indiv.fitG
                if self.selectionFitNeigh(minimize) == True:
                    #print 'better neigh!'
                    improveN = True
                    changeBit = nCount

                nCount = nCount + 1

            if improveN == False:
                if restart == True:
                    oldbit = self.oldindiv.bit
                    self.restart(fitName, minimize)
                    newbit = self.oldindiv.bit
                    diff = self.diffBits(oldbit, newbit)
                    start = time.time()
                    for i in diff:
                        self.update(i)
                    updateT = updateT + time.time() - start
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                start = time.time()
                self.update(changeBit)
                updateT = updateT + time.time() - start

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
            #print 
            #print 'current', self.oldindiv.bit, 'fit', self.oldindiv.fit
            for i in neighs:
                self.indiv.bit = np.copy(i)
                self.indiv = self.evalPop(self.indiv)
                #print 'neigh: ', self.indiv.bit, 'fit', self.indiv.fit
                if  self.selectionFit(minimize) == True:
                    improveN = True

#            self.trace.append(Struct(fitEval= self.fitEval,fit = self.oldindiv.fit))
            if improveN == False:
                if restart == True:
                    self.restart('fit', minimize, True)
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
                    self.restart(fitName, minimize, True)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
        return { 'nEvals': self.fitEval, 'sol': self.bsf.fit, 'fitG': self.bsf.fitG, 'bit':self.bsf.bit}

    def diffBits(self, a, b):
        diff = []
        for i in range(self.dim):
            if a[i] != b[i]:
                diff.append(i)
        return diff

    def restart(self, fitName, minimize, evaluate):

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
            if evaluate == True:
                self.oldindiv = self.evalPop(self.oldindiv)
        else :
            self.oldindiv = self.initIndivNeigh(self.dim)
            if evaluate == True:
                self.oldindiv = self.evalPopNeigh(self.oldindiv, fitName, minimize)

#        print 'after restart'
#        print self.oldindiv.bit, self.oldindiv.fit, self.oldindiv.fitG
#        print self.bsf.bit, self.bsf.fit, self.bsf.fitG

    def walk(self, fitName, minimize, evaluate, length):
        # update the bsf solution
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

        flipBits = random.sample(xrange(self.dim), length)
        for i in flipBits:
            if self.oldindiv.bit[i] == '1':
                self.oldindiv.bit[i] = '0'
            else:
                self.oldindiv.bit[i] = '1'
        return flipBits


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
        self.WA = [] # array representing Walsh terms
        for k in self.model.w.keys(): 
            if self.model.w[k] != 0:
                a = []
                for i in range(self.dim):
                    if k[i] == '1':
                        a.append(i)
                self.WA.append( Struct(arr = a, w = self.model.w[k]) )
#        self.printWA()

    def binCount(self, arr, bit):
        """
        count the one bit of union self.WA[i].arr and bit
        """
        s = 0
        for i in arr:
            if bit[i] == '1':
                s = s + 1
        return s

    def genFitBest(self,minimize):
        """
        generate the index of best neigh according to sumArr only (surrogate of fitness)
        """
        # check improving move
        improve = False
        self.improveA = []
        for i in range(self.dim):
            if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i]<-self.threshold):
                self.improveA.append(i) 
                improve = True

        if improve == False:
            return False, None, self.dim

        for i in self.improveA:
            if i == self.improveA[0]:
                best = self.sumArr[i]
                bestI = i
            elif (best<self.sumArr[i] - self.threshold and minimize == True) or (best>self.sumArr[i] + self.threshold and minimize == False): # seek for max S
                best = self.sumArr[i]
                bestI = i
                    
        return True, bestI, self.dim

    def genFitBest2(self,minimize):
        """
        generate the index of best distance 2 neighborhoods according to sumArr only (surrogate of fitness)

        return: 1) whether there is an improving move in distance 2 neigh
                2) the index of best distance 1 neigh for taking the next move
                3) the number of evaluations consumed by this step
        """
        evalCount = self.dim
        self.improveA = []
        neighImprove = []

        # checking the distance 1 neigh
        for i in range(self.dim):
            if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i]<-self.threshold):
                self.improveA.append(i) 
                neighImprove.append(Struct(index =[i], val = self.sumArr[i])) # add dist 1 neigh into consideration as well

#        print self.sumArr
#        print self.improveA

        # checking the distance 2 neigh, remember to preserve context
        for i in range(self.dim):
            self.mirrorParam() # everything is reset, pretending nothing happened
            self.updateFake(i)
            if i in self.improveAfake:
                self.improveAfake.remove(i)
            if i in self.Inter:
                for j in self.Inter[i].arr:
                    evalCount = evalCount + 1

                    self.sumArrFake[j] = self.sumArrFake[j]+self.sumArr[i]
                    if (minimize == True and self.sumArrFake[j] > self.threshold) or (minimize == False and self.sumArrFake[j]<-self.threshold):
                        if j not in self.improveAfake:
                            self.improveAfake.append(j)
                    elif j in self.improveAfake:
                        self.improveAfake.remove(j)

            for j in self.improveAfake:
                neighImprove.append(Struct(index =[i,j], val = self.sumArrFake[j]))

#        for i in range(len(neighImprove)):
#            print neighImprove[i].index, neighImprove[i].val

#        for i in neighImprove:
#            print i.index, i.val
               
        if not neighImprove:
            #return False, None, self.dim*self.dim
            return False, None, evalCount

        for i in range(len(neighImprove)):
            if i == 0:
                best = neighImprove[i].val
                bestI = neighImprove[i].index
            elif ( best<neighImprove[i].val - self.threshold and minimize == True) or ( best>neighImprove[i].val + self.threshold and minimize == False ): # seek for max S
                best = neighImprove[i].val
                bestI = neighImprove[i].index

        bestIlist = []
        for i in range(len(neighImprove)):
            if abs(best - neighImprove[i].val) < self.threshold:
                candI = neighImprove[i].index
                if candI not in bestIlist:
                    bestIlist.append(candI)

        #print 'bestIlist',bestIlist
        bestI = random.choice(bestIlist)
#        print 'bestList', bestIlist
#        print 'bestI', bestI
        if type(bestI) is int:
            # make a consistent interface
            bestI = [bestI]
                    
        #return True, bestI, self.dim*self.dim
        return True, bestI, evalCount

    def genFitBestsm(self,minimize):
        """
        generate the index of best neigh according to sumArr only (surrogate of fitness)
        """
        #TODO need to update the threshold
        # check improving move
        improve = False
        self.Buffer = []
        for i in range(self.dim):
            if (minimize == True and self.sumArr[i] > 0) or (minimize == False and self.sumArr[i]<0):
                self.Buffer.append(i) 
                improve = True

        if improve == False:
            return False, None

        for i in self.Buffer:
            if i == self.Buffer[0]:
                best = self.sumArr[i]
                bestI = i
            elif (best<self.sumArr[i] and minimize == True) or (best>self.sumArr[i] and minimize == False): # seek for max S
                best = self.sumArr[i]
                bestI = i

        self.P = [bestI]

        # produce buffer list (the independent improving set)
        for i in [ j for j in self.Buffer if j != bestI ]:
            if i not in self.Inter:
                self.P.append(i)
            else :
                inter = False
                for j in self.P:
                    if j in self.Inter[i].arr:
                        inter = True
                if inter == False:
                    self.P.append(i)

        return True, bestI

    def updateFitBestsm(self, minimize):
        #TODO need to update the threshold
        for i in self.P:
            self.Buffer.remove(i)
            if i in self.Inter:
                for j in self.Inter[i].arr:
                    if ((minimize == True and self.sumArr[j] > 0) or (minimize == False and self.sumArr[j]<0)): 
                        if j not in self.Buffer:
                            self.Buffer.append(j)
                    elif j in self.Buffer:
                        self.Buffer.remove(j)

        if not self.Buffer:
            return False, None
        #print improveA

        for i in self.Buffer:
            if i == self.Buffer[0]:
                best = self.sumArr[i]
                bestI = i
            elif (best<self.sumArr[i] and minimize == True) or (best>self.sumArr[i] and minimize == False): # seek for max S
                best = self.sumArr[i]
                bestI = i

        self.P = [bestI]

        for i in [ j for j in self.Buffer if j != bestI]:
            if i not in self.Inter:
                self.P.append(i)
            else :
                inter = False
                for j in self.P:
                    if j in self.Inter[i].arr:
                        inter = True
                if inter == False:
                    self.P.append(i)

        return True, bestI

    def updateFitBest(self, p, minimize):
        self.improveA.remove(p)
        evalCount = 0
        if p in self.Inter:
            for i in self.Inter[p].arr: 
                evalCount = evalCount + 1
                if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i]< self.threshold ):
                    if i not in self.improveA:
                        self.improveA.append(i)
                elif i in self.improveA:
                    self.improveA.remove(i)

        if not self.improveA:
            return False, None, evalCount

        for i in self.improveA:
            if i == self.improveA[0]:
                best = self.sumArr[i]
                bestI = i
            elif (best<self.sumArr[i] and minimize == True) or (best>self.sumArr[i] and minimize == False): 
                best = self.sumArr[i]
                bestI = i
                    
        return True, bestI, evalCount

    def updateFitBest2(self, P, minimize):
        """ 
        generate the index of best distance 2 neighborhoods according to sumArr only (surrogate of fitness), by performing partial updates
        """
        neighImprove = []
        evalCount = 0
        
        for p in P:
            if p in self.improveA:
                self.improveA.remove(p)
            if p in self.Inter:
                for i in self.Inter[p].arr: 
                    evalCount = evalCount + 1
                    if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i]<-self.threshold):
                        if i not in self.improveA:
                            self.improveA.append(i)
                    elif i in self.improveA:
                        self.improveA.remove(i)
        p = P[-1]

        for i in self.improveA:
            """ add distance 1 neigh under consideration """
            neighImprove.append(Struct(index=[i], val = self.sumArr[i]))

        # checking the distance 2 neigh, remember to preserve context
        for i in [k for k in range(self.dim) if k!=p]:
            self.mirrorParam() # everything is reset, pretending nothing happened
            self.updateFake(i)
            if i in self.improveAfake:
                self.improveAfake.remove(i)
            if i in self.Inter:
                for j in self.Inter[i].arr:
                    evalCount = evalCount + 1
                    self.sumArrFake[j] = self.sumArrFake[j]+self.sumArr[i]
                    if (minimize == True and self.sumArrFake[j] > self.threshold) or (minimize == False and self.sumArrFake[j]<-self.threshold):
                        if j not in self.improveAfake:
                            self.improveAfake.append(j)
                    elif j in self.improveAfake:
                        self.improveAfake.remove(j)

            for j in self.improveAfake:
                neighImprove.append(Struct(index =[i,j], val = self.sumArrFake[j]))

#        for i in range(len(neighImprove)):
#            print neighImprove[i].index, neighImprove[i].val
        
        if not neighImprove:
            #return False, None, evalCount + self.dim * (self.dim-1)
            return False, None, evalCount 

        for i in range(len(neighImprove)):
            if i == 0:
                best = neighImprove[i].val
                bestI = neighImprove[i].index[0]
            elif (best<neighImprove[i].val - self.threshold and minimize == True) or (best>neighImprove[i].val + self.threshold and minimize == False): # seek for max S
                best = neighImprove[i].val
                bestI = neighImprove[i].index[0]

        bestIlist = []
        for i in range(len(neighImprove)):
            if abs(best - neighImprove[i].val) < self.threshold:
                candI = neighImprove[i].index
                if candI not in bestIlist:
                    bestIlist.append(candI)

        bestI = random.choice(bestIlist)
        if type(bestI) is int:
            # make a consistent interface
            bestI = [bestI]

        #return True, bestI, evalCount + self.dim * (self.dim-1)
        return True, bestI, evalCount 

    def genMeanBest(self,minimize):
        """
        generate the index of best neigh according to {S_p(X)-2/N \Sigma_{i=1}^{N}C_{ip}(X)} only (surrogate of fitness)
        """
        # check improving move 
        improve = False
        self.improveA = []
        for i in range(self.dim):
            if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i]<-self.threshold):
                self.improveA.append(i)
                improve = True

        if improve == False:
            return False, None, self.dim

        # find the best value
        for i in self.improveA:
            if i == self.improveA[0]:
                best = self.SC[i]
                bestI = i
            elif ( best<self.SC[i] - self.threshold and minimize == True ) or ( best>self.SC[i] + self.threshold and minimize == False ): # seek for max S
                best = self.SC[i]
                bestI = i
        return True, bestI, self.dim

    def genMeanBest2(self,minimize):
        """
        generate the index of best distance 2 neigh according to {S_p(X)-2/N \Sigma_{i=1}^{N}C_{ip}(X)} only (surrogate of fitness)

        return: 1) whether there is an improving move in distance 2 neigh
                2) the index of best distance 1 neigh for taking the next move
                3) the number of evaluations consumed by this step
        """
        # check improving move 
        evalCount = self.dim
        self.improveA = []
        neighImprove = []

        # check dist 1 neigh
        for i in range(self.dim):
            if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i]<-self.threshold):
                self.improveA.append(i)
                neighImprove.append(Struct(index =[i], val = self.SC[i])) # add dist 1 neigh into consideration as well

        # checking the distance 2 neigh, remember to preserve context
        for i in range(self.dim):
            self.mirrorSCParam() # everything is reset, pretending nothing happened
            self.updateFake(i)
            self.updateSCfake(i)
            if i in self.improveAfake:
                self.improveAfake.remove(i)
            if i in self.Inter:
                for j in self.Inter[i].arr:
                    evalCount = evalCount + 1
                    self.SCfake[j] = self.SCfake[j]+self.SC[i]
                    if (minimize == True and self.SCfake[j] > self.threshold) or (minimize == False and self.SCfake[j]<-self.threshold):
                        if j not in self.improveAfake:
                            self.improveAfake.append(j)
                    elif j in self.improveAfake:
                        self.improveAfake.remove(j)

            for j in self.improveAfake:
                neighImprove.append(Struct(index =[i,j], val = self.SCfake[j]))

        if not neighImprove:
            #return False, None, self.dim*self.dim
            return False, None, evalCount

        for i in range(len(neighImprove)):
            if i == 0:
                best = neighImprove[i].val
                bestI = neighImprove[i].index
            elif ( best<neighImprove[i].val - self.threshold and minimize == True) or ( best>neighImprove[i].val + self.threshold and minimize == False ): # seek for max S
                best = neighImprove[i].val
                bestI = neighImprove[i].index

        bestIlist = []
        for i in range(len(neighImprove)):
            if abs(best - neighImprove[i].val) < self.threshold:
                candI = neighImprove[i].index
                if candI not in bestIlist:
                    bestIlist.append(candI)

        #print 'bestIlist',bestIlist
        bestI = random.choice(bestIlist)
#        print 'bestList', bestIlist
#        print 'bestI', bestI
        if type(bestI) is int:
            # make a consistent interface
            bestI = [bestI]
                    
        #return True, bestI, self.dim*self.dim
        return True, bestI, evalCount

    def updateMeanBest(self, p, minimize):
        self.improveA.remove(p)
        evalCount = 0
        if p in self.Inter:
            for i in self.Inter[p].arr:
                evalCount = evalCount + 1
                if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i]<-self.threshold):
                    if i not in self.improveA:
                        self.improveA.append(i)
                elif i in self.improveA:
                    self.improveA.remove(i)

        if not self.improveA:
            return False, None, evalCount

        for i in self.improveA:
            if i == self.improveA[0]:
                best = self.SC[i]
                bestI = i
            elif (best<self.SC[i] - self.threshold and minimize == True) or (best>self.SC[i]+self.threshold and minimize == False): # seek for max S
                best = self.SC[i]
                bestI = i
                    
        return True, bestI, evalCount

    def updateMeanBest2(self, P, minimize):
        evalCount = 0
        neighImprove = []

        for p in P:
            if p in self.improveA:
                self.improveA.remove(p)
            if p in self.Inter:
                for i in self.Inter[p].arr: 
                    evalCount = evalCount + 1
                    if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i]<-self.threshold):
                        if i not in self.improveA:
                            self.improveA.append(i)
                    elif i in self.improveA:
                        self.improveA.remove(i)
        p = P[-1]

        for i in self.improveA:
            """ add distance 1 neigh under consideration """
            neighImprove.append(Struct(index=[i], val = self.SC[i]))

        # checking the distance 2 neigh, remember to preserve context
        for i in [k for k in range(self.dim) if k!=p]:
            self.mirrorSCParam() # everything is reset, pretending nothing happened
            self.updateFake(i)
            self.updateSCfake(i)
            if i in self.improveAfake:
                self.improveAfake.remove(i)
            if i in self.Inter:
                for j in self.Inter[i].arr:
                    evalCount = evalCount + 1
                    self.SCfake[j] = self.SCfake[j]+self.SC[i]
                    if (minimize == True and self.SCfake[j] > self.threshold) or (minimize == False and self.SCfake[j]<-self.threshold):
                        if j not in self.improveAfake:
                            self.improveAfake.append(j)
                    elif j in self.improveAfake:
                        self.improveAfake.remove(j)

            for j in self.improveAfake:
                neighImprove.append(Struct(index =[i,j], val = self.SCfake[j]))

        if not neighImprove:
            #return False, None, self.dim*self.dim
            return False, None, evalCount

        for i in range(len(neighImprove)):
            if i == 0:
                best = neighImprove[i].val
                bestI = neighImprove[i].index
            elif ( best<neighImprove[i].val - self.threshold and minimize == True) or ( best>neighImprove[i].val + self.threshold and minimize == False ): # seek for max S
                best = neighImprove[i].val
                bestI = neighImprove[i].index

        bestIlist = []
        for i in range(len(neighImprove)):
            if abs(best - neighImprove[i].val) < self.threshold:
                candI = neighImprove[i].index
                if candI not in bestIlist:
                    bestIlist.append(candI)

        #print 'bestIlist',bestIlist
        bestI = random.choice(bestIlist)
#        print 'bestList', bestIlist
#        print 'bestI', bestI
        if type(bestI) is int:
            # make a consistent interface
            bestI = [bestI]
                    
        #return True, bestI, self.dim*self.dim
        return True, bestI, evalCount


    def initWal(self):
        """ 
        1. 
        compute the sum array for the first time, according to the initial solution
        2. 
        initialize a full matrix, Coincidence, C_ij, which is a symmetric one, and diagonal is empty (as we don't want to recompute SumArr), empty whenever i >= j.
        Two options for implementing C matrix
            a*. precompute C matrix 
            b. construct it on the fly using dict()
        3. 
        initialize interaction list for each variable
        4.
        initialize a dict of interaction structure, where interactive bits and the index of WAS (walsh including sign)
        """
        self.sumArr = np.zeros(self.dim)
        self.WAS = np.tile(Struct(arr = [], w = 0), len(self.model.w.keys()))# Walsh coefficients with sign, represented in Array
        self.lookup = dict()
        self.infectBit = dict()
        self.C = np.zeros((self.dim,self.dim)) # coincidence matrix
        self.Inter = dict()
        self.InterCount = np.zeros(self.dim)

#        for i in range(self.dim):
#            self.Inter.append(Set()) # coincidence matrix

        #self.C = dict() # coincidence matrix
        for i in range(len(self.WA)):
            W = int(math.pow(-1, self.binCount(self.WA[i].arr, self.indiv.bit))) * self.WA[i].w
            self.WAS[i] = Struct(arr = self.WA[i].arr, w = W)
            comb = self.genComb(len(self.WA[i].arr))
            #print i, self.WA[i].arr, comb

            for j in self.WA[i].arr:
                self.sumArr[j] = self.sumArr[j] + W
                if len(self.WA[i].arr)>1: # for at least order Walsh terms
                    if j not in self.Inter: # the entry of i doesn't exist yet
                        self.Inter[j] = Struct(arr=Set(), WI=Set())

                    for k in self.WA[i].arr:
                        if k != j:
                            self.Inter[j].arr.add(k)
                    self.Inter[j].WI.add(i)

                # add list of order >= 3 Walsh terms for the purpose of updating C matrix
                if len(self.WA[i].arr) >= 3:
                    if j not in self.infectBit: 
                        self.infectBit[j] = [Struct(arr=self.WA[i].arr, WI=i)]
                    else :
                        self.infectBit[j].append(Struct(arr=self.WA[i].arr, WI=i))

            for j in comb: # for each list in comb
                j0 = self.WA[i].arr[int(j[0])]
                j1 = self.WA[i].arr[int(j[1])]
#                if (j0, j1) in self.C.keys():
#                    self.C[j0,j1] = self.C[j0,j1] + W
#                elif W != 0:
#                    self.C[j0,j1] = W
                self.C[j0,j1] = self.C[j0,j1] + W
        
        for i in range(self.dim):
            if i in self.Inter:
                self.InterCount[i] = len(self.Inter[i].arr)
            else :
                self.InterCount[i] = 0

#        print self.Inter
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

    def initSC(self):
        # compute the SC array
        self.SC = np.zeros(self.dim)
        self.Z = np.zeros(self.dim)
        self.orderC = np.zeros((self.dim,self.dim))

        for p in range(self.dim):

            phi = np.zeros(self.model.k+1)
            if p in self.Inter:
                for i in self.Inter[p].WI:
                    order = len(self.WAS[i].arr)
                    phi[order-1] = phi[order-1] + self.WAS[i].w

            self.Z[p] = self.sumArr[p]
            for i in range(1, self.model.k+1):
                if phi[i] != 0:
                    self.Z[p] = self.Z[p] + i * phi[i]

            self.SC[p] = self.sumArr[p] - 2/float(self.dim) * self.Z[p]
            #self.SC[i] = self.sumArr[i] - 2/float(self.dim) * self.compCsum(i)

        for i in range(len(self.WAS)):
            lenArr = len(self.WAS[i].arr)
            comb = self.genComb(lenArr)
            for j in comb:
                j0 = self.WAS[i].arr[int(j[0])]
                j1 = self.WAS[i].arr[int(j[1])]
                self.orderC[j0,j1] = self.orderC[j0,j1] + lenArr * self.WAS[i].w

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

    def compPhisum(self,p):
        """
        \varphi_{p,i}^{\prime}(x) = \Sigma_{order j terms, that touches bit p}
        """
        phi = np.zeros(self.model.k+1)
        if p in self.Inter:
            for i in self.Inter[p].WI:
                order = len(self.WAS[i].arr)
                phi[order-1] = phi[order-1] + self.WAS[i].w

        s = self.sumArr[p]
        for i in range(1, self.model.k+1):
            if phi[i] != 0:
                s = s + i * phi[i]
        return s

    def compCsum(self,p):
        """
        \sigma_{i=1}^{N} C_{ip}: be careful with C_{ii}, i \in N
        """
        s = 0

        for i in range(p):
            s = s + self.C[i,p]

        for i in range(p+1, self.dim):
            s = s + self.C[p,i]

        s = s + self.sumArr[p]

        return s

    def update(self, p):
        """
        By keeping track of coincidence matrix, 
        Cij stands for S_i(y_j) = S_i(x) - C_ij
        partially update the Sum Array and self.WAS, given the bit which is changed
        """
        self.sumArr[p] = - self.sumArr[p]
        
        if p in self.Inter:
            for i in self.Inter[p].arr:
                if i < p:
                    self.sumArr[i] = self.sumArr[i] - 2*self.C[i,p]
                    self.C[i,p] = - self.C[i,p]
                else:
                    self.sumArr[i] = self.sumArr[i] - 2*self.C[p,i]
                    self.C[p,i] = - self.C[p,i]

        # update the rest of elements in C matrix
        if p in self.infectBit.keys():
            for i in self.infectBit[p]:
                arr = copy.deepcopy(i.arr)
                arr.remove(p)
                comb = self.genComb(len(arr))
                for k in range(len(comb)):
                    k0 = arr[int(comb[k][0])]
                    k1 = arr[int(comb[k][1])]
                    self.C[k0,k1] = self.C[k0,k1] - 2 * self.WAS[i.WI].w

    def updateFake(self, p):
        """
        The fake version, the updates are made to the mirror data structures
        By keeping track of coincidence matrix, 
        Cij stands for S_i(y_j) = S_i(x) - C_ij
        partially update the Sum Array and self.WAS, given the bit which is changed
        """
        self.sumArrFake[p] = - self.sumArrFake[p]
        
        if p in self.Inter:
            for i in self.Inter[p].arr:
                if i < p:
                    self.sumArrFake[i] = self.sumArrFake[i] - 2*self.C[i,p]
                else:
                    self.sumArrFake[i] = self.sumArrFake[i] - 2*self.C[p,i]

    def updateWAS(self,p):
        if p in self.Inter:
            for i in self.Inter[p].WI:
                self.WAS[i].w = - self.WAS[i].w

    def updateSC(self, p):
        self.SC[p] = - self.SC[p]
        self.Z[p] = - self.Z[p]

        #update Z array
        if p in self.Inter:
            for i in self.Inter[p].arr:
                if i < p :
                    self.Z[i] = self.Z[i]  - 2* self.orderC[i,p]
                    self.orderC[i,p] = - self.orderC[i,p]
                else :
                    self.Z[i] = self.Z[i]  - 2* self.orderC[p,i]
                    self.orderC[p,i] = - self.orderC[p,i]
                self.SC[i] = self.sumArr[i] - 2/float(self.dim) * self.Z[i]

        if p in self.infectBit.keys():
            for i in self.infectBit[p]:
                arr = copy.deepcopy(i.arr)
                arr.remove(p)
                lenArr = len(arr)
                comb = self.genComb(lenArr)
                for k in range(len(comb)):
                    k0 = arr[int(comb[k][0])]
                    k1 = arr[int(comb[k][1])]
                    self.orderC[k0,k1] = self.orderC[k0,k1] - 2 * (lenArr + 1)* self.WAS[i.WI].w

    def updateSCfake(self, p):
        self.SCfake[p] = - self.SCfake[p]
        self.Zfake[p] = - self.Zfake[p]

        #update Z array
        if p in self.Inter:
            for i in self.Inter[p].arr:
                if i < p :
                    self.Zfake[i] = self.Zfake[i]  - 2* self.orderC[i,p]
                else :
                    self.Zfake[i] = self.Zfake[i]  - 2* self.orderC[p,i]
                self.SCfake[i] = self.sumArrFake[i] - 2/float(self.dim) * self.Zfake[i]

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

    def printInter(self):
        for i in self.Inter:
            print i,self.Inter[i].arr

    def mirrorParam(self):
        """ create a copy of all data structures required for update """
        self.sumArrFake = copy.deepcopy(self.sumArr)
        self.improveAfake = copy.deepcopy(self.improveA)
#        self.Cfake = copy.deepcopy(self.C)
#        self.WASfake = copy.deepcopy(self.WAS)

    def mirrorSCParam(self):
        """ create a copy of all data structures required for SC update """
        self.sumArrFake = copy.deepcopy(self.sumArr)
        self.Zfake = copy.deepcopy(self.Z)
        self.SCfake = copy.deepcopy(self.SC)
        self.improveAfake = copy.deepcopy(self.improveA)

    def bitImpact(self):
        self.transWal()
        self.indiv = self.initIndiv(self.dim)
        self.fitEval = 0
        self.initWal()
        rep = 100
        fitChange = np.zeros(self.dim)
        for i in range(rep):
            self.indiv = self.initIndiv(self.dim)
            self.indiv = self.evalPop(self.indiv)
#            print 'fx', self.indiv.bit, self.indiv.fit
#            print
            fx = self.indiv.fit
            for j in range(self.dim):
                self.newIndiv = copy.deepcopy(self.indiv)
                if self.newIndiv.bit[j] == '1':
                    self.newIndiv.bit[j] = '0'
                else:
                    self.newIndiv.bit[j] = '1'
                self.newIndiv = self.evalPop(self.newIndiv)
#                print 'fy', self.newIndiv.bit, self.newIndiv.fit
                fy = self.newIndiv.fit
                fitChange[j] = fitChange[j] + np.abs(fx-fy)
#                print fitChange
#                print 
            fitChange[j] = fitChange[j]/float(rep)


#        for i in zip(fitChange,self.InterCount):
#            print i
        print "%.2f" % (np.corrcoef([fitChange,self.InterCount])[1,0])
        return {'nEvals': None, 'sol': None, 'bit': None, 'init': None, 'update': None}
