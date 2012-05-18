""" 
Local Search module:
    This local search module gets the problem (in terms of objective function),
    and return a solution found by hill climber
"""

import numpy as np
import WalshAnalysis as wal
import itertools as it
import nkLandscape as nk
import tool as tl
import random
import string
import math
import copy
import sys
import os
import time
import pdb
from sets import Set
from operator import itemgetter
import individual 

# import for Cython
from cpython cimport bool
cimport cython
from cython.parallel import prange, parallel, threadid
from libcpp.vector cimport vector
from libcpp.set cimport set
from libc.stdlib cimport malloc, free
from cython.operator cimport dereference as deref, preincrement as inc

class Struct:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

cdef class LocalSearch:
    cdef object func
    cdef object model
    cdef int MaxFit
    cdef int dim
    cdef int fitEval
    cdef list oldpop
    cdef object bsf

    def __init__(self, model, MaxFit, dim):
        self.func = model.compFit
        self.model = model
        self.MaxFit = MaxFit
        self.dim = dim

    cpdef run(self, fitName, minimize, restart, compM, beamWidth=1):
        if compM == 'bf':
            if fitName =='fit': 
                return self.runFit(minimize,restart)
            else:
                return self.runNeigh(fitName, minimize,restart)
        elif compM == 'walWalk':
            if fitName == 'fit':
                return self.runFitSwalk(fitName, minimize, restart)
            elif fitName == 'mean':
                return self.runMeanSCwalk(fitName, minimize, restart)
        elif compM == 'walWalkNext':
            if fitName == 'fit':
                return self.runFitSwalkNext(fitName, minimize, restart)
            elif fitName == 'mean':
                return self.runMeanSCwalkNext(fitName, minimize, restart)
        elif compM == 'walRest':
            if fitName == 'fit':
                return self.runFitSrest(fitName, minimize, restart)
            elif fitName == 'mean':
                return self.runMeanSCrest(fitName, minimize, restart)
        elif compM == 'walRestNext':
            if fitName == 'fit':
                return self.runFitSrestNext(fitName, minimize, restart)
            elif fitName == 'mean':
                return self.runMeanSCrestNext(fitName, minimize, restart)
        elif compM == 'BeamWalkNext':
            if fitName == 'fit':
                return self.runBeamFitSwalkNext(fitName, minimize, restart, beamWidth)
            elif fitName == 'mean':
                return self.runBeamMeanSCwalkNext(fitName, minimize, restart, beamWidth)
        elif compM == 'BeamWalk':
            if fitName == 'fit':
                return self.runBeamFitSwalk(fitName, minimize, restart, beamWidth)
            elif fitName == 'mean':
                return self.runBeamMeanSCwalk(fitName, minimize, restart, beamWidth)
        elif compM == 'supm':
            if fitName == 'fit':
                return self.runFitsm(fitName, minimize, restart)
            elif fitName == 'mean':
                return self.runMeansm(fitName, minimize, restart)
        elif compM == 'bitImp':
            return self.bitImpact()
        elif compM == 'walSearch':
            return self.runWalSearch(fitName, minimize, restart)
        elif compM == 'checkOptWal':
            self.checkOptWal()
            return None
        elif compM == 'checkHyper':
            return self.checkHyper()
        elif compM == 'checkHyperRank':
            return self.checkHyperRank()
        elif compM == 'hyperSearch' or compM == 'hyperSqSearch' or compM == 'hyperWalSearch':
            if fitName == 'fit':
                return self.hyperSearchFit(fitName, minimize, restart)
            elif fitName == 'mean':
                return self.hyperSearchMean(fitName, minimize, restart)


    def checkHyperRank(self):
        """
        examine the rank of optimum hyperplane in those list of hyperplanes associated with each subfunction
        """
        self.model.transWal()
        bit,fit = tl.compFit(self.model)
        a = sorted(zip(bit,fit), key=lambda a_entry: a_entry[1]) 
        optBit = a[0][0]
        optFit = a[0][1]
        print 'opti\n',optBit, optFit
        rank = 0

        for i in range(self.dim):
            subBit = self.model.neighs[i][:]
            subBit.append(i)
            subBit.sort()

            optTpl = []
            for j in subBit:
                if optBit[j] == '1':
                    optTpl.append(j)

            # check every template that matches the subfunction
            seqBits = nk.genSeqBits(len(subBit))
            print 
            schFitArr = []
            for j in seqBits:
                schFit = 0
                # convert bit string to array representation
                schTpl = []
                for k in range(len(j)):
                    if j[k] == '1':
                        schTpl.append(subBit[k])

                for k in self.model.WA:
                    subset = True
                    for l in k.arr:
                        if l not in subBit:
                            subset = False
                            break
                    if subset == True:
                        schFit = schFit + int(math.pow(-1, self.binCountArr(k.arr, schTpl))) * k.w
                schFitArr.append(Struct(fit=schFit,arr=schTpl))
                #print subBit, j, schFit
            print subBit

            schFitSort = sorted(schFitArr, key=lambda i: abs(i.fit))
            # check the rank of optimum solution in the  list of hyperplane associated with a subfunction
            for j in range(len(schFitSort)):
                if schFitSort[j].arr == optTpl:
                    rank = rank + j
                    print j

        print 'rank',rank

    def checkHyper(self):
        """
        examine the fitness of one particular bit over all hyperplanes associated with each subfunction
        """
        bit,fit = tl.compFit(self.model)
        a = sorted(zip(bit,fit), key=lambda a_entry: a_entry[1]) 
        optBit = a[0][0]
        optFit = a[0][1]
        print 'opti\n',optBit, optFit
        #for i in range(len(a)): 
#        for i in range(10): 
#            print '%s\t%.3f' %(a[i][0],a[i][1])
        
        rep = 10
        for i in range(rep):
            sol = self.genSolProp(self.model.sumFitA)
            hamDist = 0
            # compute the hamming distance
            for i in range(self.dim):
                if sol[i] != optBit[i]:
                    hamDist = hamDist + 1
            print 'Hyper solution\n', sol, self.func(sol), hamDist

        randSol = self.initIndiv(self.dim)
        hamDistRand = 0
        for i in range(self.dim):
            if randSol.bit[i] != optBit[i]:
                hamDistRand = hamDistRand + 1
        print 'Random Solution\t', self.func(randSol.bit), hamDistRand

        return {'nEvals': 0, 'sol': self.func(sol), 'bit': hamDist, 'init': self.func(randSol.bit), 'update': hamDistRand}
       


    def genSolProp(self, sumFitA):
#        for i in range(self.dim):
#            print '%d\tOne: %.2f\tZero: %.2f' %(i, sumFitA[i].one, sumFitA[i].zero)
#        print 
        sol = []
        for i in range(self.dim):
            if abs(sumFitA[i].zero) < self.threshold and abs(sumFitA[i].one) < self.threshold:
                if random.random() < 0.5: 
                    sol.append('0')
                else:
                    sol.append('1')
            else:
                if random.random() < sumFitA[i].zero / (sumFitA[i].one + sumFitA[i].zero + 0.0):
                    sol.append('0')
                else:
                    sol.append('1')
        return sol
        
##     def hyperSearchFit(self,fitName, minimize, restart):
##         """ 
##         performing hyper search using the probability generated by Hyperplane analysis
##         """

##         self.fitEval = 0
##         start = os.times()[0]
##         self.oldindiv = Struct( fit = 0, bit = self.genSolProp(self.model.sumFitA) )
## #        self.oldindiv = self.initIndiv(self.dim)
##         self.oldindiv = self.evalPop(self.oldindiv)
##         self.initWal(self.model)

##         self.bsf = copy.deepcopy(self.oldindiv)
##         #self.model.WA = []

##         init = False
##         updateT = 0
##         walkLen = 10
##         initT = os.times()[0] - start
##         start = os.times()[0]
##         while self.fitEval < self.MaxFit:
##             if init == False:
##                 improveN, bestI, evalCount = self.genFitBest(minimize)
##                 init = True
##             else:
##                 improveN, bestI, evalCount = self.oldindiv.updateFitBest(bestI,minimize)
##             self.fitEval = self.fitEval + evalCount

##             if improveN == False:
##                 if restart == True:
##                     updateT = updateT + os.times()[0] - start
##                     startR = os.times()[0]
##                     self.oldindiv = self.evalPop(self.oldindiv)
## #                    oldbit = self.oldindiv.bit
## #                    self.fitEval = self.fitEval - 1
## #                    self.hyperRestart(fitName, minimize, False)
## #                    diff = self.diffBits(oldbit, self.oldindiv.bit)

##                     diff, self.oldindiv = self.walk(fitName, minimize,False, walkLen, self.oldindiv)
##                     init = False
##                     for i in diff:
##                         self.oldindiv.update(i)
##                         self.oldindiv.updateWAS(i)
##                     initT = initT + os.times()[0] - startR
##                     start = os.times()[0]
##                 else:
##                     self.oldindiv = self.evalPop(self.oldindiv)
##                     return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit}
##             else : # improveN is TRUE 
##                 self.oldindiv.update(bestI)
##                 self.oldindiv.updateWAS(bestI)
##                 if self.oldindiv.bit[bestI] == '1':
##                     self.oldindiv.bit[bestI] = '0'
##                 else:
##                     self.oldindiv.bit[bestI] = '1'
##         self.bsf = self.evalPop(self.bsf)
##         updateT = updateT + os.times()[0] - start
##         return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit,'init':initT, 'update':updateT} 

    def checkOptWal(self):
        """
        check the sorted Walsh signs for all Walsh coefficients
        """
        bit,fit = tl.compFit(self.model)
        a = sorted(zip(bit,fit), key=lambda a_entry: a_entry[1]) 
        self.model.transWal()
        
        fit = []
        allCount = []
        for i in range(2):
            optBit = a[i][0]
            optFit = a[i][1]
            print 'Top', i+1, 'solution', optBit, optFit
            self.WAsort = sorted(self.model.WA, key=lambda i: abs(i.w), reverse=True)
            WAS = []
            negCount = 0
            for i in self.WAsort:
                temp = int(math.pow(-1, self.binCount(i.arr, optBit))) * i.w
                if temp  < -self.threshold:
                    negCount = negCount + 1
                WAS.append(int(math.pow(-1, self.binCount(i.arr, optBit))) * i.w)

            for i in range(50):
                if i != 0:
                    print '%.4f' %(WAS[i]), self.WAsort[i].arr

            #print 'negative / All non-zero:\n%d/%d\n' %(negCount, len(WAS)-1)
            fit.append(optFit)
            allCount.append(negCount)

        print 'Correlation: fitness / negetive Count', (np.corrcoef([fit,allCount])[1,0])
#        i = len(a)/2
#        optBit = a[i][0]
#        optFit = a[i][1]
#        print 'Ave solution', optBit, optFit
#        self.WAsort = sorted(self.WA, key=lambda i: abs(i.w), reverse=True)
#        WAS = []
#        negCount = 0
#        for i in self.WAsort:
#            temp = int(math.pow(-1, self.binCount(i.arr, optBit))) * i.w
#            if temp  < -self.threshold:
#                negCount = negCount + 1
#            WAS.append(int(math.pow(-1, self.binCount(i.arr, optBit))) * i.w)
#
##        for i in range(len(WAS)):
##            if i != 0:
##                print '%.4f' %(WAS[i])
#
#        print 'Negative / All non-zero:\n%d/%d\n' %(negCount, len(WAS)-1)
#
#        i = len(a) - 1
#        optBit = a[i][0]
#        optFit = a[i][1]
#        print 'Worst solution', optBit, optFit
#        self.WAsort = sorted(self.model.WA, key=lambda i: abs(i.w), reverse=True)
#        WAS = []
#        negCount = 0
#        for i in self.WAsort:
#            temp = int(math.pow(-1, self.binCount(i.arr, optBit))) * i.w
#            if temp  < -self.threshold:
#                negCount = negCount + 1
#            WAS.append(int(math.pow(-1, self.binCount(i.arr, optBit))) * i.w)
#
##        for i in range(len(WAS)):
##            if i != 0:
##                print '%.4f' %(WAS[i])
#        print 'Negative / All non-zero:\n%d/%d\n' %(negCount, len(WAS)-1)

    def bitImpact(self):
        self.model.transWal()
        self.indiv = self.initIndiv(self.dim)
        self.fitEval = 0
        self.initWal(self.model)
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

##     def runFitsm(self,fitName, minimize, restart):
##         """ 
##         steepest descent local search running on S, with supermove enable
##         """
##         self.fitEval = 0

##         start = os.times()[0]
##         self.model.transWal()
##         self.oldindiv = self.initIndiv(self.dim)
##         self.oldindiv = self.evalPop(self.oldindiv)
##         self.initWal(self.model)
##         #self.printInter()
        
##         self.bsf = copy.deepcopy(self.oldindiv)
        
##         self.model.WA = []
##         #self.printInter()
        
##         init = False
##         updateT = 0
##         initT = os.times()[0] - start
##         start = os.times()[0]
##         while self.fitEval < self.MaxFit:
##             if init == False:
##                 improveN, bestI = self.genFitBestsm(minimize)
##                 init = True
##             else:
##                 improveN, bestI = self.oldindiv.updateFitBestsm(minimize)

##             if improveN == False:
##                 if restart == True:
##                     updateT = updateT + os.times()[0] - start
##                     startR = os.times()[0]
##                     oldbit = self.oldindiv.bit
##                     self.oldindiv = self.evalPop(self.oldindiv)
## #                    print 'restart'

##                     self.fitEval = self.fitEval - 1
##                     self.restart(fitName, minimize, False)
##                     init = False
##                     diff = self.diffBits(oldbit, self.oldindiv.bit)
                    
##                     for i in diff:
##                         self.oldindiv.update(i)
##                         self.oldindiv.updateWAS(i)
##                     initT = initT + os.times()[0] - startR
##                     start = os.times()[0]
##                 else:
##                     self.oldindiv = self.evalPop(self.oldindiv)
##                     return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit}
##             else : # improveN is TRUE 
## #                print self.P
##                 self.fitEval = self.fitEval + len(self.P) * self.dim
##                 for i in self.P:
##                     self.oldindiv.update(i)
##                     self.oldindiv.updateWAS(i)
##                     #takenBits[bestI] = True
##     #                updateCT = updateCT + os.times()[0] - start
##     #                bestBitsCount[bestI] = bestBitsCount[bestI] + 1
##                     if self.oldindiv.bit[i] == '1':
##                         self.oldindiv.bit[i] = '0'
##                     else:
##                         self.oldindiv.bit[i] = '1'

##         self.bsf = self.evalPop(self.bsf)
##         updateT = updateT + os.times()[0] - start
##         return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit,'init':initT, 'update':updateT}

    def runFitSwalk(self,fitName, minimize, restart):
        """ 
        steepest descent local search running on S
        """
        self.fitEval = 0
        start = os.times()[0]
        self.model.transWal()
        
        self.oldindiv = individual.Individual(n=self.dim)
        self.oldindiv.init()
        self.oldindiv = self.evalPop(self.oldindiv)
        self.oldindiv.initWal(self.model)
        self.model.initInter()
        self.bsf = copy.deepcopy(self.oldindiv)
        self.oldindiv.genImproveS(minimize)
        self.model.WA = []

        initC = 1
        updateC = 0

        descT = 0
        pertT = 0
        updatePertT = 0
        updateT = 0
        self.fitEval = 0
        walkLen = 10

        initT = os.times()[0] - start

        while self.fitEval < self.MaxFit:
            start = os.times()[0]
            improveN, bestI = self.oldindiv.steepFitDesc(minimize)
            descT = descT + os.times()[0] - start
            
            if improveN == False:
                initC = initC + 1
                if restart == True:
                    start = os.times()[0]
                    diff, self.oldindiv = self.walk(fitName, minimize,False, walkLen, self.oldindiv)
                    pertT = pertT + os.times()[0] - start

                    start = os.times()[0]
                    for i in diff:
                        self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[i]
                        self.oldindiv.update(i)
                        self.oldindiv.updateWAS(i)
                        self.oldindiv.updatePertImprS(i, minimize)
                    updatePertT = updatePertT + os.times()[0] - start

                    self.fitEval = self.fitEval + len(diff)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                start = os.times()[0]
                self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[bestI] 
                self.oldindiv.update(bestI)
                self.oldindiv.updateWAS(bestI)
                self.oldindiv.updateImprS(bestI, minimize)
                self.fitEval = self.fitEval + 1
                updateT = updateT + os.times()[0] - start
                updateC = updateC + 1
                if self.oldindiv.bit[bestI] == '1':
                    self.oldindiv.bit[bestI] = '0'
                else:
                    self.oldindiv.bit[bestI] = '1'

        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit, 'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'initC':initC, 'updateC':updateC}

    def runFitSwalkNext(self,fitName, minimize, restart):
        """ 
        next descent local search running on S
        """
        self.fitEval = 0
        start = os.times()[0]
        self.model.transWal()

        self.oldindiv = individual.Individual(n=self.dim)
        self.oldindiv.init()
        self.oldindiv = self.evalPop(self.oldindiv)
        self.oldindiv.initWal(self.model)
        self.model.initInter()
        self.bsf = copy.deepcopy(self.oldindiv)
        self.oldindiv.genImproveS(minimize)
        self.model.WA = []

        initC = 1
        updateC = 0

        descT = 0
        pertT = 0
        updatePertT = 0
        updateT = 0
        self.fitEval = 0
        walkLen = 10

        initT = os.times()[0] - start

        while self.fitEval < self.MaxFit:
            start = os.times()[0]
            improveN, bestI = self.oldindiv.nextDesc()
            descT = descT + os.times()[0] - start
            
            if improveN == False:
                initC = initC + 1
                if restart == True:

                    start = os.times()[0]
                    diff, self.oldindiv = self.walk(fitName, minimize,False, walkLen, self.oldindiv)
                    pertT = pertT + os.times()[0] - start

                    start = os.times()[0]
                    for i in diff:
                        self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[i]
                        self.oldindiv.update(i)
                        self.oldindiv.updateWAS(i)
                        self.oldindiv.updatePertImprS(i, minimize)
                    updatePertT = updatePertT + os.times()[0] - start

                    self.fitEval = self.fitEval + len(diff)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                start = os.times()[0]
                self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[bestI]
                self.oldindiv.update(bestI)
                self.oldindiv.updateWAS(bestI)
                self.oldindiv.updateImprS(bestI, minimize)
                self.fitEval = self.fitEval + 1
                updateT = updateT + os.times()[0] - start
                updateC = updateC + 1
                if self.oldindiv.bit[bestI] == '1':
                    self.oldindiv.bit[bestI] = '0'
                else:
                    self.oldindiv.bit[bestI] = '1'
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit, 'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'initC':initC, 'updateC':updateC}

    def runFitS2(self,fitName, minimize, restart):
        """ 
        steepest descent local search running on S
        """
        self.fitEval = 0
        start = os.times()[0]
        self.model.transWal()
        self.oldindiv = self.initIndiv(self.dim)
        self.oldindiv = self.evalPop(self.oldindiv)
        self.initWal(self.model)
        self.bsf = copy.deepcopy(self.oldindiv)
        self.model.WA = []

        walkLen = 10
        init = False
        updateT = 0
        initT = os.times()[0] - start
        start = os.times()[0]

##     def runFitS2walk(self,fitName, minimize, restart):
##         """ 
##         steepest descent local search running on S
##         """
##         self.fitEval = 0
##         start = os.times()[0]
##         self.model.transWal()
##         self.oldindiv = self.initIndiv(self.dim)
##         self.oldindiv = self.evalPop(self.oldindiv)
##         self.initWal(self.model)
##         self.bsf = copy.deepcopy(self.oldindiv)
##         self.model.WA = []

##         walkLen = 10
##         init = False
##         updateT = 0
##         initT = os.times()[0] - start
##         start = os.times()[0]

##         while self.fitEval < self.MaxFit:
##             if init == False:
##                 improveN, bestI, evalCount = self.genFitBest2(minimize)
##                 init = True
##             else:
##                 improveN, bestI, evalCount = self.oldindiv.updateFitBest2(bestI,minimize)
##             self.fitEval = self.fitEval + evalCount

## #            print 'oldindiv',self.oldindiv.bit
## #            print 'improveA',self.improveA
## #            print improveN, bestI, self.fitEval
##             print self.oldindiv.bit
##             print improveN, bestI, self.fitEval
##             print 
## #            pdb.set_trace()
        
##             if improveN == False:
##                 if restart == True:
##                     updateT = updateT + os.times()[0] - start
##                     startR = os.times()[0]
##                     self.oldindiv = self.evalPop(self.oldindiv)
##                     diff, self.oldindiv = self.walk(fitName, minimize,False, walkLen, self.oldindiv)
##                     init = False

## #                    print self.bsf.fit
## #                    pdb.set_trace()

##                     for i in diff:
##                         self.oldindiv.update(i)
##                         self.oldindiv.updateWAS(i)
##                     initT = initT + os.times()[0] - startR
##                     start = os.times()[0]
##                 else:
##                     self.oldindiv = self.evalPop(self.oldindiv)
##                     return {'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.bsf.bit,'init':initT, 'update':updateT}
##             else : # improveN is TRUE 
##                 for i in bestI:
##                     self.oldindiv.update(i)
##                     self.oldindiv.updateWAS(i)
##                     if self.oldindiv.bit[i] == '1':
##                         self.oldindiv.bit[i] = '0'
##                     else:
##                         self.oldindiv.bit[i] = '1'

##         if fitName == 'fit' and minimize == True :
##             if self.bsf.fit > self.oldindiv.fit:
##                 self.bsf = copy.deepcopy(self.oldindiv)
##         elif fitName == 'fit' and minimize == False :
##             if self.bsf.fit < self.oldindiv.fit:
##                 self.bsf = copy.deepcopy(self.oldindiv)
##         updateT = updateT + os.times()[0] - start
##         return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit,'init':initT, 'update':updateT}

    def runFitWal(self,fitName, minimize, restart):
        """
        Walsh Analysis for speeding up steepest descent local search
        """
        self.oldindiv = self.initIndiv(self.dim)
        self.fitEval = 0
        
        self.model.transWal()
        self.indiv = copy.deepcopy(self.oldindiv)
        self.initWal(self.model)
#        self.model.initInter()

        self.oldindiv = self.evalPop(self.oldindiv)
        self.bsf = copy.deepcopy(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)

        self.model.WA = []
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
                self.indiv.fit = oldFit - 2*self.oldindiv.sumArr[nCount]
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
                        self.oldindiv.update(i)
                else:
#                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit,'trace':self.trace}
                    #print 'compPSum', compPSumT
                    #print 'update', updateT
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                self.oldindiv.update(changeBit)
                lastChangeBit = changeBit

            if init == False:
                init = True

        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit}


    def runFitSrest(self,fitName, minimize, restart):
        """ 
        steepest descent local search running on S
        """
        self.fitEval = 0
        start = os.times()[0]
        self.model.transWal()

        self.oldindiv = individual.Individual(n=self.dim)
        self.oldindiv.init()
        self.oldindiv = self.evalPop(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)
        self.oldindiv.initWal(self.model)
        self.model.initInter()
        self.bsf = copy.deepcopy(self.oldindiv)
        self.oldindiv.genImproveS(minimize)
        self.model.WA = []

        initC = 1
        updateC = 0

        descT = 0
        pertT = 0
        updatePertT = 0
        updateT = 0
        self.fitEval = 0
        walkLen = 10

        initT = os.times()[0] - start

        while self.fitEval < self.MaxFit:
            start = os.times()[0]
            improveN, bestI = self.oldindiv.steepFitDesc(minimize)
            descT = descT + os.times()[0] - start
            
            if improveN == False:
                initC = initC + 1
                if restart == True:
                    start = os.times()[0]
                    oldbit = self.oldindiv.bit
                    oldfit = self.oldindiv.fit
                    self.restart(fitName, minimize, False)
                    pertT = pertT + os.times()[0] - start

                    start = os.times()[0]
                    diff = self.diffBits(oldbit, self.oldindiv.bit)
                    self.oldindiv.fit = oldfit
                    for i in diff:
                        self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[i]
                        self.oldindiv.update(i)
                        self.oldindiv.updateWAS(i)
                        self.oldindiv.updatePertImprS(i, minimize)
                    updatePertT = updatePertT + os.times()[0] - start

                    self.fitEval = self.fitEval + len(diff)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                start = os.times()[0]

                self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[bestI] 

                self.oldindiv.update(bestI)
                self.oldindiv.updateWAS(bestI)
                self.oldindiv.updateImprS(bestI, minimize)
                self.fitEval = self.fitEval + 1
                updateT = updateT + os.times()[0] - start
                updateC = updateC + 1
                if self.oldindiv.bit[bestI] == '1':
                    self.oldindiv.bit[bestI] = '0'
                else:
                    self.oldindiv.bit[bestI] = '1'
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit, 'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'initC':initC, 'updateC':updateC}


    def runFitSrestNext(self,fitName, minimize, restart):
        """ 
        next descent local search running on S
        """
        self.fitEval = 0
        start = os.times()[0]
        self.model.transWal()
        self.oldindiv = individual.Individual(n=self.dim)
        self.oldindiv.init()
        self.oldindiv = self.evalPop(self.oldindiv)
        self.oldindiv.initWal(self.model)
        self.model.initInter()
        self.bsf = copy.deepcopy(self.oldindiv)
        self.oldindiv.genImproveS(minimize)
        self.model.WA = []

        initC = 1
        updateC = 0

        descT = 0
        pertT = 0
        updatePertT = 0
        updateT = 0
        self.fitEval = 0

        initT = os.times()[0] - start

        while self.fitEval < self.MaxFit:
            start = os.times()[0]
            improveN, bestI = self.oldindiv.nextDesc()
            descT = descT + os.times()[0] - start
            
            if improveN == False:
                initC = initC + 1
                if restart == True:
                    start = os.times()[0]
                    oldbit = self.oldindiv.bit
                    oldfit = self.oldindiv.fit
                    self.restart(fitName, minimize, False)
                    pertT = pertT + os.times()[0] - start

                    start = os.times()[0]
                    diff = self.diffBits(oldbit, self.oldindiv.bit)
                    self.oldindiv.fit = oldfit
                    for i in diff:
                        self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[i]
                        self.oldindiv.update(i)
                        self.oldindiv.updateWAS(i)
                        self.oldindiv.updatePertImprS(i, minimize)
                    updatePertT = updatePertT + os.times()[0] - start

                    self.fitEval = self.fitEval + len(diff)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                start = os.times()[0]

                self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[bestI] 

                self.oldindiv.update(bestI)
                self.oldindiv.updateWAS(bestI)
                self.oldindiv.updateImprS(bestI, minimize)
                self.fitEval = self.fitEval + 1
                updateT = updateT + os.times()[0] - start
                updateC = updateC + 1
                if self.oldindiv.bit[bestI] == '1':
                    self.oldindiv.bit[bestI] = '0'
                else:
                    self.oldindiv.bit[bestI] = '1'

        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit, 'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'initC':initC, 'updateC':updateC}

    def runMeanSCwalk(self,fitName, minimize, restart):
        """ 
        steepest descent local search with respect to mean of neighs by Walsh Analysis
        """
        self.fitEval = 0
        
        start = os.times()[0]
        self.model.transWal()
        self.oldindiv = individual.Individual(n=self.dim,neigh=True)
        self.oldindiv.init()
        self.oldindiv = self.evalPop(self.oldindiv)
        self.model.initInter()
        self.oldindiv.initWal(self.model)
        self.oldindiv.initSC()
        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (sum(self.oldindiv.sumArr))
        self.oldindiv.genImproveSC(minimize)
        self.bsf = copy.deepcopy(self.oldindiv)
        self.model.WA = []

        initC = 0
        updateC = 0

        descT = 0
        pertT = 0
        updatePertT = 0
        updateT = 0
        self.fitEval = 0
        walkLen = 10

        initT = os.times()[0] - start
        while self.fitEval < self.MaxFit:
            start = os.times()[0]
            improveN, bestI = self.oldindiv.steepMeanDesc(minimize)
            descT = descT + os.times()[0] - start

            if improveN == False:
                initC = initC + 1
                if restart == True:

                    start = os.times()[0]
                    diff, self.oldindiv = self.walk(fitName, minimize, False, walkLen, self.oldindiv)
                    pertT = pertT + os.times()[0] - start

                    start = os.times()[0]
                    for i in diff:
                        self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[i] 
                        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))

                        self.oldindiv.update(i)
                        self.oldindiv.updateSC(i)
                        self.oldindiv.updateWAS(i)
                        self.oldindiv.updatePertImprSC(i, minimize)
                    updatePertT = updatePertT + os.times()[0] - start

                    self.fitEval = self.fitEval + len(diff)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                start = os.times()[0]

                self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[bestI] 
                self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))

                self.oldindiv.update(bestI)
                self.oldindiv.updateSC(bestI)
                self.oldindiv.updateWAS(bestI)
                self.oldindiv.updateImprSC(bestI, minimize)
                self.fitEval = self.fitEval + 1
                updateT = updateT + os.times()[0] - start
                updateC = updateC + 1
                if self.oldindiv.bit[bestI] == '1':
                    self.oldindiv.bit[bestI] = '0'
                else:
                    self.oldindiv.bit[bestI] = '1'

        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'fitG': self.bsf.fitG, 'bit':self.bsf.bit,'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'initC':initC, 'updateC':updateC}

    def runMeanSCwalkNext(self,fitName, minimize, restart):
        """ 
        next descent local search with respect to mean of neighs by Walsh Analysis
        """
        self.fitEval = 0
        
        start = os.times()[0]
        self.model.transWal()
        self.oldindiv = individual.Individual(n=self.dim,neigh=True)
        self.oldindiv.init()
        self.oldindiv = self.evalPop(self.oldindiv)
        self.model.initInter()
        self.oldindiv.initWal(self.model)
        self.oldindiv.initSC()
        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (sum(self.oldindiv.sumArr))
        self.oldindiv.genImproveSC(minimize)
        self.bsf = copy.deepcopy(self.oldindiv)
        self.model.WA = []

        initC = 0
        updateC = 0

        descT = 0
        pertT = 0
        updatePertT = 0
        updateT = 0
        self.fitEval = 0
        walkLen = 10

        initT = os.times()[0] - start
        while self.fitEval < self.MaxFit:
            start = os.times()[0]
            improveN, bestI = self.oldindiv.nextDesc()
            descT = descT + os.times()[0] - start

            if improveN == False:
                initC = initC + 1
                if restart == True:

                    start = os.times()[0]
                    diff, self.oldindiv = self.walk(fitName, minimize,False, walkLen, self.oldindiv)
                    pertT = pertT + os.times()[0] - start
                    
                    start = os.times()[0]
                    for i in diff:
                        self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[i] 
                        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))

                        self.oldindiv.update(i)
                        self.oldindiv.updateSC(i)
                        self.oldindiv.updateWAS(i)
                        self.oldindiv.updatePertImprSC(i, minimize)
                    updatePertT = updatePertT + os.times()[0] - start

                    self.fitEval = self.fitEval + len(diff)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                start = os.times()[0]
                self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[bestI] 
                self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))
                self.oldindiv.update(bestI)
                self.oldindiv.updateSC(bestI)
                self.oldindiv.updateWAS(bestI)
                self.oldindiv.updateImprSC(bestI, minimize)
                self.fitEval = self.fitEval + 1
                updateT = updateT + os.times()[0] - start
                updateC = updateC + 1
                if self.oldindiv.bit[bestI] == '1':
                    self.oldindiv.bit[bestI] = '0'
                else:
                    self.oldindiv.bit[bestI] = '1'

        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'fitG': self.bsf.fitG, 'bit':self.bsf.bit,'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'initC':initC, 'updateC':updateC}

    def runMeanSCrest(self,fitName, minimize, restart):
        """ 
        steepest descent local search with respect to mean of neighs by Walsh Analysis
        """
        self.fitEval = 0
        
        start = os.times()[0]
        self.model.transWal()
        self.oldindiv = individual.Individual(n=self.dim,neigh=True)
        self.oldindiv.init()
        self.oldindiv = self.evalPop(self.oldindiv)
        self.model.initInter()
        self.oldindiv.initWal(self.model)
        self.oldindiv.initSC()
        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (sum(self.oldindiv.sumArr))
        self.oldindiv.genImproveSC(minimize)
        self.bsf = copy.deepcopy(self.oldindiv)
        self.model.WA = []

        initC = 0
        updateC = 0

        descT = 0
        pertT = 0
        updatePertT = 0
        updateT = 0
        self.fitEval = 0

        initT = os.times()[0] - start
        while self.fitEval < self.MaxFit:
            start = os.times()[0]
            improveN, bestI = self.oldindiv.steepMeanDesc(minimize)
            descT = descT + os.times()[0] - start

            if improveN == False:
                initC = initC + 1
                if restart == True:
                    start = os.times()[0]
                    oldbit = self.oldindiv.bit
                    oldfit = self.oldindiv.fit
                    self.restart(fitName, minimize, False)
                    pertT = pertT + os.times()[0] - start

                    start = os.times()[0]
                    diff = self.diffBits(oldbit, self.oldindiv.bit)
                    self.oldindiv.fit = oldfit  
                    for i in diff:
                        self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[i] 
                        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))

                        self.oldindiv.update(i)
                        self.oldindiv.updateSC(i)
                        self.oldindiv.updateWAS(i)
                        self.oldindiv.updatePertImprS(i, minimize)
                    updatePertT = updatePertT + os.times()[0] - start

                    self.fitEval = self.fitEval + len(diff)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                start = os.times()[0]

                self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[bestI] 
                self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))

                self.oldindiv.update(bestI)
                self.oldindiv.updateSC(bestI)
                self.oldindiv.updateWAS(bestI)
                self.oldindiv.updateImprSC(bestI, minimize)
                self.fitEval = self.fitEval + 1
                updateT = updateT + os.times()[0] - start
                updateC = updateC + 1
                if self.oldindiv.bit[bestI] == '1':
                    self.oldindiv.bit[bestI] = '0'
                else:
                    self.oldindiv.bit[bestI] = '1'

        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'fitG': self.bsf.fitG, 'bit':self.bsf.bit,'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'initC':initC, 'updateC':updateC}

    def runMeanSCrestNext(self,fitName, minimize, restart):
        """ 
        steepest descent local search with respect to mean of neighs by Walsh Analysis
        """
        self.fitEval = 0
        
        start = os.times()[0]
        self.model.transWal()
        self.oldindiv = individual.Individual(n=self.dim,neigh=True)
        self.oldindiv.init()
        self.oldindiv = self.evalPop(self.oldindiv)
        self.model.initInter()
        self.oldindiv.initWal(self.model)
        self.oldindiv.initSC()
        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (sum(self.oldindiv.sumArr))
        self.oldindiv.genImproveSC(minimize)
        self.bsf = copy.deepcopy(self.oldindiv)
        self.model.WA = []

        initC = 0
        updateC = 0

        descT = 0
        pertT = 0
        updatePertT = 0
        updateT = 0
        self.fitEval = 0

        initT = os.times()[0] - start
        while self.fitEval < self.MaxFit:
            start = os.times()[0]
            improveN, bestI = self.oldindiv.nextDesc()
            descT = descT + os.times()[0] - start

            if improveN == False:
                initC = initC + 1
                if restart == True:
                    start = os.times()[0]
                    oldbit = self.oldindiv.bit
                    oldfit = self.oldindiv.fit
                    self.restart(fitName, minimize, False)
                    pertT = pertT + os.times()[0] - start

                    start = os.times()[0]
                    diff = self.diffBits(oldbit, self.oldindiv.bit)
                    self.oldindiv.fit = oldfit  
                    for i in diff:
                        self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[i] 
                        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))

                        self.oldindiv.update(i)
                        self.oldindiv.updateSC(i)
                        self.oldindiv.updateWAS(i)
                        self.oldindiv.updatePertImprS(i, minimize)
                    updatePertT = updatePertT + os.times()[0] - start

                    self.fitEval = self.fitEval + len(diff)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                start = os.times()[0]

                self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[bestI] 
                self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))

                self.oldindiv.update(bestI)
                self.oldindiv.updateSC(bestI)
                self.oldindiv.updateWAS(bestI)
                self.oldindiv.updateImprSC(bestI, minimize)
                self.fitEval = self.fitEval + 1
                updateT = updateT + os.times()[0] - start
                updateC = updateC + 1
                if self.oldindiv.bit[bestI] == '1':
                    self.oldindiv.bit[bestI] = '0'
                else:
                    self.oldindiv.bit[bestI] = '1'

        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'fitG': self.bsf.fitG, 'bit':self.bsf.bit,'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'initC':initC, 'updateC':updateC}
    
    def hyperSearchMean(self,fitName, minimize, restart):
        """ 
        performing hyper search using the probability generated by Hyperplane analysis
        """

        self.fitEval = 0
        
        start = os.times()[0]
        self.oldindiv = Struct( fit = 0, fitG = 0, bit = self.genSolProp(self.model.sumFitA) )
        self.oldindiv = self.evalPop(self.oldindiv)
        self.initWal(self.model)
        
        self.initSC()
        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (sum(self.oldindiv.sumArr))

        self.bsf = copy.deepcopy(self.oldindiv)
        #self.model.WA = []

        init = False
        updateT = 0
        walkLen = 10
        initT = os.times()[0] - start
        start = os.times()[0]
        while self.fitEval < self.MaxFit:
            if init == False:
                improveN, bestI, evalCount = self.genMeanBest(minimize)
                init = True
            else :
                improveN, bestI, evalCount = self.oldindiv.updateMeanBest(bestI,minimize)
            self.fitEval = self.fitEval + self.dim

            if improveN == False:
                if restart == True:
                    updateT = updateT + os.times()[0] - start
                    startR = os.times()[0]
                    self.oldindiv = self.evalPop(self.oldindiv)
                    self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))
                    
                    diff, self.oldindiv = self.walk( fitName, minimize, False, walkLen , self.oldindiv, self.oldindiv)
                    init = False
                    
                    for i in diff:
                        self.oldindiv.update(i)

                        self.oldindiv.updateSC(i)
                        self.oldindiv.updateWAS(i)
                    initT = initT + os.times()[0] - startR
                    start = os.times()[0]
                else:
                    self.oldindiv = self.evalPop(self.oldindiv)
                    self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))
                    self.fitEval = self.fitEval - 1
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                self.oldindiv.update(bestI)
                self.oldindiv.updateSC(bestI)

                self.oldindiv.updateWAS(bestI)
                if self.oldindiv.bit[bestI] == '1':
                    self.oldindiv.bit[bestI] = '0'
                else:
                    self.oldindiv.bit[bestI] = '1'

        self.bsf = self.evalPop(self.bsf)
        self.fitEval = self.fitEval - 1
        diff = self.diffBits(self.bsf.bit, self.oldindiv.bit)
        for i in diff:
            self.oldindiv.update(i)
        self.bsf.fitG = self.bsf.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))
        updateT = updateT + os.times()[0] - start
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'fitG': self.bsf.fitG, 'bit':self.bsf.bit,'init':initT, 'update':updateT}

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


    def runWalSearch(self,fitName, minimize, restart):
        """
        initial assignment according Walsh coefficients
        """
        self.combLoopup = dict()
        self.model.transWal()
        self.WAsort = sorted(self.model.WA, key=lambda i: abs(i.w), reverse=True)
        sol = []
        determineBit = []
        for i in range(self.dim) :
            sol.append('0')
            determineBit.append(False)
#        self.printWAsort()
        poss = []
        for i in range(1,len(self.WAsort)):
            arr = self.WAsort[i].arr
            print 
            print 'arr',arr
            
            if (minimize == True and self.WAsort[i].w < -self.threshold) or (minimize == False and self.WAsort[i].w > self.threshold):
                odd = False
            else :
                odd = True
            iters = self.genPossBit(odd, arr)

            
            # reduce and reconstruct the possible cases
            if len(poss)!=0:
                copyPoss = copy.copy(poss)
                for j in range(len(copyPoss)):
                    print 'j.a',copyPoss[j].a
                    join = list(Set(arr).intersection(copyPoss[j].a)) #[filter(lambda x: x in arr, sublist) for sublist in j.a]
                    print 'join', join
                    # for every bit in jion set, both arrs should be identical
                    tempArr = []
                    if len(join)!=0:
                        conflict = True 
                        for ii in iters:
                            for jj in copyPoss[j].i:
                                itent = True
                                for k in join:
                                    # check whether all bits in intersection bits are the same
                                    print 'ii', ii, '\tjj', jj, '\tk', k
                                    if bool(k in ii) ^ bool(k in jj): 
                                        # bits in the overlapping position are not identical
                                        itent = False
                                        break
                                    print 'identical', itent
                                if itent == True:
                                    # reconstruct the possible bitstring
                                    tempArr.append(list(Set(ii).union(Set(jj))))
                                    conflict = False
                                    print 'tempArr', tempArr
                        if conflict == False:
                            poss.pop(j) # TODO may be problematic
                            print 'join arr', list(Set(copyPoss[j].a).union(Set(arr)))
                            poss.append( Struct(i=copy.deepcopy(tempArr), a=list(Set(copyPoss[j].a).union(Set(arr))) ))
                    else:
                        poss.append(Struct(i=iters, a=arr))

                    if len(poss[-1].i) == 1:
                        for k in poss[-1].i[0]:
                            sol[k] = '1'
                        for k in poss[-1].a:
                            determineBit[k] = True  

                    print 'determineBit', determineBit
                    if False not in determineBit:
                        sol =  ''.join(sol)
                        print sol
                        print '%.3e' %(self.func(sol))
                        sys.exit()
            else:
                poss.append(Struct(i=iters, a=arr))

            if len(arr)==1:
                sol[arr[0]] = str(int(odd))
                determineBit[arr[0]] = True
            print 'len',len(poss)

            if len(arr) == 1:
                if odd == False :
                    sol[arr[0]] = '0'
                else:
                    sol[arr[0]] = '1'

            pdb.set_trace()
             
        sol =  ''.join(sol)
        print sol
        print  '%.3e' %(self.func(sol))
        print 
        return { 'nEvals': 1, 'sol': None, 'bit': None}

    cdef runBeamFitSwalkNext(self,fitName, minimize, restart, beamWidth):
        """ 
        BEAM next descent local search running on S, 
        where beamWidth is the number of beams
        """
        self.fitEval = 0
        start = os.times()[0]
        self.model.transWal()
        self.model.initInter()
        self.oldpop = []
        for i in xrange(beamWidth):
            self.oldpop.append(individual.Individual(n=self.dim))
            self.oldpop[i].init()
            self.oldpop[i] = self.evalPop(self.oldpop[i])
            self.oldpop[i].initWal(self.model)
            self.oldpop[i].genImproveS(minimize)

        # self.bsf is a single individual recording the best-so-far solution
        #self.bsf = copy.deepcopy(self.oldpop[self.bestIndiv(self.oldpop, minimize, fitName)])
        self.bsf = self.oldpop[self.bestIndiv(self.oldpop, minimize, fitName)]
        self.model.WA = [] # release memory

        initC = 1
        updateC = 0

        descT = 0
        pertT = 0
        updatePertT = 0
        updateT = 0
        self.fitEval = 0
        walkLen = 10

        initT = os.times()[0] - start

        while self.fitEval < self.MaxFit:
            print self.fitEval
            for i in xrange(len(self.oldpop)):
                start = os.times()[0]
                improveN, bestI = self.oldpop[i].nextDesc()
                descT = descT + os.times()[0] - start

                if improveN == False:
                    initC = initC + 1
                    if restart == True:
                        start = os.times()[0]
                        diff, self.oldpop[i] = self.walk(fitName, minimize, False, walkLen, self.oldpop[i])

                        pertT = pertT + os.times()[0] - start

                        start = os.times()[0]
                        for j in diff:
                            #                            self.oldpop[i].fit = self.oldpop[i].fit - 2*self.oldpop[i].sumArr[j]
                            print 'begin sumArr'
                            self.oldpop[i].updateSumArr(j)
                            print 'begin update'
                            self.oldpop[i].update(j)
                            print 'begin updateWAS'
                            self.oldpop[i].updateWAS(j)
                            print 'begin updatePertImprS'
                            self.oldpop[i].updatePertImprS(j, minimize)
                            print 'end'
                        updatePertT = updatePertT + os.times()[0] - start

                        self.fitEval = self.fitEval + len(diff)
                    else:
                        return { 'nEvals': self.fitEval, 'sol': self.oldpop[i].fit, 'bit':self.oldpop[i].bit}
                else : # improveN is TRUE 
                    start = os.times()[0]
                    #                    self.oldpop[i].fit = self.oldpop[i].fit - 2*self.oldpop[i].sumArr[bestI]
                    print 'begin sumArr'
                    self.oldpop[i].updateSumArr(bestI)
                    print 'begin update'
                    self.oldpop[i].update(bestI)
                    print 'begin updateWAS'
                    self.oldpop[i].updateWAS(bestI)
                    print 'begin sumArrImprS'
                    self.oldpop[i].updateImprS(bestI, minimize)
                    print 'end'
                    self.fitEval = self.fitEval + 1
                    updateT = updateT + os.times()[0] - start
                    updateC = updateC + 1
                    self.oldpop[i].flip(bestI)
                    ## if self.oldpop[i].bit[bestI] == '1':
                    ##     self.oldpop[i].bit[bestI] = '0'
                    ## else:
                    ##     self.oldpop[i].bit[bestI] = '1'
        print 'return'
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit, 'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'initC':initC, 'updateC':updateC}

    def runBeamFitSwalk(self,fitName, minimize, restart, beamWidth):
        """ 
        BEAM steepest descent local search running on S, 
        where beamWidth is the number of beams
        """

        self.fitEval = 0
        start = os.times()[0]
        self.model.transWal()
        self.model.initInter()
        self.oldpop = []
        for i in xrange(beamWidth):
            self.oldpop.append(individual.Individual(n=self.dim)) 
            self.oldpop[i].init()
            self.oldpop[i] = self.evalPop(self.oldpop[i])
            self.oldpop[i].initWal(self.model)
            self.oldpop[i].genImproveS(minimize)

        # self.bsf is a single individual recording the best-so-far solution
        self.bsf = copy.deepcopy(self.oldpop[self.bestIndiv(self.oldpop, minimize, fitName)])
        self.model.WA = [] # release memory

        initC = 1
        updateC = 0

        descT = 0
        pertT = 0
        updatePertT = 0
        updateT = 0
        self.fitEval = 0
        walkLen = 10

        initT = os.times()[0] - start

        while self.fitEval < self.MaxFit:
            for i in xrange(len(self.oldpop)):
                start = os.times()[0]
                improveN, bestI = self.oldpop[i].steepFitDesc(minimize)
                descT = descT + os.times()[0] - start
                
                if improveN == False:
                    initC = initC + 1
                    if restart == True:
                        start = os.times()[0]
                        diff, self.oldpop[i] = self.walk(fitName, minimize, False, walkLen, self.oldpop[i])


                        pertT = pertT + os.times()[0] - start

                        start = os.times()[0]
                        for j in diff:
                            self.oldpop[i].fit = self.oldpop[i].fit - 2*self.oldpop[i].sumArr[j]
                            self.oldpop[i].update(j)
                            self.oldpop[i].updateWAS(j)
                            self.oldpop[i].updatePertImprS(j, minimize)
                        updatePertT = updatePertT + os.times()[0] - start

                        self.fitEval = self.fitEval + len(diff)
                    else:
                        return { 'nEvals': self.fitEval, 'sol': self.oldpop[i].fit, 'bit':self.oldpop[i].bit}
                else : # improveN is TRUE 
                    start = os.times()[0]
                    self.oldpop[i].fit = self.oldpop[i].fit - 2*self.oldpop[i].sumArr[bestI]
                    self.oldpop[i].update(bestI)
                    self.oldpop[i].updateWAS(bestI)
                    self.oldpop[i].updateImprS(bestI, minimize)
                    self.fitEval = self.fitEval + 1
                    updateT = updateT + os.times()[0] - start
                    updateC = updateC + 1
                    if self.oldpop[i].bit[bestI] == '1':
                        self.oldpop[i].bit[bestI] = '0'
                    else:
                        self.oldpop[i].bit[bestI] = '1'
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit, 'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'initC':initC, 'updateC':updateC}

    def runBeamMeanSCwalkNext(self,fitName, minimize, restart, beamWidth):
        """ 
        BEAM next descent local search running on SC, 
        where beamWidth is the number of beams
        """

        self.fitEval = 0
        start = os.times()[0]
        self.model.transWal()
        self.model.initInter()
        self.oldpop = []
        for i in xrange(beamWidth):
            self.oldpop.append(individual.Individual(n=self.dim, neigh=True)) 
            self.oldpop[i].init()
            self.oldpop[i] = self.evalPop(self.oldpop[i])
            self.oldpop[i].initWal(self.model)
            self.oldpop[i].initSC()
            self.oldpop[i].fitG = self.oldpop[i].fit - 2/float(self.dim) * (sum(self.oldpop[i].sumArr))
            self.oldpop[i].genImproveSC(minimize)

        # self.bsf is a single individual recording the best-so-far solution
        self.bsf = copy.deepcopy(self.oldpop[self.bestIndiv(self.oldpop, minimize, fitName)])
        self.model.WA = [] # release memory

        initC = 1
        updateC = 0

        descT = 0
        pertT = 0
        updatePertT = 0
        updateT = 0
        self.fitEval = 0
        walkLen = 10

        initT = os.times()[0] - start

        while self.fitEval < self.MaxFit:
            for i in xrange(len(self.oldpop)):
                start = os.times()[0]
                improveN, bestI = self.oldpop[i].nextDesc()
                descT = descT + os.times()[0] - start
                
                if improveN == False:
                    initC = initC + 1
                    if restart == True:
                        start = os.times()[0]
                        diff, self.oldpop[i] = self.walk(fitName, minimize, False, walkLen, self.oldpop[i])

                        pertT = pertT + os.times()[0] - start

                        start = os.times()[0]
                        for j in diff:
                            self.oldpop[i].fit = self.oldpop[i].fit - 2*self.oldpop[i].sumArr[j]
                            self.oldpop[i].fitG = self.oldpop[i].fit - 2/float(self.dim) * (np.sum(self.oldpop[i].sumArr))
                            self.oldpop[i].update(j)
                            self.oldpop[i].updateSC(j)
                            self.oldpop[i].updateWAS(j)
                            self.oldpop[i].updatePertImprSC(j, minimize)
                        updatePertT = updatePertT + os.times()[0] - start

                        self.fitEval = self.fitEval + len(diff)
                    else:
                        return { 'nEvals': self.fitEval, 'sol': self.oldpop[i].fit, 'fitG': self.oldpop[i].fitG, 'bit':self.oldpop[i].bit}
                else : # improveN is TRUE 
                    start = os.times()[0]
                    self.oldpop[i].fit = self.oldpop[i].fit - 2*self.oldpop[i].sumArr[bestI]
                    self.oldpop[i].fitG = self.oldpop[i].fit - 2/float(self.dim) * (np.sum(self.oldpop[i].sumArr))
                    self.oldpop[i].update(bestI)
                    self.oldpop[i].updateSC(bestI)
                    self.oldpop[i].updateWAS(bestI)
                    self.oldpop[i].updateImprSC(bestI, minimize)
                    self.fitEval = self.fitEval + 1
                    updateT = updateT + os.times()[0] - start
                    updateC = updateC + 1
                    if self.oldpop[i].bit[bestI] == '1':
                        self.oldpop[i].bit[bestI] = '0'
                    else:
                        self.oldpop[i].bit[bestI] = '1'
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit,'fitG': self.bsf.fitG, 'bit':self.bsf.bit, 'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'initC':initC, 'updateC':updateC}

    def runBeamMeanSCwalk(self,fitName, minimize, restart, beamWidth):
        """ 
        BEAM steepest descent local search running on SC, 
        where beamWidth is the number of beams
        """

        self.fitEval = 0
        start = os.times()[0]
        self.model.transWal()
        self.model.initInter()
        self.oldpop = []
        for i in xrange(beamWidth):
            self.oldpop.append(individual.Individual(n=self.dim, neigh=True)) 
            self.oldpop[i].init()
            self.oldpop[i] = self.evalPop(self.oldpop[i])
            self.oldpop[i].initWal(self.model)
            self.oldpop[i].initSC()
            self.oldpop[i].fitG = self.oldpop[i].fit - 2/float(self.dim) * (sum(self.oldpop[i].sumArr))
            self.oldpop[i].genImproveSC(minimize)

        # self.bsf is a single individual recording the best-so-far solution
        self.bsf = copy.deepcopy(self.oldpop[self.bestIndiv(self.oldpop, minimize, fitName)])
        self.model.WA = [] # release memory

        initC = 1
        updateC = 0

        descT = 0
        pertT = 0
        updatePertT = 0
        updateT = 0
        self.fitEval = 0
        walkLen = 10

        initT = os.times()[0] - start

        while self.fitEval < self.MaxFit:
            for i in xrange(len(self.oldpop)):
                start = os.times()[0]
                improveN, bestI = self.oldpop[i].steepMeanDesc(minimize)
                descT = descT + os.times()[0] - start
                
                if improveN == False:
                    initC = initC + 1
                    if restart == True:
                        start = os.times()[0]
                        diff, self.oldpop[i] = self.walk(fitName, minimize, False, walkLen, self.oldpop[i])

                        pertT = pertT + os.times()[0] - start

                        start = os.times()[0]
                        for j in diff:
                            self.oldpop[i].fit = self.oldpop[i].fit - 2*self.oldpop[i].sumArr[j]
                            self.oldpop[i].fitG = self.oldpop[i].fit - 2/float(self.dim) * (np.sum(self.oldpop[i].sumArr))
                            self.oldpop[i].update(j)
                            self.oldpop[i].updateSC(j)
                            self.oldpop[i].updateWAS(j)
                            self.oldpop[i].updatePertImprSC(j, minimize)
                        updatePertT = updatePertT + os.times()[0] - start

                        self.fitEval = self.fitEval + len(diff)
                    else:
                        return { 'nEvals': self.fitEval, 'sol': self.oldpop[i].fit, 'fitG': self.oldpop[i].fitG, 'bit':self.oldpop[i].bit}
                else : # improveN is TRUE 
                    start = os.times()[0]
                    self.oldpop[i].fit = self.oldpop[i].fit - 2*self.oldpop[i].sumArr[bestI]
                    self.oldpop[i].fitG = self.oldpop[i].fit - 2/float(self.dim) * (np.sum(self.oldpop[i].sumArr))
                    self.oldpop[i].update(bestI)
                    self.oldpop[i].updateSC(bestI)
                    self.oldpop[i].updateWAS(bestI)
                    self.oldpop[i].updateImprSC(bestI, minimize)
                    self.fitEval = self.fitEval + 1
                    updateT = updateT + os.times()[0] - start
                    updateC = updateC + 1
                    if self.oldpop[i].bit[bestI] == '1':
                        self.oldpop[i].bit[bestI] = '0'
                    else:
                        self.oldpop[i].bit[bestI] = '1'
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit,'fitG': self.bsf.fitG, 'bit':self.bsf.bit, 'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'initC':initC, 'updateC':updateC}

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

        if fitName == 'fit':
            self.oldindiv.init()
            if evaluate == True:
                self.oldindiv = self.evalPop(self.oldindiv)
        else :
            self.oldindiv.init()
            if evaluate == True:
                self.oldindiv = self.evalPopNeigh(self.oldindiv, fitName, minimize)

    def hyperRestart(self, fitName, minimize, evaluate):
        """
        instead of random restart, generate the restart point according to probability
        """
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
            self.oldindiv = Struct( fit = 0, bit = self.genSolProp(self.model.sumFitA) )
            if evaluate == True:
                self.oldindiv = self.evalPop(self.oldindiv)
        else :
            self.oldindiv = Struct( fit = 0, fitG = 0, bit = self.genSolProp(self.model.sumFitA) )
            if evaluate == True:
                self.oldindiv = self.evalPopNeigh(self.oldindiv, fitName, minimize)


    def walk(self, fitName, minimize, evaluate, length, oldindiv):
        # update the bsf solution
        if fitName == 'fit' and minimize == True :
            if self.bsf.fit > oldindiv.fit:
                self.bsf = copy.deepcopy(oldindiv)
        elif fitName == 'fit' and minimize == False :
            if self.bsf.fit < oldindiv.fit:
                self.bsf = copy.deepcopy(oldindiv)
        elif minimize == True :
            if self.bsf.fitG > oldindiv.fitG:
                self.bsf = copy.deepcopy(oldindiv)
        elif minimize == False :
            if self.bsf.fitG < oldindiv.fitG:
                self.bsf = copy.deepcopy(oldindiv)

        flipBits = random.sample(xrange(self.dim), length)
        for i in flipBits:
            oldindiv.flip(i)
            ## if oldindiv.bit[i] == '1':
            ##     oldindiv.bit[i] = '0'
            ## else:
            ##     oldindiv.bit[i] = '1'
        return flipBits, oldindiv

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
        #return copy.deepcopy(indiv)
        return indiv

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

    def binCount(self, arr, bit):
        """
        count the one bit of union self.model.WA[i].arr and bit
        """
        s = 0
        for i in arr:
            if bit[i] == '1':
                s = s + 1
        return s

    def binCountArr(self, a1, a2):
        """
        count the number of one bits appearing in both a1 and a2
        """
        s = 0
        for i in a1:
            if i in a2:
                s = s + 1

        return s

#    def initWal(self):
#        """ 
#        1. 
#        compute the sum array for the first time, according to the initial solution
#        2. 
#        initialize a full matrix, Coincidence, C_ij, which is a symmetric one, and diagonal is empty (as we don't want to recompute SumArr), empty whenever i >= j.
#        Two options for implementing C matrix
#            a. precompute C matrix 
#            b. construct it on the fly using dict()
#        3. 
#        initialize interaction list for each variable
#        4.
#        initialize a dict of interaction structure, where interactive bits and the index of WAS (walsh including sign)
#        """
#        self.oldindiv.sumArr = np.zeros(self.dim)
#        self.WAS = np.tile(Struct(arr = [], w = 0), len(self.model.w.keys()))# Walsh coefficients with sign, represented in Array
#        self.lookup = dict()
#        self.infectBit = dict()
#        self.C = np.zeros((self.dim,self.dim)) # coincidence matrix
#        self.Inter = dict()
#
##        self.InterCount = np.zeros(self.dim)
#
##        for i in range(self.dim):
##            self.Inter.append(Set()) # coincidence matrix
#
#        #self.C = dict() # coincidence matrix
#        for i in range(len(self.model.WA)):
#            W = int(math.pow(-1, self.binCount(self.model.WA[i].arr, self.indiv.bit))) * self.model.WA[i].w
#            self.WAS[i] = Struct(arr = self.model.WA[i].arr, w = W)
#            comb = self.genComb(len(self.model.WA[i].arr))
#            #print i, self.model.WA[i].arr, comb
#
#            for j in self.model.WA[i].arr:
#                self.oldindiv.sumArr[j] = self.oldindiv.sumArr[j] + W
#                if len(self.model.WA[i].arr)>1: # for at least order Walsh terms
#                    if j not in self.Inter: # the entry of i doesn't exist yet
#                        self.Inter[j] = Struct(arr=Set(), WI=Set())
#
#                    for k in self.model.WA[i].arr:
#                        if k != j:
#                            self.Inter[j].arr.add(k)
#                    self.Inter[j].WI.add(i)
#
#                # add list of order >= 3 Walsh terms for the purpose of updating C matrix
#                if len(self.model.WA[i].arr) >= 3:
#                    if j not in self.infectBit: 
#                        self.infectBit[j] = [Struct(arr=self.model.WA[i].arr, WI=i)]
#                    else :
#                        self.infectBit[j].append(Struct(arr=self.model.WA[i].arr, WI=i))
#
#            for j in comb: # for each list in comb
#                j0 = self.model.WA[i].arr[int(j[0])]
#                j1 = self.model.WA[i].arr[int(j[1])]
#                self.C[j0,j1] = self.C[j0,j1] + W
        

#    def initWal(self, indiv):
#        """ 
#        1. 
#        compute the sum array for the first time, according to the initial solution
#        2. 
#        initialize a full matrix, Coincidence, C_ij, which is a symmetric one, and diagonal is empty (as we don't want to recompute SumArr), empty whenever i >= j.
#        Two options for implementing C matrix
#            a. precompute C matrix 
#            b. construct it on the fly using dict()
#        3. 
#        initialize interaction list for each variable
#        4.
#        initialize a dict of interaction structure, where interactive bits and the index of WAS (walsh including sign)
#        """
#        self.oldindiv.sumArr = np.zeros(self.dim)
#        self.WAS = np.tile(Struct(arr = [], w = 0), len(self.model.w.keys()))# Walsh coefficients with sign, represented in Array
#        self.lookup = dict()
#        self.infectBit = dict()
#        self.C = np.zeros((self.dim,self.dim)) # coincidence matrix
#        self.Inter = dict()
#
#        for i in range(len(self.model.WA)):
#            W = int(math.pow(-1, self.binCount(self.model.WA[i].arr, self.indiv.bit))) * self.model.WA[i].w
#            self.WAS[i] = Struct(arr = self.model.WA[i].arr, w = W)
#            comb = self.genComb(len(self.model.WA[i].arr))
#
#            for j in self.model.WA[i].arr:
#                self.oldindiv.sumArr[j] = self.oldindiv.sumArr[j] + W
#                if len(self.model.WA[i].arr)>1: # for at least order Walsh terms
#                    if j not in self.Inter: # the entry of i doesn't exist yet
#                        self.Inter[j] = Struct(arr=Set(), WI=Set())
#
#                    for k in self.model.WA[i].arr:
#                        if k != j:
#                            self.Inter[j].arr.add(k)
#                    self.Inter[j].WI.add(i)
#
#                # add list of order >= 3 Walsh terms for the purpose of updating C matrix
#                if len(self.model.WA[i].arr) >= 3:
#                    if j not in self.infectBit: 
#                        self.infectBit[j] = [Struct(arr=self.model.WA[i].arr, WI=i)]
#                    else :
#                        self.infectBit[j].append(Struct(arr=self.model.WA[i].arr, WI=i))
#
#            for j in comb: # for each list in comb
#                j0 = self.model.WA[i].arr[int(j[0])]
#                j1 = self.model.WA[i].arr[int(j[1])]
#                self.C[j0,j1] = self.C[j0,j1] + W
        


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
            print k, '%.3f' % (self.model.w[k])

    def printWsort(self):
        """
        print all the sorted Walsh terms
        """
#        sorted(self.model.w.values())
        a = sorted(self.model.w.iteritems(), key=itemgetter(1), reverse=True)
        for i in a:
            print i[0], '%.3f' %(i[1])

    def printWAsort(self):
        for i in self.WAsort:
            print i.arr,'\t\t%.3f' %(i.w)

    def printWA(self):
        """
        print all walsh terms with array representation
        """
        for i in range(len(self.model.WA)):
           print self.model.WA[i].arr, self.model.WA[i].w

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
        self.oldindiv.sumArrFake = copy.deepcopy(self.oldindiv.sumArr)
        self.Cfake = copy.deepcopy(self.C)
        self.WASfake = copy.deepcopy(self.WAS)

    def genCombOne(self, odd, order):
        """ generate the number of possible ones, 
            given whether it should be odd or not,
            and the order of Walsh coefficients
        """
        if odd == True:
            if (odd, order) not in self.combLoopup:
                self.combLoopup[odd, order] = 2*np.array(range((order+1)/2)) + 1
        else:
            if (odd, order) not in self.combLoopup:
                self.combLoopup[odd, order] = 2*np.array(range((order)/2+1))

        return copy.deepcopy(self.combLoopup[odd, order])

    def genPossBit(self, odd, arr):
        comb = self.genCombOne(odd, len(arr))
        iters = []
        for i in comb:
            #print 'comb', i
            for j in it.combinations(arr, i):
           #     print 'j', list(j)
                iters.append(list(j))
           # print 'temp', temp
#            iters.append(copy.deepcopy(temp))
#            a = self.combinations(arr, i)
#            pdb.set_trace()
            #iters.append(self.combinations(arr, i))
        return iters
    
    def bestIndiv(self, pop, minimize, fitName):
        """ 
        return the index of best individual in the given pop,
        with respect to minimize and fitName
        """
        for i in xrange(len(pop)): 
            if fitName == 'fit':
                val = pop[i].fit
            elif fitName == 'mean':
                val = pop[i].fitG
            else :
                print 'Unknown fitName', fitName
                sys.exit(-1)

            if i==0: # first element, initialization
                best = val
                I = 0

            if (minimize == True and val < best) or (minimize == False and val > best): 
                best = val
                I = i

        return I
    
    def printPop(self, pop, fitName):
        for i in xrange(len(pop)):
            if fitName == 'fit':
                print i, pop[i].fit
            elif fitName == 'mean':
                print i, pop[i].fitG

    def printIndiv(self, indiv, fitName):
        if fitName == 'fit':
            print indiv.fit
        elif fitName == 'mean':
            print indiv.fitG

