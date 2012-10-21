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
cimport cython
from cython.parallel import prange, parallel, threadid
from libcpp.vector cimport vector
from libcpp.set cimport set
from libc.stdlib cimport malloc, free
from cython.operator cimport dereference as deref, preincrement as inc

cdef extern from "string.h":
    char *strcpy(char *dest, char *src)

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
    cdef object oldindiv
    cdef object tempIndiv
    cdef object indiv

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
        elif compM == 'bfUpdate':
            if fitName == 'fit':
                return self.runFitUpdate(fitName, minimize, restart)
        elif compM == 'walWalk':
            if fitName == 'fit':
                return self.runFitSwalk(fitName, minimize, restart)
            elif fitName == 'mean':
                return self.runMeanSCwalk(fitName, minimize, restart)
            elif fitName == 'switchF':
                return self.runSwitchFwalk('mean', minimize, restart)
            elif fitName == 'switchAvg':
                return self.runSwitchAvgwalk('mean', minimize, restart)
            elif fitName == 'combF':
                return self.runCombFwalk('mean', minimize, restart)
            elif fitName == 'combAvg':
                return self.runCombAvgwalk('mean', minimize, restart)
            elif fitName == 'fitCheck':
                return self.fitCheck('fit', minimize, restart)
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
        elif compM == 'walRestFlip':
            if fitName == 'fit':
                return self.runFitSrestFlip(fitName, minimize, restart)
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


    def runFitSwalk(self,fitName, minimize, restart):
        """
        steepest descent local search running on S
        """

        self.fitEval = 0
        start = time.time()
        self.model.transWal()
        self.oldindiv = individual.Individual(n=self.dim)
        self.oldindiv.init()
        self.oldindiv = self.evalPop(self.oldindiv)
        self.oldindiv.initWal(self.model)
        # self.model.initInter()
        self.bsf = individual.Individual(oldIndiv=self.oldindiv)
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

        traceEval = []
        traceFit = []

        initT = time.time() - start

        while self.fitEval < self.MaxFit:
            start = time.time()
            improveN, bestI = self.oldindiv.steepFitDesc(minimize)
            descT = descT + time.time() - start

            if improveN == False:
                initC = initC + 1
                if restart == True:
                    start = time.time()
                    diff, self.oldindiv = self.walk(fitName, minimize,False, walkLen, self.oldindiv)
                    pertT = pertT + time.time() - start


                    start = time.time()
                    for i in diff:
                        self.oldindiv.updateEval(i)
                        #                        self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[i]
                        # self.oldindiv.printSumArr()
                        self.oldindiv.update(i)
                        self.oldindiv.updateWAS(i)
                        self.oldindiv.updatePertImprS(i, minimize)
                    updatePertT = updatePertT + time.time() - start

                    # print 'bsf', self.bsf.fit
                    self.fitEval = self.fitEval + len(diff)

                    traceEval.append(self.fitEval)
                    traceFit.append(self.bsf.fit)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE
                start = time.time()
                # self.oldindiv.printSumArr()
                self.oldindiv.updateEval(bestI)
                self.oldindiv.update(bestI)
                self.oldindiv.updateWAS(bestI)
                self.oldindiv.updateImprS(bestI, minimize)
                self.fitEval = self.fitEval + 1
                updateT = updateT + time.time() - start
                updateC = updateC + 1
                self.oldindiv.flip(bestI)

                # print 'indiv', self.oldindiv.fit
                # self.tempIndiv = self.evalPop(self.oldindiv)
                # print 'real fit', self.tempIndiv.fit
                # print

                # if self.oldindiv.bit[bestI] == '1':
                #     self.oldindiv.bit[bestI] = '0'
                # else:
                #     self.oldindiv.bit[bestI] = '1'

        self.oldindiv.destructorWal(fitName)
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit, 'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'initC':initC, 'updateC':updateC, 'traceEval':traceEval, 'traceFit':traceFit}

    def fitCheck(self,fitName, minimize, restart):
        """
        steepest descent local search running on S, check the
        distribution of sum of Walsh terms in different orders
        """
        self.fitEval = 0
        start = time.time()
        self.model.transWal()

        self.oldindiv = individual.Individual(n=self.dim)
        self.oldindiv.init()
        self.oldindiv = self.evalPop(self.oldindiv)
        self.oldindiv.initWal(self.model)
        # self.model.initInter()
        self.bsf = individual.Individual(oldIndiv=self.oldindiv)
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
        traceEval = []
        traceFit = []
        initT = time.time() - start

        while self.fitEval < self.MaxFit:
            start = time.time()
            improveN, bestI = self.oldindiv.steepFitDesc(minimize)
            descT = descT + time.time() - start

            if improveN == False:
                initC = initC + 1
                if restart == True:
                    # check what traps local search before issuing the random walk
                    self.oldindiv.checkWalshSum()

                    start = time.time()
                    diff, self.oldindiv = self.walk(fitName, minimize,False, walkLen, self.oldindiv)
                    pertT = pertT + time.time() - start


                    start = time.time()
                    for i in diff:
                        self.oldindiv.updateEval(i)
                        self.oldindiv.update(i)
                        self.oldindiv.updateWAS(i)
                        self.oldindiv.updatePertImprS(i, minimize)
                    updatePertT = updatePertT + time.time() - start

                    self.fitEval = self.fitEval + len(diff)

                    traceEval.append(self.fitEval)
                    traceFit.append(self.bsf.fit)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE
                start = time.time()
                self.oldindiv.updateEval(bestI)
                self.oldindiv.update(bestI)
                self.oldindiv.updateWAS(bestI)
                self.oldindiv.updateImprS(bestI, minimize)
                self.fitEval = self.fitEval + 1
                updateT = updateT + time.time() - start
                updateC = updateC + 1
                self.oldindiv.flip(bestI)

        self.oldindiv.destructorWal(fitName)
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit, 'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'initC':initC, 'updateC':updateC, 'traceEval':traceEval, 'traceFit':traceFit}

    def runFitSwalkNext(self,fitName, minimize, restart):
        """
        next descent local search running on S
        """
        self.fitEval = 0
        start = time.time()
        self.model.transWal()

        self.oldindiv = individual.Individual(n=self.dim)
        self.oldindiv.init()
        self.oldindiv = self.evalPop(self.oldindiv)
        self.oldindiv.initWal(self.model)
        # self.model.initInter()
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

        initT = time.time() - start

        while self.fitEval < self.MaxFit:
            start = time.time()
            improveN, bestI = self.oldindiv.nextDesc()
            descT = descT + time.time() - start

            if improveN == False:
                initC = initC + 1
                if restart == True:

                    start = time.time()
                    diff, self.oldindiv = self.walk(fitName, minimize,False, walkLen, self.oldindiv)
                    pertT = pertT + time.time() - start

                    start = time.time()
                    for i in diff:
                        self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[i]
                        self.oldindiv.update(i)
                        self.oldindiv.updateWAS(i)
                        self.oldindiv.updatePertImprS(i, minimize)
                    updatePertT = updatePertT + time.time() - start

                    self.fitEval = self.fitEval + len(diff)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE
                start = time.time()
                self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[bestI]
                self.oldindiv.update(bestI)
                self.oldindiv.updateWAS(bestI)
                self.oldindiv.updateImprS(bestI, minimize)
                self.fitEval = self.fitEval + 1
                updateT = updateT + time.time() - start
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
        start = time.time()
        self.model.transWal()
        self.oldindiv = self.initIndiv(self.dim)
        self.oldindiv = self.evalPop(self.oldindiv)
        self.initWal(self.model)
        self.bsf = copy.deepcopy(self.oldindiv)
        self.model.WA = []

        walkLen = 10
        init = False
        updateT = 0
        initT = time.time() - start
        start = time.time()

##     def runFitS2walk(self,fitName, minimize, restart):
##         """
##         steepest descent local search ;running on S
##         """
##         self.fitEval = 0
##         start = time.time()
##         self.model.transWal()
##         self.oldindiv = self.initIndiv(self.dim)
##         self.oldindiv = self.evalPop(self.oldindiv)
##         self.initWal(self.model)
##         self.bsf = copy.deepcopy(self.oldindiv)
##         self.model.WA = []

##         walkLen = 10
##         init = False
##         updateT = 0
##         initT = time.time() - start
##         start = time.time()

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
##                     updateT = updateT + time.time() - start
##                     startR = time.time()
##                     self.oldindiv = self.evalPop(self.oldindiv)
##                     diff, self.oldindiv = self.walk(fitName, minimize,False, walkLen, self.oldindiv)
##                     init = False

## #                    print self.bsf.fit
## #                    pdb.set_trace()

##                     for i in diff:
##                         self.oldindiv.update(i)
##                         self.oldindiv.updateWAS(i)
##                     initT = initT + time.time() - startR
##                     start = time.time()
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
##         updateT = updateT + time.time() - start
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


    cdef runFitSrest(self,fitName, minimize, restart):
        """
        steepest descent local search running on S
        """
        # print 'runFitrest'
        self.fitEval = 0
        start = time.time()

        self.model.transWal()
        self.oldindiv = individual.Individual(n=self.dim)
        self.oldindiv.init()
        self.oldindiv = self.evalPop(self.oldindiv)
        self.oldindiv.initWal(self.model)
        # self.model.genInter()
        self.bsf = individual.Individual(oldIndiv=self.oldindiv)
        # print 'initWal'
        self.oldindiv.genImproveS(minimize)
        self.model.WA = []
        initC = 1
        updateC = 0
        # print 'init oldindiv', self.oldindiv.bit, self.oldindiv.fit
        # print 'init improve', self.oldindiv.improveA
        # print 'bit', self.oldindiv.bit, 'fit', self.oldindiv.fit
        # self.oldindiv.printSumArr()        
        
        descT = 0
        pertT = 0
        updatePertT = 0
        updateT = 0
        self.fitEval = 0
        walkLen = 10

        traceEval = []
        traceFit = []


        initT = time.time() - start

        while self.fitEval < self.MaxFit:
            start = time.time()
            improveN, bestI = self.oldindiv.steepFitDesc(minimize)
            # print 'bit', self.oldindiv.bit, 'fit', self.oldindiv.fit
            # self.oldindiv.printSumArr()
            # print 'bestI', bestI
            # print 'improveA', self.oldindiv.improveA
                        
            # print 'steep'
            descT = descT + time.time() - start
            # print 'oldindiv', self.oldindiv.bit, self.oldindiv.fit


            if improveN == False:
                initC = initC + 1
                if restart == True:
                    start = time.time()
                    oldbit = self.oldindiv.bit
                    oldfit = self.oldindiv.fit
                    self.restart(fitName, minimize, False)
                    # print 'after restart bit', self.oldindiv.bit
                    # print 'restart', 'bsf', self.bsf.fit, '\n'
                    
                    pertT = pertT + time.time() - start

                    start = time.time()
                    diff = self.diffBits(oldbit, self.oldindiv.bit)
                    self.oldindiv.fit = oldfit
                    for i in diff:
                        # self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[i]# TODO: need to count it in the next experiment
                        self.oldindiv.updateEval(i)
                        self.oldindiv.update(i)
                        self.oldindiv.updateWAS(i)
                        self.oldindiv.printSumArr()
                        self.oldindiv.updatePertImprS(i, minimize)
                    updatePertT = updatePertT + time.time() - start # TODO: need to count it in the next experiment
                    # self.fitEval = self.fitEval + len(diff) # TODO: need to count it in the next experiment
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE
                start = time.time()
                # self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[bestI]
                # print 'updateEval'
                self.oldindiv.updateEval(bestI)
                # print 'update'
                self.oldindiv.update(bestI)
                # print 'updateWAS'
                self.oldindiv.updateWAS(bestI)
                # print 'updateImprS'
                self.oldindiv.updateImprS(bestI, minimize)
                self.fitEval = self.fitEval + 1
                updateT = updateT + time.time() - start
                updateC = updateC + 1
                self.oldindiv.flip(bestI)
        # print 'dest'
        self.oldindiv.destructorWal(fitName)
        # print 'init', initC, 'update', updateC
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit, 'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'initC':initC, 'updateC':updateC, 'traceEval':traceEval, 'traceFit':traceFit}

    def runFitSrestFlip(self,fitName, minimize, restart):
        """
        steepest descent local search running on S, with tracing bit-flips enabled
        """
        self.fitEval = 0
        start = time.time()
        self.model.transWal()
        self.oldindiv = individual.Individual(n=self.dim)
        self.oldindiv.init()
        self.oldindiv = self.evalPop(self.oldindiv)
        self.oldindiv.initWal(self.model)
        # self.model.initInter()
        self.bsf = individual.Individual(oldIndiv=self.oldindiv)
        self.oldindiv.genImproveS(minimize)
        self.model.WA = []
        # print 'init', self.bsf.fit
        initC = 1
        updateC = 0

        # print 'init oldindiv', self.oldindiv.bit, self.oldindiv.fit
        # print 'init improve', self.oldindiv.improveA

        descT = 0
        pertT = 0
        updatePertT = 0
        updateT = 0
        self.fitEval = 0
        walkLen = 10

        traceEval = []
        traceFit = []
        traceFlip = np.zeros(self.dim)

        initT = time.time() - start

        while self.fitEval < self.MaxFit:
            # print 'fitEval', self.fitEval
            start = time.time()
            improveN, bestI = self.oldindiv.steepFitDesc(minimize)
            descT = descT + time.time() - start
            # print 'oldindiv', self.oldindiv.bit, self.oldindiv.fit

            if improveN == False:
                initC = initC + 1
                if restart == True:
                    start = time.time()
                    oldbit = self.oldindiv.bit
                    oldfit = self.oldindiv.fit
                    self.restart(fitName, minimize, False)
                    # print 'restart', 'bsf', self.bsf.fit, '\n'

                    pertT = pertT + time.time() - start

                    start = time.time()
                    diff = self.diffBits(oldbit, self.oldindiv.bit)
                    self.oldindiv.fit = oldfit
                    for i in diff:
                        # self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[i]# TODO: need to count it in the next experiment
                        self.oldindiv.updateEval(i)
                        self.oldindiv.update(i)
                        self.oldindiv.updateWAS(i)
                        self.oldindiv.updatePertImprS(i, minimize)
                    updatePertT = updatePertT + time.time() - start # TODO: need to count it in the next experiment
                    # self.fitEval = self.fitEval + len(diff) # TODO: need to count it in the next experiment
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE
                # print 'bestI', bestI
                # print 'improveA', self.oldindiv.improveA
                start = time.time()
                # self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[bestI]
                self.oldindiv.updateEval(bestI)
                self.oldindiv.update(bestI)
                self.oldindiv.updateWAS(bestI)
                self.oldindiv.updateImprS(bestI, minimize)
                self.fitEval = self.fitEval + 1
                updateT = updateT + time.time() - start
                updateC = updateC + 1
                self.oldindiv.flip(bestI)
                traceFlip[bestI] = traceFlip[bestI] + 1

        self.oldindiv.destructorWal(fitName)

        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit, 'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'initC':initC, 'updateC':updateC, 'traceEval':traceEval, 'traceFit':traceFit, 'traceFlip':traceFlip}


    def runFitSrestNext(self,fitName, minimize, restart):
        """
        next descent local search running on S
        """
        self.fitEval = 0
        start = time.time()
        self.model.transWal()
        self.oldindiv = individual.Individual(n=self.dim)
        self.oldindiv.init()
        self.oldindiv = self.evalPop(self.oldindiv)
        self.oldindiv.initWal(self.model)
        # self.model.initInter()
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

        initT = time.time() - start

        while self.fitEval < self.MaxFit:
            start = time.time()
            improveN, bestI = self.oldindiv.nextDesc()
            descT = descT + time.time() - start

            if improveN == False:
                initC = initC + 1
                if restart == True:
                    start = time.time()
                    oldbit = self.oldindiv.bit
                    oldfit = self.oldindiv.fit
                    self.restart(fitName, minimize, False)
                    pertT = pertT + time.time() - start

                    start = time.time()
                    diff = self.diffBits(oldbit, self.oldindiv.bit)
                    self.oldindiv.fit = oldfit
                    for i in diff:
                        self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[i]
                        self.oldindiv.update(i)
                        self.oldindiv.updateWAS(i)
                        self.oldindiv.updatePertImprS(i, minimize)
                    updatePertT = updatePertT + time.time() - start

                    self.fitEval = self.fitEval + len(diff)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE
                start = time.time()

                self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[bestI]

                self.oldindiv.update(bestI)
                self.oldindiv.updateWAS(bestI)
                self.oldindiv.updateImprS(bestI, minimize)
                self.fitEval = self.fitEval + 1
                updateT = updateT + time.time() - start
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

        # print('init start')

        start = time.time()
        self.model.transWal()
        # self.model.initInter()

        self.oldindiv = individual.Individual(n=self.dim,neigh=True)
        self.oldindiv.init()
        self.oldindiv = self.evalPop(self.oldindiv)
        self.oldindiv.initWal(self.model)
        self.oldindiv.initSC()
        self.oldindiv.compFitG()
        self.oldindiv.genImproveSC(minimize)

        # self.bsf = copy.deepcopy(self.oldindiv)
        self.bsf = individual.Individual(oldIndiv=self.oldindiv)
        self.model.WA = []
        # print 'bsf', self.bsf.bit, self.bsf.fit, self.bsf.fitG
        # print('init end')

        initC = 0
        updateC = 0

        descT = 0
        pertT = 0
        updatePertT = 0
        updateT = 0
        self.fitEval = 0
        walkLen = 10

        traceEval = []
        traceFit = []
        traceFitG = []

        initT = time.time() - start
        while self.fitEval < self.MaxFit:
            start = time.time()
            improveN, bestI = self.oldindiv.steepMeanDesc(minimize)
            descT = descT + time.time() - start

            if improveN == False:
                initC = initC + 1
                if restart == True:

                    start = time.time()
                    diff, self.oldindiv = self.walk(fitName, minimize, False, walkLen, self.oldindiv)
                    pertT = pertT + time.time() - start


                    start = time.time()
                    # print('walk begin')
                    for i in diff:
                        # self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[i]
                        # self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))
                        self.oldindiv.updateEval(i)
                        self.oldindiv.update(i)
                        self.oldindiv.updateSC(i)
                        self.oldindiv.updateWAS(i)
                        self.oldindiv.updatePertImprSC(i, minimize)
                        self.oldindiv.compFitG()

                    # print 'bsf', self.bsf.bit, self.bsf.fit, self.bsf.fitG
                    updatePertT = updatePertT + time.time() - start
                    # print('walk end')
                    self.fitEval = self.fitEval + len(diff)

                    traceEval.append(self.fitEval)
                    traceFit.append(self.bsf.fit)
                    traceFitG.append(self.bsf.fitG)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE
                start = time.time()

                # print('flip begin')
                # self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[bestI]
                # self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))
                self.oldindiv.updateEval(bestI)
                self.oldindiv.update(bestI)
                self.oldindiv.updateSC(bestI)
                self.oldindiv.updateWAS(bestI)
                self.oldindiv.updateImprSC(bestI, minimize)
                self.fitEval = self.fitEval + 1
                self.oldindiv.compFitG()
                self.oldindiv.flip(bestI)
                # print('flip end')
                updateT = updateT + time.time() - start
                updateC = updateC + 1

                # print 'indiv', self.oldindiv.bit, self.oldindiv.fit, self.oldindiv.fitG
                # self.tempIndiv = self.evalPopNeigh(self.oldindiv, fitName, minimize)
                # print 'real fit', self.oldindiv.bit, self.tempIndiv.fit, self.tempIndiv.fitG
                # print

                # if self.oldindiv.bit[bestI] == '1':
                #     self.oldindiv.bit[bestI] = '0'
                # else:
                #     self.oldindiv.bit[bestI] = '1'

        # print('dest begin')
        self.oldindiv.destructorWal(fitName)
        # print('dest end')
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'fitG': self.bsf.fitG, 'bit':self.bsf.bit,'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'initC':initC, 'updateC':updateC, 'traceEval':traceEval, 'traceFit':traceFit, 'traceFitG':traceFitG}

    def runSwitchFwalk(self,fitName, minimize, restart):
        """
        steepest descent local search with respect to mean of neighs by Walsh Analysis
        """
        self.fitEval = 0

        start = time.time()
        self.model.transWal()
        # self.model.initInter()

        self.oldindiv = individual.Individual(n=self.dim,neigh=True)
        self.oldindiv.init()
        self.oldindiv = self.evalPop(self.oldindiv)
        self.oldindiv.initWal(self.model)
        self.oldindiv.initSC()
        self.oldindiv.compFitG()
        self.oldindiv.genImproveS(minimize)
        self.oldindiv.genImproveSC(minimize)

        # self.bsf = copy.deepcopy(self.oldindiv)
        self.bsf = individual.Individual(oldIndiv=self.oldindiv)
        self.model.WA = []
        # print 'bsf', self.bsf.bit, self.bsf.fit, self.bsf.fitG
        # print('init end')

        platC = 0                       # hit a plateau in the original space
        restC = 0                       # issue a restart when both of method fail
        updateC = 0                     # the number ascent steep taken by LS

        descT = 0
        pertT = 0
        updatePertT = 0
        updateT = 0
        self.fitEval = 0
        walkLen = 10

        traceEval = []
        traceFit = []
        traceFitG = []

        initT = time.time() - start
        while self.fitEval < self.MaxFit:
            start = time.time()
            improveN, bestI = self.oldindiv.steepFitDesc(minimize)
            descT = descT + time.time() - start

            if improveN == False:
                platC = platC + 1
                improveMeanN, bestI = self.oldindiv.steepMeanDesc(minimize)
            if improveN==False and improveMeanN == False:
                restC = restC + 1
                if restart == True:

                    start = time.time()
                    diff, self.oldindiv = self.walk(fitName, minimize, False, walkLen, self.oldindiv)
                    pertT = pertT + time.time() - start


                    start = time.time()
                    # print('walk begin')
                    for i in diff:
                        # self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[i]
                        # self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))
                        self.oldindiv.updateEval(i)
                        self.oldindiv.update(i)
                        self.oldindiv.updateSC(i)
                        self.oldindiv.updateWAS(i)
                        self.oldindiv.updatePertImprS(i, minimize)
                        self.oldindiv.updatePertImprSC(i, minimize)
                        self.oldindiv.compFitG()

                    # print 'bsf', self.bsf.bit, self.bsf.fit, self.bsf.fitG
                    updatePertT = updatePertT + time.time() - start
                    # print('walk end')
                    self.fitEval = self.fitEval + len(diff)

                    traceEval.append(self.fitEval)
                    traceFit.append(self.bsf.fit)
                    traceFitG.append(self.bsf.fitG)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE
                start = time.time()

                # print('flip begin')
                # self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[bestI]
                # self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))
                self.oldindiv.updateEval(bestI)
                self.oldindiv.update(bestI)
                self.oldindiv.updateSC(bestI)
                self.oldindiv.updateWAS(bestI)
                self.oldindiv.updateImprS(bestI, minimize)
                self.oldindiv.updateImprSC(bestI, minimize)
                self.fitEval = self.fitEval + 1
                self.oldindiv.compFitG()
                self.oldindiv.flip(bestI)
                # print('flip end')
                updateT = updateT + time.time() - start
                updateC = updateC + 1

                # print 'indiv', self.oldindiv.bit, self.oldindiv.fit, self.oldindiv.fitG
                # self.tempIndiv = self.evalPopNeigh(self.oldindiv, fitName, minimize)
                # print 'real fit', self.oldindiv.bit, self.tempIndiv.fit, self.tempIndiv.fitG
                # print

                # if self.oldindiv.bit[bestI] == '1':
                #     self.oldindiv.bit[bestI] = '0'
                # else:
                #     self.oldindiv.bit[bestI] = '1'

        # print('dest begin')
        self.oldindiv.destructorWal(fitName)
        # print('dest end')
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'fitG': self.bsf.fitG, 'bit':self.bsf.bit,'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'platC':platC, 'restC':restC, 'updateC':updateC, 'traceEval':traceEval, 'traceFit':traceFit, 'traceFitG':traceFitG}

    def runCombFwalk(self,fitName, minimize, restart):
        self.fitEval = 0

        start = time.time()
        self.model.transWal()
        # self.model.initInter()

        self.oldindiv = individual.Individual(n=self.dim,neigh=True)
        self.oldindiv.init()
        self.oldindiv = self.evalPop(self.oldindiv)
        self.oldindiv.initWal(self.model)
        self.oldindiv.initSC()
        self.oldindiv.compFitG()
        self.oldindiv.genImproveS(minimize)
        self.oldindiv.genImproveSC(minimize)

        self.bsf = individual.Individual(oldIndiv=self.oldindiv)
        self.model.WA = []

        platC = 0                       # hit a plateau in the original space
        restC = 0                       # issue a restart when both of method fail
        updateC = 0                     # the number ascent steep taken by LS

        descT = 0
        pertT = 0
        updatePertT = 0
        updateT = 0
        self.fitEval = 0
        walkLen = 10

        useF = True
        switch = False

        traceEval = []
        traceFit = []
        traceFitG = []

        initT = time.time() - start
        while self.fitEval < self.MaxFit:
            start = time.time()
            if useF == True:
                improveN, bestI = self.oldindiv.steepFitDesc(minimize)
            else:
                improveN, bestI = self.oldindiv.steepMeanDesc(minimize)
            switch=False
            # print 'improveN', improveN
            # print 'useF', useF
            # print 'switch', switch
            # print
            descT = descT + time.time() - start

            if improveN == False:
                platC = platC + 1
                useF = not useF
                if useF == True:
                    improveN, bestI = self.oldindiv.steepFitDesc(minimize)
                else:
                    improveN, bestI = self.oldindiv.steepMeanDesc(minimize)

                switch = True
                # print 'platC', platC
                # print 'improveN', improveN
                # print 'switch', switch
                # print

            if improveN==False and switch==True:
                restC = restC + 1
                # print 'restC', restC
                # print 'useF', useF

                if restart == True:

                    start = time.time()
                    diff, self.oldindiv = self.walk(fitName, minimize, False, walkLen, self.oldindiv)
                    pertT = pertT + time.time() - start


                    start = time.time()
                    # print('walk begin')
                    for i in diff:
                        # self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[i]
                        # self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))
                        self.oldindiv.updateEval(i)
                        self.oldindiv.update(i)
                        self.oldindiv.updateSC(i)
                        self.oldindiv.updateWAS(i)
                        self.oldindiv.updatePertImprS(i, minimize)
                        self.oldindiv.updatePertImprSC(i, minimize)
                        self.oldindiv.compFitG()

                    # print 'bsf', self.bsf.bit, self.bsf.fit, self.bsf.fitG
                    updatePertT = updatePertT + time.time() - start
                    # print('walk end')
                    self.fitEval = self.fitEval + len(diff)

                    traceEval.append(self.fitEval)
                    traceFit.append(self.bsf.fit)
                    traceFitG.append(self.bsf.fitG)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE
                start = time.time()

                # print('flip begin')
                # self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[bestI]
                # self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))
                self.oldindiv.updateEval(bestI)
                self.oldindiv.update(bestI)
                self.oldindiv.updateSC(bestI)
                self.oldindiv.updateWAS(bestI)
                self.oldindiv.updateImprS(bestI, minimize)
                self.oldindiv.updateImprSC(bestI, minimize)
                self.fitEval = self.fitEval + 1
                self.oldindiv.compFitG()
                self.oldindiv.flip(bestI)
                # print('flip end')
                updateT = updateT + time.time() - start
                updateC = updateC + 1

                # print 'indiv', self.oldindiv.bit, self.oldindiv.fit, self.oldindiv.fitG
                # self.tempIndiv = self.evalPopNeigh(self.oldindiv, fitName, minimize)
                # print 'real fit', self.oldindiv.bit, self.tempIndiv.fit, self.tempIndiv.fitG
                # print

                # if self.oldindiv.bit[bestI] == '1':
                #     self.oldindiv.bit[bestI] = '0'
                # else:
                #     self.oldindiv.bit[bestI] = '1'

        # print('dest begin')
        self.oldindiv.destructorWal(fitName)
        # print('dest end')
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'fitG': self.bsf.fitG, 'bit':self.bsf.bit,'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'platC':platC, 'restC':restC, 'updateC':updateC, 'traceEval':traceEval, 'traceFit':traceFit, 'traceFitG':traceFitG}

    def runSwitchAvgwalk(self,fitName, minimize, restart):
        """
        steepest descent local search with respect to mean of neighs by Walsh Analysis
        """
        self.fitEval = 0

        start = time.time()
        self.model.transWal()
        # self.model.initInter()

        self.oldindiv = individual.Individual(n=self.dim,neigh=True)
        self.oldindiv.init()
        self.oldindiv = self.evalPop(self.oldindiv)
        self.oldindiv.initWal(self.model)
        self.oldindiv.initSC()
        self.oldindiv.compFitG()
        self.oldindiv.genImproveS(minimize)
        self.oldindiv.genImproveSC(minimize)

        # self.bsf = copy.deepcopy(self.oldindiv)
        self.bsf = individual.Individual(oldIndiv=self.oldindiv)
        self.model.WA = []
        # print 'bsf', self.bsf.bit, self.bsf.fit, self.bsf.fitG
        # print('init end')

        platC = 0                       # hit a plateau in the original space
        restC = 0                       # issue a restart when both of method fail
        updateC = 0                     # the number ascent steep taken by LS

        descT = 0
        pertT = 0
        updatePertT = 0
        updateT = 0
        self.fitEval = 0
        walkLen = 10

        traceEval = []
        traceFit = []
        traceFitG = []

        initT = time.time() - start
        while self.fitEval < self.MaxFit:
            start = time.time()
            improveMeanN, bestI = self.oldindiv.steepMeanDesc(minimize)

            descT = descT + time.time() - start

            if improveMeanN == False:
                # if improveMeanN==True:
                platC = platC + 1
                improveN, bestI = self.oldindiv.steepFitDesc(minimize)
                #     print improveMeanN
            if  improveMeanN == False and improveN==False:
                restC = restC + 1
                if restart == True:

                    start = time.time()
                    diff, self.oldindiv = self.walk(fitName, minimize, False, walkLen, self.oldindiv)
                    pertT = pertT + time.time() - start

                    start = time.time()
                    # print('walk begin')
                    for i in diff:
                        # self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[i]
                        # self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))
                        self.oldindiv.updateEval(i)
                        self.oldindiv.update(i)
                        self.oldindiv.updateSC(i)
                        self.oldindiv.updateWAS(i)
                        self.oldindiv.updatePertImprS(i, minimize)
                        self.oldindiv.updatePertImprSC(i, minimize)
                        self.oldindiv.compFitG()

                    # print 'bsf', self.bsf.bit, self.bsf.fit, self.bsf.fitG
                    updatePertT = updatePertT + time.time() - start
                    # print('walk end')
                    self.fitEval = self.fitEval + len(diff)

                    traceEval.append(self.fitEval)
                    traceFit.append(self.bsf.fit)
                    traceFitG.append(self.bsf.fitG)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE
                start = time.time()

                # print('flip begin')
                # self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[bestI]
                # self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))
                self.oldindiv.updateEval(bestI)
                self.oldindiv.update(bestI)
                self.oldindiv.updateSC(bestI)
                self.oldindiv.updateWAS(bestI)
                self.oldindiv.updateImprS(bestI, minimize)
                self.oldindiv.updateImprSC(bestI, minimize)
                self.fitEval = self.fitEval + 1
                self.oldindiv.compFitG()
                self.oldindiv.flip(bestI)
                # print('flip end')
                updateT = updateT + time.time() - start
                updateC = updateC + 1

                # print 'indiv', self.oldindiv.bit, self.oldindiv.fit, self.oldindiv.fitG
                # self.tempIndiv = self.evalPopNeigh(self.oldindiv, fitName, minimize)
                # print 'real fit', self.oldindiv.bit, self.tempIndiv.fit, self.tempIndiv.fitG
                # print

                # if self.oldindiv.bit[bestI] == '1':
                #     self.oldindiv.bit[bestI] = '0'
                # else:
                #     self.oldindiv.bit[bestI] = '1'

        # print('dest begin')
        self.oldindiv.destructorWal(fitName)
        # print('dest end')
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'fitG': self.bsf.fitG, 'bit':self.bsf.bit,'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'platC':platC, 'restC':restC, 'updateC':updateC, 'traceEval':traceEval, 'traceFit':traceFit, 'traceFitG':traceFitG}

    def runCombAvgwalk(self,fitName, minimize, restart):
        self.fitEval = 0
        start = time.time()
        self.model.transWal()
        # self.model.initInter()

        self.oldindiv = individual.Individual(n=self.dim,neigh=True)
        self.oldindiv.init()
        self.oldindiv = self.evalPop(self.oldindiv)
        self.oldindiv.initWal(self.model)
        self.oldindiv.initSC()
        self.oldindiv.compFitG()
        self.oldindiv.genImproveS(minimize)
        self.oldindiv.genImproveSC(minimize)

        self.bsf = individual.Individual(oldIndiv=self.oldindiv)
        self.model.WA = []

        platC = 0                       # hit a plateau in the original space
        restC = 0                       # issue a restart when both of method fail
        updateC = 0                     # the number ascent steep taken by LS

        descT = 0
        pertT = 0
        updatePertT = 0
        updateT = 0
        self.fitEval = 0
        walkLen = 10

        useF = False
        switch = False

        traceEval = []
        traceFit = []
        traceFitG = []

        initT = time.time() - start
        while self.fitEval < self.MaxFit:
            start = time.time()
            if useF == True:
                improveN, bestI = self.oldindiv.steepFitDesc(minimize)
            else:
                improveN, bestI = self.oldindiv.steepMeanDesc(minimize)
            switch=False
            # print 'improveN', improveN
            # print 'useF', useF
            # print 'switch', switch
            # print
            descT = descT + time.time() - start

            if improveN == False:
                platC = platC + 1
                useF = not useF
                if useF == True:
                    improveN, bestI = self.oldindiv.steepFitDesc(minimize)
                else:
                    improveN, bestI = self.oldindiv.steepMeanDesc(minimize)

                switch = True
                # print 'platC', platC
                # print 'improveN', improveN
                # print 'switch', switch
                # print

            if improveN==False and switch==True:
                restC = restC + 1
                # print 'restC', restC
                # print 'useF', useF

                if restart == True:

                    start = time.time()
                    diff, self.oldindiv = self.walk(fitName, minimize, False, walkLen, self.oldindiv)
                    pertT = pertT + time.time() - start


                    start = time.time()
                    # print('walk begin')
                    for i in diff:
                        # self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[i]
                        # self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))
                        self.oldindiv.updateEval(i)
                        self.oldindiv.update(i)
                        self.oldindiv.updateSC(i)
                        self.oldindiv.updateWAS(i)
                        self.oldindiv.updatePertImprS(i, minimize)
                        self.oldindiv.updatePertImprSC(i, minimize)
                        self.oldindiv.compFitG()

                    # print 'bsf', self.bsf.bit, self.bsf.fit, self.bsf.fitG
                    updatePertT = updatePertT + time.time() - start
                    # print('walk end')
                    self.fitEval = self.fitEval + len(diff)

                    traceEval.append(self.fitEval)
                    traceFit.append(self.bsf.fit)
                    traceFitG.append(self.bsf.fitG)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE
                start = time.time()

                # print('flip begin')
                # self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[bestI]
                # self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))
                self.oldindiv.updateEval(bestI)
                self.oldindiv.update(bestI)
                self.oldindiv.updateSC(bestI)
                self.oldindiv.updateWAS(bestI)
                self.oldindiv.updateImprS(bestI, minimize)
                self.oldindiv.updateImprSC(bestI, minimize)
                self.fitEval = self.fitEval + 1
                self.oldindiv.compFitG()
                self.oldindiv.flip(bestI)
                # print('flip end')
                updateT = updateT + time.time() - start
                updateC = updateC + 1

                # print 'indiv', self.oldindiv.bit, self.oldindiv.fit, self.oldindiv.fitG
                # self.tempIndiv = self.evalPopNeigh(self.oldindiv, fitName, minimize)
                # print 'real fit', self.oldindiv.bit, self.tempIndiv.fit, self.tempIndiv.fitG
                # print

                # if self.oldindiv.bit[bestI] == '1':
                #     self.oldindiv.bit[bestI] = '0'
                # else:
                #     self.oldindiv.bit[bestI] = '1'

        # print('dest begin')
        self.oldindiv.destructorWal(fitName)
        # print('dest end')
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'fitG': self.bsf.fitG, 'bit':self.bsf.bit,'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'platC':platC, 'restC':restC, 'updateC':updateC, 'traceEval':traceEval, 'traceFit':traceFit, 'traceFitG':traceFitG}


    def runMeanSCwalkNext(self,fitName, minimize, restart):
        """
        next descent local search with respect to mean of neighs by Walsh Analysis
        """
        self.fitEval = 0

        start = time.time()
        self.model.transWal()
        self.oldindiv = individual.Individual(n=self.dim,neigh=True)
        self.oldindiv.init()
        self.oldindiv = self.evalPop(self.oldindiv)
        # self.model.initInter()
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

        initT = time.time() - start
        while self.fitEval < self.MaxFit:
            start = time.time()
            improveN, bestI = self.oldindiv.nextDesc()
            descT = descT + time.time() - start

            if improveN == False:
                initC = initC + 1
                if restart == True:

                    start = time.time()
                    diff, self.oldindiv = self.walk(fitName, minimize,False, walkLen, self.oldindiv)
                    pertT = pertT + time.time() - start

                    start = time.time()
                    for i in diff:
                        self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[i]
                        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))

                        self.oldindiv.update(i)
                        self.oldindiv.updateSC(i)
                        self.oldindiv.updateWAS(i)
                        self.oldindiv.updatePertImprSC(i, minimize)
                    updatePertT = updatePertT + time.time() - start

                    self.fitEval = self.fitEval + len(diff)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE
                start = time.time()
                self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[bestI]
                self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))
                self.oldindiv.update(bestI)
                self.oldindiv.updateSC(bestI)
                self.oldindiv.updateWAS(bestI)
                self.oldindiv.updateImprSC(bestI, minimize)
                self.fitEval = self.fitEval + 1
                updateT = updateT + time.time() - start
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

        start = time.time()
        self.model.transWal()
        self.oldindiv = individual.Individual(n=self.dim,neigh=True)
        self.oldindiv.init()
        self.oldindiv = self.evalPop(self.oldindiv)
        # self.model.initInter()
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

        initT = time.time() - start
        while self.fitEval < self.MaxFit:
            start = time.time()
            improveN, bestI = self.oldindiv.steepMeanDesc(minimize)
            descT = descT + time.time() - start

            if improveN == False:
                initC = initC + 1
                if restart == True:
                    start = time.time()
                    oldbit = self.oldindiv.bit
                    oldfit = self.oldindiv.fit
                    self.restart(fitName, minimize, False)
                    pertT = pertT + time.time() - start

                    start = time.time()
                    diff = self.diffBits(oldbit, self.oldindiv.bit)
                    self.oldindiv.fit = oldfit
                    for i in diff:
                        self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[i]
                        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))

                        self.oldindiv.update(i)
                        self.oldindiv.updateSC(i)
                        self.oldindiv.updateWAS(i)
                        self.oldindiv.updatePertImprS(i, minimize)
                    updatePertT = updatePertT + time.time() - start

                    self.fitEval = self.fitEval + len(diff)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE
                start = time.time()

                self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[bestI]
                self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))

                self.oldindiv.update(bestI)
                self.oldindiv.updateSC(bestI)
                self.oldindiv.updateWAS(bestI)
                self.oldindiv.updateImprSC(bestI, minimize)
                self.fitEval = self.fitEval + 1
                updateT = updateT + time.time() - start
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

        start = time.time()
        self.model.transWal()
        self.oldindiv = individual.Individual(n=self.dim,neigh=True)
        self.oldindiv.init()
        self.oldindiv = self.evalPop(self.oldindiv)
        # self.model.initInter()
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

        initT = time.time() - start
        while self.fitEval < self.MaxFit:
            start = time.time()
            improveN, bestI = self.oldindiv.nextDesc()
            descT = descT + time.time() - start

            if improveN == False:
                initC = initC + 1
                if restart == True:
                    start = time.time()
                    oldbit = self.oldindiv.bit
                    oldfit = self.oldindiv.fit
                    self.restart(fitName, minimize, False)
                    pertT = pertT + time.time() - start

                    start = time.time()
                    diff = self.diffBits(oldbit, self.oldindiv.bit)
                    self.oldindiv.fit = oldfit
                    for i in diff:
                        self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[i]
                        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))

                        self.oldindiv.update(i)
                        self.oldindiv.updateSC(i)
                        self.oldindiv.updateWAS(i)
                        self.oldindiv.updatePertImprS(i, minimize)
                    updatePertT = updatePertT + time.time() - start

                    self.fitEval = self.fitEval + len(diff)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE
                start = time.time()

                self.oldindiv.fit = self.oldindiv.fit - 2*self.oldindiv.sumArr[bestI]
                self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))

                self.oldindiv.update(bestI)
                self.oldindiv.updateSC(bestI)
                self.oldindiv.updateWAS(bestI)
                self.oldindiv.updateImprSC(bestI, minimize)
                self.fitEval = self.fitEval + 1
                updateT = updateT + time.time() - start
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

        start = time.time()
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
        initT = time.time() - start
        start = time.time()
        while self.fitEval < self.MaxFit:
            if init == False:
                improveN, bestI, evalCount = self.genMeanBest(minimize)
                init = True
            else :
                improveN, bestI, evalCount = self.oldindiv.updateMeanBest(bestI,minimize)
            self.fitEval = self.fitEval + self.dim

            if improveN == False:
                if restart == True:
                    updateT = updateT + time.time() - start
                    startR = time.time()
                    self.oldindiv = self.evalPop(self.oldindiv)
                    self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.oldindiv.sumArr))

                    diff, self.oldindiv = self.walk( fitName, minimize, False, walkLen , self.oldindiv, self.oldindiv)
                    init = False

                    for i in diff:
                        self.oldindiv.update(i)

                        self.oldindiv.updateSC(i)
                        self.oldindiv.updateWAS(i)
                    initT = initT + time.time() - startR
                    start = time.time()
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
        updateT = updateT + time.time() - start
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'fitG': self.bsf.fitG, 'bit':self.bsf.bit,'init':initT, 'update':updateT}

    cdef runFit(self, minimize,restart):
        """
        brute force approach for running BILS
        """
        cdef list improveA
        # print 'init'
        # self.oldindiv = self.initIndiv(self.dim)
        self.oldindiv = individual.Individual(n=self.dim)
        self.oldindiv.init()
        self.fitEval = 0
        # print 'init'
        self.oldindiv = self.evalPop(self.oldindiv)
        # print 'init'
        # self.bsf = copy.deepcopy(self.oldindiv)
        self.bsf = individual.Individual(oldIndiv=self.oldindiv)
        # self.trace = [Struct(fitEval= self.fitEval,fit = self.oldindiv.fit)]
        initC = 1
        updateC = 0
        # print 'init oldindiv', self.oldindiv.bit, self.oldindiv.fit

        updateT = 0

        start = time.time()

        while self.fitEval < self.MaxFit:

            # print 'fitEval', self.fitEval
            # neighs = self.neighbors()
            #print
            # print 'current', self.oldindiv.bit, 'fit', self.oldindiv.fit
            improveA = []
            bestFit = self.oldindiv.fit   # keep track of the best fitness
            # print 'oldindiv', self.oldindiv.bit, self.oldindiv.fit
            for i in xrange(self.dim):
                self.indiv = individual.Individual(oldIndiv=self.oldindiv)
                self.indiv.flip(i)
                self.indiv = self.evalPop(self.indiv)
                # print self.indiv.bit
            # for i in neighs:
            #     self.indiv.bit = np.copy(i)
            #     self.indiv = self.evalPop(self.indiv)
            #     #print 'neigh: ', self.indiv.bit, 'fit', self.indiv.fit

                improveA, bestFit = self.selectionFit(minimize, improveA, bestFit, i)
                # print  'i',i,'best', bestFit, 'current', self.indiv.fit, 'bsf', self.bsf.fit, 'len', len(improveA)

                # if self.fitEval>=self.MaxFit:
                #     break

#            self.trace.append(Struct(fitEval= self.fitEval,fit = self.oldindiv.fit))
            if len(improveA)==0: # issue restart
                updateT = updateT + time.time() - start
                initC = initC + 1
                if restart == True:
                    self.restart('fit', minimize, True)
                    # print 'restart', 'bsf', self.bsf.fit, '\n'
                else:
                    return {'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit}
                start = time.time()
            else: # randomly take an best-improvement move
                updateC = updateC + 1
                # print 'improve'
                bestI = random.choice(improveA)
                # print 'bestI', bestI, 'eval', self.fitEval
                # print 'improveA', improveA
                self.oldindiv.flip(bestI)
                self.oldindiv = self.evalPop(self.oldindiv)
                # print self.oldindiv.fit
                self.fitEval = self.fitEval + 1

        updateT = updateT + time.time() - start

        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit, 'initC':initC, 'updateC':updateC,  'updateT':updateT}



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

    cdef runFitUpdate(self, fitName, minimize, restart):
        # TODO: change the function type cdef
        """
        run local search using fit function, with paritial update
        """
        # cdef double* subFitArr
        # print 'runFitUpdate'
        self.fitEval = 0
        start = time.time()
        self.oldindiv = individual.Individual(n=self.dim)
        self.oldindiv.init()

        # print 'genInter'
        self.model.genInter()
        self.model.genListSubFunc()

        # self.model.printListSubFunc()
        
        self.oldindiv = self.evalPop(self.oldindiv)
        # print ''
        # print 'genU'
        self.model.genU() 
        self.oldindiv.initBfUpdate(self.oldindiv, self.evalPop, minimize, self.model)
        self.bsf = individual.Individual(oldIndiv=self.oldindiv)
        # print 'bit', self.oldindiv.bit, 'fit', self.oldindiv.fit
        # print 'init SumArr'
        self.oldindiv.printSumArr()
        
        while self.fitEval < self.MaxFit:
            # print 'desc eval', self.fitEval, 'max', self.MaxFit
            improveN, bestI = self.oldindiv.steepFitDescPartialUpdate(minimize)
            # print 'bit', self.oldindiv.bit, 'fit', self.oldindiv.fit
            # self.oldindiv.printSumArr()
            # print 'bestI', bestI
            # print 'improveA', self.oldindiv.improveA
                            

            if improveN == False:
                if restart == True:
                    oldbit = self.oldindiv.bit
                    oldfit = self.oldindiv.fit
                    self.restart(fitName, minimize, False)
                    # print 'after restart bit', self.oldindiv.bit
                    # print 'restart', 'bsf', self.bsf.fit, '\n'
                    
                    diff = self.diffBits(oldbit, self.oldindiv.bit)
                    self.oldindiv.fit = oldfit
                    self.oldindiv.bit = oldbit
                    # print 'old bit', self.oldindiv.bit, self.oldindiv.fit
                    for i in diff:
                        # take the move
                        # print 'eval'
                        self.oldindiv.updateEvalPartialUpdate(i)
                        # updating data structures to reflect current solution, namely the first derivate 
                        # print 'sumArr'
                        self.oldindiv.updateSumArr(i, self.oldindiv)
                        

                        # print 'pertimprs'
                        self.oldindiv.updatePertImprSpartialUpdate(i, minimize)
                        self.oldindiv.flip(i)
                        self.oldindiv.printSumArr()
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit}
            else:
                # print 'eval'
                self.oldindiv.updateEvalPartialUpdate(bestI)
                # print 'sumArr'
                self.oldindiv.updateSumArr(bestI, self.oldindiv)
                # print 'imprs'
                self.oldindiv.updateImprSpartialUpdate(bestI, minimize)
                self.fitEval = self.fitEval + 1
                self.oldindiv.flip(bestI)
                
        self.oldindiv.destructorBfUpdate()
        #        print 'return'
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit}

        
        

    cdef runBeamFitSwalkNext(self,fitName, minimize, restart, beamWidth):
        """
        BEAM next descent local search running on S,
        where beamWidth is the number of beams
        """
        self.fitEval = 0
        start = time.time()
        self.model.transWal()
        # self.model.initInter()
        self.oldpop = []
        for i in xrange(beamWidth):
            self.oldpop.append(individual.Individual(n=self.dim))
            self.oldpop[i].init()
            self.oldpop[i] = self.evalPop(self.oldpop[i])
            self.oldpop[i].initWal(self.model)
            self.oldpop[i].genImproveS(minimize)

        # self.bsf is a single individual recording the best-so-far solution
        #self.bsf = copy.deepcopy(self.oldpop[self.bestIndiv(self.oldpop, minimize, fitName)])
        self.bsf = individual.Individual(oldIndiv=self.oldpop[self.bestIndiv(self.oldpop, minimize, fitName)])
        #print self.bsf.fit
        self.model.WA = [] # release memory

        initC = 1
        updateC = 0
        descT = 0
        pertT = 0
        updatePertT = 0
        updateT = 0
        self.fitEval = 0
        walkLen = 10

        initT = time.time() - start

        while self.fitEval < self.MaxFit:
            for i in xrange(len(self.oldpop)):
                start = time.time()
                improveN, bestI = self.oldpop[i].nextDesc()
                descT = descT + time.time() - start
#                print 'bestI', bestI
#                print 'improveA', self.oldpop[i].improveA


                if improveN == False:
                    initC = initC + 1
                    if restart == True:
                        start = time.time()
                        diff, self.oldpop[i] = self.walk(fitName, minimize, False, walkLen, self.oldpop[i])
#                        print 'walk',self.bsf.fit, self.fitEval

                        pertT = pertT + time.time() - start

                        start = time.time()
                        for j in diff:
                            #                            self.oldpop[i].fit = self.oldpop[i].fit - 2*self.oldpop[i].sumArr[j]
                            self.oldpop[i].updateEval(j)
                            self.oldpop[i].update(j)
                            self.oldpop[i].updateWAS(j)
                            self.oldpop[i].updatePertImprS(j, minimize)
                        updatePertT = updatePertT + time.time() - start

                        self.fitEval = self.fitEval + len(diff)
                    else:
                        return { 'nEvals': self.fitEval, 'sol': self.oldpop[i].fit, 'bit':self.oldpop[i].bit}
                else : # improveN is TRUE
                    start = time.time()
                    #                    self.oldpop[i].fit = self.oldpop[i].fit - 2*self.oldpop[i].sumArr[bestI]
                    self.oldpop[i].updateEval(bestI)
                    self.oldpop[i].update(bestI)
                    self.oldpop[i].updateWAS(bestI)
                    self.oldpop[i].updateImprS(bestI, minimize)
                    self.fitEval = self.fitEval + 1
#                    print 'decs', self.bsf.fit, self.fitEval
                    updateT = updateT + time.time() - start
                    updateC = updateC + 1
                    self.oldpop[i].flip(bestI)
                    ## if self.oldpop[i].bit[bestI] == '1':
                    ##     self.oldpop[i].bit[bestI] = '0'
                    ## else:
                    ##     self.oldpop[i].bit[bestI] =
        for i in xrange(len(self.oldpop)):
            self.oldpop[i].destructorWal(fitName)
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit, 'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'initC':initC, 'updateC':updateC}

    def runBeamFitSwalk(self,fitName, minimize, restart, beamWidth):
        """
        BEAM steepest descent local search running on S,
        where beamWidth is the number of beams
        """
        self.fitEval = 0
        start = time.time()
        self.model.transWal()
        # self.model.initInter()
        self.oldpop = []
        for i in xrange(beamWidth):
            self.oldpop.append(individual.Individual(n=self.dim))
            self.oldpop[i].init()
            self.oldpop[i] = self.evalPop(self.oldpop[i])
            self.oldpop[i].initWal(self.model)
            self.oldpop[i].genImproveS(minimize)

        # self.bsf is a single individual recording the best-so-far solution
        #self.bsf = copy.deepcopy(self.oldpop[self.bestIndiv(self.oldpop, minimize, fitName)])
        self.bsf = individual.Individual(oldIndiv=self.oldpop[self.bestIndiv(self.oldpop, minimize, fitName)])
        #print self.bsf.fit
        self.model.WA = [] # release memory

        initC = 1
        updateC = 0
        descT = 0
        pertT = 0
        updatePertT = 0
        updateT = 0
        self.fitEval = 0
        walkLen = 10

        initT = time.time() - start

        while self.fitEval < self.MaxFit:
            for i in xrange(len(self.oldpop)):
                start = time.time()
                improveN, bestI = self.oldpop[i].steepFitDesc(minimize)
                descT = descT + time.time() - start

                if improveN == False:
                    initC = initC + 1
                    if restart == True:
                        start = time.time()
                        diff, self.oldpop[i] = self.walk(fitName, minimize, False, walkLen, self.oldpop[i])
         #               print 'walk',self.bsf.fit, self.fitEval

                        pertT = pertT + time.time() - start

                        start = time.time()
                        for j in diff:
                            #                            self.oldpop[i].fit = self.oldpop[i].fit - 2*self.oldpop[i].sumArr[j]
                            self.oldpop[i].updateEval(j)
                            self.oldpop[i].update(j)
                            self.oldpop[i].updateWAS(j)
                            self.oldpop[i].updatePertImprS(j, minimize)
                        updatePertT = updatePertT + time.time() - start

                        self.fitEval = self.fitEval + len(diff)
                    else:
                        return { 'nEvals': self.fitEval, 'sol': self.oldpop[i].fit, 'bit':self.oldpop[i].bit}
                else : # improveN is TRUE
                    start = time.time()
                    #                    self.oldpop[i].fit = self.oldpop[i].fit - 2*self.oldpop[i].sumArr[bestI]
                    self.oldpop[i].updateEval(bestI)
                    self.oldpop[i].update(bestI)
                    self.oldpop[i].updateWAS(bestI)
                    self.oldpop[i].updateImprS(bestI, minimize)
                    self.fitEval = self.fitEval + 1
          #          print 'decs', self.bsf.fit, self.fitEval
                    updateT = updateT + time.time() - start
                    updateC = updateC + 1
                    self.oldpop[i].flip(bestI)
                    ## if self.oldpop[i].bit[bestI] == '1':
                    ##     self.oldpop[i].bit[bestI] = '0'
                    ## else:
                    ##     self.oldpop[i].bit[bestI] = '1'
        for i in xrange(len(self.oldpop)):
            self.oldpop[i].destructorWal(fitName)
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit, 'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'initC':initC, 'updateC':updateC}

    def runBeamMeanSCwalkNext(self,fitName, minimize, restart, beamWidth):
        """
        BEAM next descent local search running on SC,
        where beamWidth is the number of beams
        """
        self.fitEval = 0
        start = time.time()
        self.model.transWal()
        # self.model.initInter()
        self.oldpop = []
        for i in xrange(beamWidth):
            self.oldpop.append(individual.Individual(n=self.dim, neigh=True))
            self.oldpop[i].init()
            self.oldpop[i] = self.evalPop(self.oldpop[i])
            self.oldpop[i].initWal(self.model)
            self.oldpop[i].initSC()
            self.oldpop[i].compFitG()
            self.oldpop[i].genImproveSC(minimize)

        self.bsf = individual.Individual(oldIndiv=self.oldpop[self.bestIndiv(self.oldpop, minimize, fitName)])
        self.model.WA = [] # release memory

        initC = 1
        updateC = 0

        descT = 0
        pertT = 0
        updatePertT = 0
        updateT = 0
        self.fitEval = 0
        walkLen = 10

        initT = time.time() - start

        while self.fitEval < self.MaxFit:
            for i in xrange(len(self.oldpop)):
                start = time.time()
                improveN, bestI = self.oldpop[i].nextDesc()
                descT = descT + time.time() - start

                if improveN == False:
                    initC = initC + 1
                    if restart == True:
                        start = time.time()
                        diff, self.oldpop[i] = self.walk(fitName, minimize, False, walkLen, self.oldpop[i])
                        pertT = pertT + time.time() - start

                        start = time.time()
                        for j in diff:
                            #                      self.oldpop[i].fit = self.oldpop[i].fit - 2*self.oldpop[i].sumArr[j]
                            #                     self.oldpop[i].fitG = self.oldpop[i].fit - 2/float(self.dim) * (np.sum(self.oldpop[i].sumArr))
                            self.oldpop[i].updateEval(j)
                            self.oldpop[i].compFitG()
                            self.oldpop[i].update(j)
                            self.oldpop[i].updateSC(j)
                            self.oldpop[i].updateWAS(j)
                            self.oldpop[i].updatePertImprSC(j, minimize)
                        updatePertT = updatePertT + time.time() - start
                        self.fitEval = self.fitEval + len(diff)
                    else:
                        return { 'nEvals': self.fitEval, 'sol': self.oldpop[i].fit, 'fitG': self.oldpop[i].fitG, 'bit':self.oldpop[i].bit}
                else : # improveN is TRUE
                    start = time.time()
                    self.oldpop[i].updateEval(bestI)
                    self.oldpop[i].compFitG()
                    self.oldpop[i].update(bestI)
                    self.oldpop[i].updateSC(bestI)
                    self.oldpop[i].updateWAS(bestI)
                    self.oldpop[i].updateImprSC(bestI, minimize)
                    self.fitEval = self.fitEval + 1
                    updateT = updateT + time.time() - start
                    updateC = updateC + 1
                    self.oldpop[i].flip(bestI)
                    ## if self.oldpop[i].bit[bestI] == '1':
                    ##     self.oldpop[i].bit[bestI] = '0'
                    ## else:
                    ##     self.oldpop[i].bit[bestI] = '1'
        for i in xrange(len(self.oldpop)):
            self.oldpop[i].destructorWal(fitName)
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit,'fitG': self.bsf.fitG, 'bit':self.bsf.bit, 'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'initC':initC, 'updateC':updateC}

    def runBeamMeanSCwalk(self,fitName, minimize, restart, beamWidth):
        """
        BEAM steepest descent local search running on SC,
        where beamWidth is the number of beams
        """
        self.fitEval = 0
        start = time.time()
        self.model.transWal()
        # self.model.initInter()
        self.oldpop = []
        for i in xrange(beamWidth):
            self.oldpop.append(individual.Individual(n=self.dim, neigh=True))
            self.oldpop[i].init()
            self.oldpop[i] = self.evalPop(self.oldpop[i])
            self.oldpop[i].initWal(self.model)
            self.oldpop[i].initSC()
            self.oldpop[i].compFitG()
            self.oldpop[i].genImproveSC(minimize)

        self.bsf = individual.Individual(oldIndiv=self.oldpop[self.bestIndiv(self.oldpop, minimize, fitName)])
        self.model.WA = [] # release memory

        initC = 1
        updateC = 0

        descT = 0
        pertT = 0
        updatePertT = 0
        updateT = 0
        self.fitEval = 0
        walkLen = 10

        initT = time.time() - start

        while self.fitEval < self.MaxFit:
            for i in xrange(len(self.oldpop)):
                start = time.time()
                improveN, bestI = self.oldpop[i].steepMeanDesc(minimize)
                descT = descT + time.time() - start

                if improveN == False:
                    initC = initC + 1
                    if restart == True:
                        start = time.time()
                        diff, self.oldpop[i] = self.walk(fitName, minimize, False, walkLen, self.oldpop[i])
#                        print self.bsf.fit, self.bsf.fitG, self.fitEval
                        pertT = pertT + time.time() - start

                        start = time.time()
                        for j in diff:
                            #                      self.oldpop[i].fit = self.oldpop[i].fit - 2*self.oldpop[i].sumArr[j]
                            #                     self.oldpop[i].fitG = self.oldpop[i].fit - 2/float(self.dim) * (np.sum(self.oldpop[i].sumArr))
                            self.oldpop[i].updateEval(j)
                            self.oldpop[i].compFitG()
                            self.oldpop[i].update(j)
                            self.oldpop[i].updateSC(j)
                            self.oldpop[i].updateWAS(j)
                            self.oldpop[i].updatePertImprSC(j, minimize)
                        updatePertT = updatePertT + time.time() - start
                        self.fitEval = self.fitEval + len(diff)
                    else:
                        return { 'nEvals': self.fitEval, 'sol': self.oldpop[i].fit, 'fitG': self.oldpop[i].fitG, 'bit':self.oldpop[i].bit}
                else : # improveN is TRUE
                    start = time.time()
#                    print 'bestI', bestI
#                    self.oldpop[i].fit = self.oldpop[i].fit - 2*self.oldpop[i].sumArr[bestI]
                    self.oldpop[i].updateEval(bestI)
#                    self.oldpop[i].fitG = self.oldpop[i].fit - 2/float(self.dim) * (np.sum(self.oldpop[i].sumArr))
                    self.oldpop[i].compFitG()
                    self.oldpop[i].update(bestI)
                    self.oldpop[i].updateSC(bestI)
                    self.oldpop[i].updateWAS(bestI)
                    self.oldpop[i].updateImprSC(bestI, minimize)
                    self.fitEval = self.fitEval + 1
                    updateT = updateT + time.time() - start
                    updateC = updateC + 1
                    self.oldpop[i].flip(bestI)
                    ## if self.oldpop[i].bit[bestI] == '1':
                    ##     self.oldpop[i].bit[bestI] = '0'
                    ## else:
                    ##     self.oldpop[i].bit[bestI] = '1'
        for i in xrange(len(self.oldpop)):
            self.oldpop[i].destructorWal(fitName)
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit,'fitG': self.bsf.fitG, 'bit':self.bsf.bit, 'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'initC':initC, 'updateC':updateC}

        ## self.fitEval = 0
        ## start = time.time()
        ## self.model.transWal()
        ## self.model.initInter()
        ## self.oldpop = []
        ## for i in xrange(beamWidth):
        ##     self.oldpop.append(individual.Individual(n=self.dim, neigh=True))
        ##     self.oldpop[i].init()
        ##     self.oldpop[i] = self.evalPop(self.oldpop[i])
        ##     self.oldpop[i].initWal(self.model)
        ##     self.oldpop[i].initSC()
        ##     self.oldpop[i].fitG = self.oldpop[i].fit - 2/float(self.dim) * (sum(self.oldpop[i].sumArr))
        ##     self.oldpop[i].genImproveSC(minimize)

        ## # self.bsf is a single individual recording the best-so-far solution
        ## self.bsf = copy.deepcopy(self.oldpop[self.bestIndiv(self.oldpop, minimize, fitName)])
        ## self.model.WA = [] # release memory

        ## initC = 1
        ## updateC = 0

        ## descT = 0
        ## pertT = 0
        ## updatePertT = 0
        ## updateT = 0
        ## self.fitEval = 0
        ## walkLen = 10

        ## initT = time.time() - start

        ## while self.fitEval < self.MaxFit:
        ##     for i in xrange(len(self.oldpop)):
        ##         start = time.time()
        ##         improveN, bestI = self.oldpop[i].steepMeanDesc(minimize)
        ##         descT = descT + time.time() - start

        ##         if improveN == False:
        ##             initC = initC + 1
        ##             if restart == True:
        ##                 start = time.time()
        ##                 diff, self.oldpop[i] = self.walk(fitName, minimize, False, walkLen, self.oldpop[i])

        ##                 pertT = pertT + time.time() - start

        ##                 start = time.time()
        ##                 for j in diff:
        ##                     self.oldpop[i].fit = self.oldpop[i].fit - 2*self.oldpop[i].sumArr[j]
        ##                     self.oldpop[i].fitG = self.oldpop[i].fit - 2/float(self.dim) * (np.sum(self.oldpop[i].sumArr))
        ##                     self.oldpop[i].update(j)
        ##                     self.oldpop[i].updateSC(j)
        ##                     self.oldpop[i].updateWAS(j)
        ##                     self.oldpop[i].updatePertImprSC(j, minimize)
        ##                 updatePertT = updatePertT + time.time() - start

        ##                 self.fitEval = self.fitEval + len(diff)
        ##             else:
        ##                 return { 'nEvals': self.fitEval, 'sol': self.oldpop[i].fit, 'fitG': self.oldpop[i].fitG, 'bit':self.oldpop[i].bit}
        ##         else : # improveN is TRUE
        ##             start = time.time()
        ##             self.oldpop[i].fit = self.oldpop[i].fit - 2*self.oldpop[i].sumArr[bestI]
        ##             self.oldpop[i].fitG = self.oldpop[i].fit - 2/float(self.dim) * (np.sum(self.oldpop[i].sumArr))
        ##             self.oldpop[i].update(bestI)
        ##             self.oldpop[i].updateSC(bestI)
        ##             self.oldpop[i].updateWAS(bestI)
        ##             self.oldpop[i].updateImprSC(bestI, minimize)
        ##             self.fitEval = self.fitEval + 1
        ##             updateT = updateT + time.time() - start
        ##             updateC = updateC + 1
        ##             if self.oldpop[i].bit[bestI] == '1':
        ##                 self.oldpop[i].bit[bestI] = '0'
        ##             else:
        ##                 self.oldpop[i].bit[bestI] = '1'
        ## return {'nEvals': self.fitEval, 'sol': self.bsf.fit,'fitG': self.bsf.fitG, 'bit':self.bsf.bit, 'init':initT, 'descT':descT, 'pertT':pertT, 'updateT':updateT, 'updatePertT':updatePertT, 'initC':initC, 'updateC':updateC}

    def diffBits(self, a, b):
        diff = []
        for i in range(self.dim):
            if a[i] != b[i]:
                diff.append(i)
        return diff

    def restart(self, fitName, minimize, evaluate):
        oldbit = self.oldindiv.bit
        # print 'oldfit', self.oldindiv.fit, 'bsf fit', self.bsf.fit

        if fitName == 'fit' and minimize == True :
            if self.bsf.fit > self.oldindiv.fit:
                # print 'update'
                self.bsf = individual.Individual(oldIndiv = self.oldindiv)
        elif fitName == 'fit' and minimize == False :
            if self.bsf.fit < self.oldindiv.fit:
                # print 'update'
                self.bsf = individual.Individual(oldIndiv = self.oldindiv)
        elif minimize == True :
            if self.bsf.fitG > self.oldindiv.fitG:
                # print 'update'
                self.bsf = individual.Individual(oldIndiv = self.oldindiv)
        elif minimize == False :
            if self.bsf.fitG < self.oldindiv.fitG:
                # print 'update'
                self.bsf = individual.Individual(oldIndiv = self.oldindiv)

        self.oldindiv.init()
        diff = self.diffBits(oldbit, self.oldindiv.bit)
        self.fitEval = self.fitEval + len(diff)
        # print 'diff', len(diff)
        if fitName == 'fit':
            if evaluate == True:
                self.oldindiv = self.evalPop(self.oldindiv)
        else :
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
                self.bsf = individual.Individual(oldIndiv = oldindiv)
#                self.bsf = copy.deepcopy(oldindiv)
#               self.bsf = oldindiv
        elif fitName == 'fit' and minimize == False :
            if self.bsf.fit < oldindiv.fit:
                self.bsf = individual.Individual(oldIndiv = oldindiv)
                #self.bsf = copy.deepcopy(oldindiv)
#                self.bsf = oldindiv
        elif minimize == True :
            if self.bsf.fitG > oldindiv.fitG:
                self.bsf = individual.Individual(oldIndiv = oldindiv)
#                self.bsf = copy.deepcopy(oldindiv)
 #               self.bsf = oldindiv
        elif minimize == False :
            if self.bsf.fitG < oldindiv.fitG:
                self.bsf = individual.Individual(oldIndiv = oldindiv)
#                self.bsf = copy.deepcopy(oldindiv)
                #self.bsf = oldindiv

        flipBits = random.sample(xrange(self.dim), length)
        for i in flipBits:
#            print 'before',oldindiv.bit[i]
            oldindiv.flip(i)
#            print 'after',oldindiv.bit[i]
            ## if oldindiv.bit[i] == '1':
            ##     oldindiv.bit[i] = '0'
            ## else:
            ##     oldindiv.bit[i] = '1'
        return flipBits, oldindiv

    def neighbors(self):
        """
        obsolete
        """
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
        # print 'evalPop'
        # print self.func
        # print 'len',len(indiv.bit)
        # print 'eval', self.fitEval
        indiv.fit = self.func(indiv.bit)
        # print 'evalPop'
        # self.fitEval = self.fitEval + 1
        #return copy.deepcopy(indiv)
        # print 'evalPop'
        return indiv

    def evalPopNeigh(self, indiv, fitName, minimize):
        """ evaluate the individual itself """
        cdef char* neighStr
        indiv.fit = self.func(indiv.bit)
        self.fitEval = self.fitEval + 1
        """ evaluate all neighborhood """
        fitN = np.zeros(self.dim)
        for j in range(self.dim):
            # flip the jth bit in bit-string
            # neighStr = np.copy(indiv.bit)
            neighStr = indiv.bit
            strcpy(neighStr, indiv.bit)
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
        return indiv
        # return copy.deepcopy(indiv)

    def selectionFit(self, minimize, improveA, bestFit, i):
        threshold = self.indiv.threshold

        if (minimize == True and bestFit > self.indiv.fit + threshold) or (minimize==False and bestFit < self.indiv.fit - threshold): # no equal moves
            bestFit = self.indiv.fit
            improveA = [i]
            return improveA, bestFit
        elif (abs(bestFit - self.indiv.fit)<threshold and len(improveA)!=0):
            improveA.append(i)

        return improveA, bestFit
        # if minimize == True:
        #     # if self.oldindiv.fit > self.indiv.fit:
        #     if bestFit > self.indiv.fit + self.threshold:   # an improving move
        #         bestFit = self.indiv.fit
        #         improveA = [self.indiv]
        #         # self.oldindiv = individual.Individual(oldIndiv=self.indiv)
        #         return True
        #     elif abs(bestFit - self.indiv.fit)<self.threshold:   # equal move
        #         improveA.append(self.indiv)
        #         return False
        #     else:
        #          return False
        # else: # for maximization
        #     # if self.oldindiv.fit < self.indiv.fit:
        #     if bestFit < self.indiv.fit - self.threshold:   # an improving move
        #         bestFit = self.indiv.fit
        #         # self.oldindiv = individual.Individual(oldIndiv=self.indiv)
        #         return True
        #     else:
        #         if abs(bestFit - self.indiv.fit)<self.threshold:   # equal move
        #             improveA.append(self.indiv)
        #         return False

    def selectionFitNeigh(self, minimize):
        if minimize == True :
            if self.oldindiv.fitG > self.indiv.fitG - self.oldindiv.threshold:
                self.oldindiv = copy.deepcopy(self.indiv)
                return True
            else:
                return False
        else: # maximization
            if self.oldindiv.fitG < self.indiv.fitG + self.oldindiv.threshold:
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

