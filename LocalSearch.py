""" 
Local Search module:
    This local search module gets the problem (in terms of objective function),
    and return a solution found by hill climber
"""

import numpy as np
#import matplotlib.pyplot as plt
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

    def run(self, fitName, minimize, restart, compM):
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
        
    def hyperSearchFit(self,fitName, minimize, restart):
        """ 
        performing hyper search using the probability generated by Hyperplane analysis
        """

        self.fitEval = 0
        start = os.times()[0]
        self.oldindiv = Struct( fit = 0, bit = self.genSolProp(self.model.sumFitA) )
#        self.oldindiv = self.initIndiv(self.dim)
        self.oldindiv = self.evalPop(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)
        self.initWal()
        
        self.bsf = copy.deepcopy(self.oldindiv)
        #self.model.WA = []
        
        init = False
        updateT = 0
        walkLen = 10
        initT = os.times()[0] - start
        start = os.times()[0]
        while self.fitEval < self.MaxFit:
            if init == False:
                improveN, bestI, evalCount = self.genFitBest(minimize)
                init = True
            else:
                improveN, bestI, evalCount = self.updateFitBest(bestI,minimize)
            self.fitEval = self.fitEval + evalCount
        
            if improveN == False:
                if restart == True:
                    updateT = updateT + os.times()[0] - start
                    startR = os.times()[0]
                    self.oldindiv = self.evalPop(self.oldindiv)

#                    oldbit = self.oldindiv.bit
#                    self.fitEval = self.fitEval - 1
#                    self.hyperRestart(fitName, minimize, False)
#                    diff = self.diffBits(oldbit, self.oldindiv.bit)

                    diff = self.walk(fitName, minimize,False, walkLen)
                    init = False
                    for i in diff:
                        self.update(i)
                        self.updateWAS(i)
                    initT = initT + os.times()[0] - startR
                    start = os.times()[0]
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
        updateT = updateT + os.times()[0] - start
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit,'init':initT, 'update':updateT} 

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

    def runFitsm(self,fitName, minimize, restart):
        """ 
        steepest descent local search running on S, with supermove enable
        """
        self.fitEval = 0

        start = os.times()[0]
        self.model.transWal()
        self.oldindiv = self.initIndiv(self.dim)
        self.oldindiv = self.evalPop(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)
        self.initWal()
        #self.printInter()
        
        self.bsf = copy.deepcopy(self.oldindiv)
        
        self.model.WA = []
        #self.printInter()
        
        init = False
        updateT = 0
        initT = os.times()[0] - start
        start = os.times()[0]
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
                    updateT = updateT + os.times()[0] - start
                    startR = os.times()[0]
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
                    initT = initT + os.times()[0] - startR
                    start = os.times()[0]
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
    #                updateCT = updateCT + os.times()[0] - start
    #                bestBitsCount[bestI] = bestBitsCount[bestI] + 1
                    if self.oldindiv.bit[i] == '1':
                        self.oldindiv.bit[i] = '0'
                    else:
                        self.oldindiv.bit[i] = '1'

        self.bsf = self.evalPop(self.bsf)
        updateT = updateT + os.times()[0] - start
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit,'init':initT, 'update':updateT}

    def runFitSwalk(self,fitName, minimize, restart):
        """ 
        steepest descent local search running on S
        """
        self.fitEval = 0
        start = os.times()[0]
        self.model.transWal()
        self.oldindiv = self.initIndiv(self.dim)
        self.oldindiv = self.evalPop(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)
        self.initWal()
        self.bsf = copy.deepcopy(self.oldindiv)
        self.genImproveS(minimize)
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
            improveN, bestI = self.steepFitDesc(minimize)
            descT = descT + os.times()[0] - start
            
            if improveN == False:
                initC = initC + 1
                if restart == True:
                    start = os.times()[0]
                    diff = self.walk(fitName, minimize,False, walkLen)
                    pertT = pertT + os.times()[0] - start

                    start = os.times()[0]
                    for i in diff:
                        self.oldindiv.fit = self.oldindiv.fit - 2*self.sumArr[i]
                        self.update(i)
                        self.updateWAS(i)
                        self.updatePertImprS(i, minimize)
                    updatePertT = updatePertT + os.times()[0] - start

                    self.fitEval = self.fitEval + len(diff)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                start = os.times()[0]

                self.oldindiv.fit = self.oldindiv.fit - 2*self.sumArr[bestI] 

                self.update(bestI)
                self.updateWAS(bestI)
                self.updateImprS(bestI, minimize)
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
        self.oldindiv = self.initIndiv(self.dim)
        self.oldindiv = self.evalPop(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)
        self.initWal()
        self.bsf = copy.deepcopy(self.oldindiv)
        self.genImproveS(minimize)
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
            improveN, bestI = self.nextDesc()
            descT = descT + os.times()[0] - start
            
            if improveN == False:
                initC = initC + 1
                if restart == True:

                    start = os.times()[0]
                    diff = self.walk(fitName, minimize,False, walkLen)
                    pertT = pertT + os.times()[0] - start

                    start = os.times()[0]
                    for i in diff:
                        self.oldindiv.fit = self.oldindiv.fit - 2*self.sumArr[i]
                        self.update(i)
                        self.updateWAS(i)
                        self.updatePertImprS(i, minimize)
                    updatePertT = updatePertT + os.times()[0] - start

                    self.fitEval = self.fitEval + len(diff)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                start = os.times()[0]
                self.oldindiv.fit = self.oldindiv.fit - 2*self.sumArr[bestI]
                self.update(bestI)
                self.updateWAS(bestI)
                self.updateImprS(bestI, minimize)
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
        self.indiv = copy.deepcopy(self.oldindiv)
        self.initWal()
        self.bsf = copy.deepcopy(self.oldindiv)
        self.model.WA = []

        walkLen = 10
        init = False
        updateT = 0
        initT = os.times()[0] - start
        start = os.times()[0]

    def runFitS2walk(self,fitName, minimize, restart):
        """ 
        steepest descent local search running on S
        """
        self.fitEval = 0
        start = os.times()[0]
        self.model.transWal()
        self.oldindiv = self.initIndiv(self.dim)
        self.oldindiv = self.evalPop(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)
        self.initWal()
        self.bsf = copy.deepcopy(self.oldindiv)
        self.model.WA = []

        walkLen = 10
        init = False
        updateT = 0
        initT = os.times()[0] - start
        start = os.times()[0]

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
            print self.oldindiv.bit
            print improveN, bestI, self.fitEval
            print 
#            pdb.set_trace()
        
            if improveN == False:
                if restart == True:
                    updateT = updateT + os.times()[0] - start
                    startR = os.times()[0]
                    self.oldindiv = self.evalPop(self.oldindiv)
                    diff = self.walk(fitName, minimize,False, walkLen)
                    init = False

#                    print self.bsf.fit
#                    pdb.set_trace()

                    for i in diff:
                        self.update(i)
                        self.updateWAS(i)
                    initT = initT + os.times()[0] - startR
                    start = os.times()[0]
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
        updateT = updateT + os.times()[0] - start
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'bit':self.bsf.bit,'init':initT, 'update':updateT}

    def runFitWal(self,fitName, minimize, restart):
        """ 
        Walsh Analysis for speeding up steepest descent local search
        """
        self.oldindiv = self.initIndiv(self.dim)
        self.fitEval = 0
        
        self.model.transWal()
        self.indiv = copy.deepcopy(self.oldindiv)
        self.initWal()
#        self.initInter()

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
        start = os.times()[0]
        self.model.transWal()
        self.oldindiv = self.initIndiv(self.dim)
        self.oldindiv = self.evalPop(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)
        self.initWal()
        self.bsf = copy.deepcopy(self.oldindiv)
        self.genImproveS(minimize)
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
            improveN, bestI = self.steepFitDesc(minimize)
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
                        self.oldindiv.fit = self.oldindiv.fit - 2*self.sumArr[i]
                        self.update(i)
                        self.updateWAS(i)
                        self.updatePertImprS(i, minimize)
                    updatePertT = updatePertT + os.times()[0] - start

                    self.fitEval = self.fitEval + len(diff)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                start = os.times()[0]

                self.oldindiv.fit = self.oldindiv.fit - 2*self.sumArr[bestI] 

                self.update(bestI)
                self.updateWAS(bestI)
                self.updateImprS(bestI, minimize)
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
        self.oldindiv = self.initIndiv(self.dim)
        self.oldindiv = self.evalPop(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)
        self.initWal()
        self.bsf = copy.deepcopy(self.oldindiv)
        self.genImproveS(minimize)
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
            improveN, bestI = self.nextDesc()
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
                        self.oldindiv.fit = self.oldindiv.fit - 2*self.sumArr[i]
                        self.update(i)
                        self.updateWAS(i)
                        self.updatePertImprS(i, minimize)
                    updatePertT = updatePertT + os.times()[0] - start

                    self.fitEval = self.fitEval + len(diff)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                start = os.times()[0]

                self.oldindiv.fit = self.oldindiv.fit - 2*self.sumArr[bestI] 

                self.update(bestI)
                self.updateWAS(bestI)
                self.updateImprS(bestI, minimize)
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
        self.oldindiv = self.initIndivNeigh(self.dim)
        self.oldindiv = self.evalPop(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)
        self.initWal()
        self.initSC()
        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (sum(self.sumArr))
        self.indiv.fitG = self.oldindiv.fitG
        self.genImproveSC(minimize)
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
            improveN, bestI = self.steepMeanDesc(minimize)
            descT = descT + os.times()[0] - start

            if improveN == False:
                initC = initC + 1
                if restart == True:

                    start = os.times()[0]
                    diff = self.walk(fitName, minimize,False, walkLen)
                    pertT = pertT + os.times()[0] - start
                    
                    start = os.times()[0]
                    for i in diff:

                        self.oldindiv.fit = self.oldindiv.fit - 2*self.sumArr[i] 
                        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.sumArr))

                        self.update(i)
                        self.updateSC(i)
                        self.updateWAS(i)
                        self.updatePertImprSC(i, minimize)
                    updatePertT = updatePertT + os.times()[0] - start

                    self.fitEval = self.fitEval + len(diff)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                start = os.times()[0]

                self.oldindiv.fit = self.oldindiv.fit - 2*self.sumArr[bestI] 
                self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.sumArr))

                self.update(bestI)
                self.updateSC(bestI)
                self.updateWAS(bestI)
                self.updateImprSC(bestI, minimize)
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
        self.oldindiv = self.initIndivNeigh(self.dim)
        self.oldindiv = self.evalPop(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)
        self.initWal()
        self.initSC()
        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (sum(self.sumArr))
        self.indiv.fitG = self.oldindiv.fitG
        self.genImproveSC(minimize)
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
            improveN, bestI = self.nextDesc()
            descT = descT + os.times()[0] - start

            if improveN == False:
                initC = initC + 1
                if restart == True:

                    start = os.times()[0]
                    diff = self.walk(fitName, minimize,False, walkLen)
                    pertT = pertT + os.times()[0] - start
                    
                    start = os.times()[0]
                    for i in diff:
                        self.oldindiv.fit = self.oldindiv.fit - 2*self.sumArr[i] 
                        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.sumArr))

                        self.update(i)
                        self.updateSC(i)
                        self.updateWAS(i)
                        self.updatePertImprSC(i, minimize)
                    updatePertT = updatePertT + os.times()[0] - start

                    self.fitEval = self.fitEval + len(diff)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                start = os.times()[0]
                self.oldindiv.fit = self.oldindiv.fit - 2*self.sumArr[bestI] 
                self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.sumArr))
                self.update(bestI)
                self.updateSC(bestI)
                self.updateWAS(bestI)
                self.updateImprSC(bestI, minimize)
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
        self.oldindiv = self.initIndivNeigh(self.dim)
        self.oldindiv = self.evalPop(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)
        self.initWal()
        self.initSC()
        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (sum(self.sumArr))
        self.indiv.fitG = self.oldindiv.fitG
        self.genImproveSC(minimize)
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
            improveN, bestI = self.steepMeanDesc(minimize)
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
                        self.oldindiv.fit = self.oldindiv.fit - 2*self.sumArr[i] 
                        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.sumArr))

                        self.update(i)
                        self.updateSC(i)
                        self.updateWAS(i)
                        self.updatePertImprS(i, minimize)
                    updatePertT = updatePertT + os.times()[0] - start

                    self.fitEval = self.fitEval + len(diff)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                start = os.times()[0]

                self.oldindiv.fit = self.oldindiv.fit - 2*self.sumArr[bestI] 
                self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.sumArr))

                self.update(bestI)
                self.updateSC(bestI)
                self.updateWAS(bestI)
                self.updateImprSC(bestI, minimize)
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
        self.oldindiv = self.initIndivNeigh(self.dim)
        self.oldindiv = self.evalPop(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)
        self.initWal()
        self.initSC()
        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (sum(self.sumArr))
        self.indiv.fitG = self.oldindiv.fitG
        self.genImproveSC(minimize)
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
            improveN, bestI = self.nextDesc()
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
                        self.oldindiv.fit = self.oldindiv.fit - 2*self.sumArr[i] 
                        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.sumArr))

                        self.update(i)
                        self.updateSC(i)
                        self.updateWAS(i)
                        self.updatePertImprS(i, minimize)
                    updatePertT = updatePertT + os.times()[0] - start

                    self.fitEval = self.fitEval + len(diff)
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                start = os.times()[0]

                self.oldindiv.fit = self.oldindiv.fit - 2*self.sumArr[bestI] 
                self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.sumArr))

                self.update(bestI)
                self.updateSC(bestI)
                self.updateWAS(bestI)
                self.updateImprSC(bestI, minimize)
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
        self.indiv = copy.deepcopy(self.oldindiv)
        self.initWal()
        
        self.initSC()
        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (sum(self.sumArr))
        self.indiv.fitG = self.oldindiv.fitG

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
                improveN, bestI, evalCount = self.updateMeanBest(bestI,minimize)
            self.fitEval = self.fitEval + self.dim

            if improveN == False:
                if restart == True:
                    updateT = updateT + os.times()[0] - start
                    startR = os.times()[0]
                    self.oldindiv = self.evalPop(self.oldindiv)
                    self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (np.sum(self.sumArr))
                    
                    diff = self.walk( fitName, minimize, False, walkLen )
                    init = False
                    
                    for i in diff:
                        self.update(i)

                        self.updateSC(i)
                        self.updateWAS(i)
                    initT = initT + os.times()[0] - startR
                    start = os.times()[0]
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
        updateT = updateT + os.times()[0] - start
        return {'nEvals': self.fitEval, 'sol': self.bsf.fit, 'fitG': self.bsf.fitG, 'bit':self.bsf.bit,'init':initT, 'update':updateT}

    def runMeanWal(self,fitName, minimize, restart):
        """ 
        steepest descent local search with respect to mean of neighs by Walsh Analysis
        """
        self.oldindiv = self.initIndivNeigh(self.dim)
        self.fitEval = 0
        
        start = os.times()[0]
        self.model.transWal()
        print 'transWal', os.times()[0] - start

        self.indiv = copy.deepcopy(self.oldindiv)

        start = os.times()[0]
        self.initWal()
        print 'initWal', os.times()[0] - start

        self.oldindiv = self.evalPop(self.oldindiv)
        self.oldindiv.fitG = self.oldindiv.fit - 2/float(self.dim) * (sum(self.sumArr))
        self.bsf = copy.deepcopy(self.oldindiv)
        self.indiv = copy.deepcopy(self.oldindiv)

        self.model.WA = []
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

                start = os.times()[0]
                self.indiv.fitG = oldFitG - 2*self.sumArr[nCount] + 4/float(self.dim) * self.compPhisum(nCount)
                compPSumT = compPSumT + os.times()[0] - start
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
                    start = os.times()[0]
                    for i in diff:
                        self.update(i)
                    updateT = updateT + os.times()[0] - start
                else:
                    return { 'nEvals': self.fitEval, 'sol': self.oldindiv.fit, 'fitG': self.oldindiv.fitG, 'bit':self.oldindiv.bit}
            else : # improveN is TRUE 
                start = os.times()[0]
                self.update(changeBit)
                updateT = updateT + os.times()[0] - start

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
            self.oldindiv = self.initIndiv(self.dim)
            if evaluate == True:
                self.oldindiv = self.evalPop(self.oldindiv)
        else :
            self.oldindiv = self.initIndivNeigh(self.dim)
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

    def genImproveS(self,minimize):
        """
        generate the index of best neigh according to sumArr only (surrogate of fitness)
        """
        # check improving move
        self.improveA = []
        for i in range(self.dim):
            if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i]<-self.threshold):
                self.improveA.append(i) 

    def genFitNext(self,minimize):
        """
        generate the index of next improving neigh according to sumArr only (surrogate of fitness)
        """
        # check improving move
        improve = False
        self.improveA = []
        for i in range(self.dim):
            if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i]<-self.threshold):
                self.improveA.append(i) 
                improve = True

        if improve == False:
            return False, None

        # randomly pick an improving move, which takes only constant time 
        bestI = random.choice(self.improveA)
                    
        return True, bestI

    def genFitBest2(self,minimize):
        """
        generate the index of best distance 2 neighborhoods according to sumArr only (surrogate of fitness)

        return: 1) whether there is an improving move in distance 2 neigh
                2) the index of best distance 1 neigh for taking the next move
                3) the number of evaluations consumed by this step
        """
        improve = False
        self.improveA = []
        neighImprove = []

        # checking the distance 1 neigh
        for i in range(self.dim):
            if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i]<-self.threshold):
                self.improveA.append(i) 
                neighImprove.append(Struct(index =[i], val = self.sumArr[i])) # add dist 1 neigh into consideration as well
                improve = True

#        print self.sumArr
#        print self.improveA

        # checking the distance 2 neigh, remember to preserve context
        for i in range(self.dim):
            self.mirrorParam() # everything is reset, pretending nothing happened
            self.updateFake(i)
            #self.updateWASfake(i)
            for j in [k for k in range(self.dim) if k!=i]:
                self.sumArrFake[j] = self.sumArrFake[j]+self.sumArr[i]
                if (minimize == True and self.sumArrFake[j] > self.threshold) or (minimize == False and self.sumArrFake[j]<-self.threshold):
                    neighImprove.append(Struct(index =[i,j], val = self.sumArrFake[j]))
                    improve = True

#        for i in range(len(neighImprove)):
#            print neighImprove[i].index, neighImprove[i].val

#        for i in neighImprove:
#            print i.index, i.val
               
        if improve == False:
            #return False, None, self.dim*self.dim
            return False, None, self.dim

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
        return True, bestI, self.dim

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

        if p in self.Inter:
            for i in self.Inter[p].arr: 
                if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i]< self.threshold ):
                    if i not in self.improveA:
                        self.improveA.append(i)
                elif i in self.improveA:
                    self.improveA.remove(i)

        if not self.improveA:
            return False, None

        random.shuffle(self.improveA)

        for i in self.improveA:
            if i == self.improveA[0]:
                best = self.sumArr[i]
                bestI = i
            elif (best<self.sumArr[i] and minimize == True) or (best>self.sumArr[i] and minimize == False): 
                best = self.sumArr[i]
                bestI = i
                    
        return True, bestI

    def steepFitDesc(self, minimize):
        if not self.improveA:
            return False, None

        random.shuffle(self.improveA)

        for i in self.improveA:
            if i == self.improveA[0]:
                best = self.sumArr[i]
                bestI = i
            elif (best<self.sumArr[i] - self.threshold and minimize == True) or (best>self.sumArr[i] + self.threshold and minimize == False): 
                best = self.sumArr[i]
                bestI = i
                    
        return True, bestI

    def steepMeanDesc(self, minimize):
        if not self.improveA:
            return False, None

        random.shuffle(self.improveA)

        # find the best value
        for i in self.improveA:
            if i == self.improveA[0]:
                best = self.SC[i]
                bestI = i
            elif ( best<self.SC[i] - self.threshold and minimize == True ) or ( best>self.SC[i] + self.threshold and minimize == False ): # seek for max S
                best = self.SC[i]
                bestI = i
        return True, bestI

    def updateFitNext(self, p, minimize):
        """
        find the next improving move by the similar update trick
        """
        self.improveA.remove(p)
        if p in self.Inter:
            for i in self.Inter[p].arr: 
                if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i]< self.threshold ):
                    if i not in self.improveA:
                        self.improveA.append(i)
                elif i in self.improveA:
                    self.improveA.remove(i)

        if not self.improveA:
            return False, None, evalCount

        # randomly pick an improving move, which takes only constant time 
        bestI = random.choice(self.improveA)
                   
        return True, bestI, evalCount

    def nextDesc(self):
        """
        find the next improving move by the similar update trick
        """
        if not self.improveA:
            return False, None

        # randomly pick an improving move, which takes only constant time 
        bestI = random.choice(self.improveA)
                   
        return True, bestI

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
            for j in [k for k in range(self.dim) if k!=i]:
                self.sumArrFake[j] = self.sumArrFake[j]+self.sumArr[i]
                if (minimize == True and self.sumArrFake[j] > self.threshold) or (minimize == False and self.sumArrFake[j]<-self.threshold):
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

    def genImproveSC(self,minimize):
        """
        generate the index of best neigh according to {S_p(X)-2/N \Sigma_{i=1}^{N}C_{ip}(X)} only (surrogate of fitness)
        """
        # check improving move 
        self.improveA = []
        for i in range(self.dim):
            if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i]<self.threshold):
                self.improveA.append(i)

    def genMeanBest(self,minimize):
        """
        generate the index of best neigh according to {S_p(X)-2/N \Sigma_{i=1}^{N}C_{ip}(X)} only (surrogate of fitness)
        """
        # check improving move 
        improve = False
        self.improveA = []
        for i in range(self.dim):
            if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i]<self.threshold):
                self.improveA.append(i)
                improve = True

        if improve == False:
            return False, None

        random.shuffle(self.improveA)

        # find the best value
        for i in self.improveA:
            if i == self.improveA[0]:
                best = self.SC[i]
                bestI = i
            elif ( best<self.SC[i] - self.threshold and minimize == True ) or ( best>self.SC[i] + self.threshold and minimize == False ): # seek for max S
                best = self.SC[i]
                bestI = i
        return True, bestI

    def genMeanNext(self,minimize):
        """
        generate the index of next improving neigh according to {S_p(X)-2/N \Sigma_{i=1}^{N}C_{ip}(X)} only (surrogate of fitness)
        """
        # check improving move 
        improve = False
        self.improveA = []
        for i in range(self.dim):
            if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i]<self.threshold):
                self.improveA.append(i)
                improve = True

        if improve == False:
            return False, None

        bestI = random.choice(self.improveA)
        return True, bestI

    def updateMeanBest(self, p, minimize):
        self.improveA.remove(p)
        if p in self.Inter:
            for i in self.Inter[p].arr:
                if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i]<self.threshold):
                    if i not in self.improveA:
                        self.improveA.append(i)
                elif i in self.improveA:
                    self.improveA.remove(i)

        if not self.improveA:
            return False, None

        random.shuffle(self.improveA)

        for i in self.improveA:
            if i == self.improveA[0]:
                best = self.SC[i]
                bestI = i
            elif (best<self.SC[i] - self.threshold and minimize == True) or (best>self.SC[i]+self.threshold and minimize == False): # seek for max S
                best = self.SC[i]
                bestI = i
                    
        return True, bestI

    def updateMeanNext(self, p, minimize):
        self.improveA.remove(p)
        if p in self.Inter:
            for i in self.Inter[p].arr:
                if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i]<self.threshold):
                    if i not in self.improveA:
                        self.improveA.append(i)
                elif i in self.improveA:
                    self.improveA.remove(i)

        if not self.improveA:
            return False, None

        bestI = random.choice(self.improveA)
        return True, bestI

    def initWal(self):
        """ 
        1. 
        compute the sum array for the first time, according to the initial solution
        2. 
        initialize a full matrix, Coincidence, C_ij, which is a symmetric one, and diagonal is empty (as we don't want to recompute SumArr), empty whenever i >= j.
        Two options for implementing C matrix
            a. precompute C matrix 
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

#        self.InterCount = np.zeros(self.dim)

#        for i in range(self.dim):
#            self.Inter.append(Set()) # coincidence matrix

        #self.C = dict() # coincidence matrix
        for i in range(len(self.model.WA)):
            W = int(math.pow(-1, self.binCount(self.model.WA[i].arr, self.indiv.bit))) * self.model.WA[i].w
            self.WAS[i] = Struct(arr = self.model.WA[i].arr, w = W)
            comb = self.genComb(len(self.model.WA[i].arr))
            #print i, self.model.WA[i].arr, comb

            for j in self.model.WA[i].arr:
                self.sumArr[j] = self.sumArr[j] + W
                if len(self.model.WA[i].arr)>1: # for at least order Walsh terms
                    if j not in self.Inter: # the entry of i doesn't exist yet
                        self.Inter[j] = Struct(arr=Set(), WI=Set())

                    for k in self.model.WA[i].arr:
                        if k != j:
                            self.Inter[j].arr.add(k)
                    self.Inter[j].WI.add(i)

                # add list of order >= 3 Walsh terms for the purpose of updating C matrix
                if len(self.model.WA[i].arr) >= 3:
                    if j not in self.infectBit: 
                        self.infectBit[j] = [Struct(arr=self.model.WA[i].arr, WI=i)]
                    else :
                        self.infectBit[j].append(Struct(arr=self.model.WA[i].arr, WI=i))

            for j in comb: # for each list in comb
                j0 = self.model.WA[i].arr[int(j[0])]
                j1 = self.model.WA[i].arr[int(j[1])]
#                if (j0, j1) in self.C.keys():
#                    self.C[j0,j1] = self.C[j0,j1] + W
#                elif W != 0:
#                    self.C[j0,j1] = W
                self.C[j0,j1] = self.C[j0,j1] + W
        
#        for i in range(self.dim):
#            if i in self.Inter:
#                self.InterCount[i] = len(self.Inter[i].arr)
#            else :
#                self.InterCount[i] = 0

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
                arr = i.arr[:]
                arr.remove(p)
                comb = self.genComb(len(arr))
                for k in range(len(comb)):
                    k0 = arr[int(comb[k][0])]
                    k1 = arr[int(comb[k][1])]
                    self.C[k0,k1] = self.C[k0,k1] - 2 * self.WAS[i.WI].w

    def updateImprS(self, p, minimize):
        self.improveA.remove(p)
        if p in self.Inter:
            for i in self.Inter[p].arr: 
                if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i]< self.threshold ):
                    if i not in self.improveA:
                        self.improveA.append(i)
                elif i in self.improveA:
                    self.improveA.remove(i)

    def updatePertImprS(self, p, minimize):
        if p in self.Inter:
            for i in self.Inter[p].arr : 
                if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i]< self.threshold ):
                    if i not in self.improveA:
                        self.improveA.append(i)
                elif i in self.improveA:
                    self.improveA.remove(i)
        
        if (minimize == True and self.sumArr[p] > self.threshold) or (minimize == False and self.sumArr[p]< self.threshold ):
            if p not in self.improveA:
                self.improveA.append(p)
        elif p in self.improveA:
            self.improveA.remove(p)

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
                    self.sumArrFake[i] = self.sumArrFake[i] - 2*self.Cfake[i,p]
#                    self.Cfake[i,p] = - self.Cfake[i,p]
                else:
                    self.sumArrFake[i] = self.sumArrFake[i] - 2*self.Cfake[p,i]
#                    self.Cfake[p,i] = - self.Cfake[p,i]

#        # update the rest of elements in C matrix
#        if p in self.infectBit.keys():
#            for i in self.infectBit[p]:
#                arr = copy.deepcopy(i.arr)
#                arr.remove(p)
#                comb = self.genComb(len(arr))
#                for k in range(len(comb)):
#                    k0 = arr[int(comb[k][0])]
#                    k1 = arr[int(comb[k][1])]
#                    self.Cfake[k0,k1] = self.Cfake[k0,k1] - 2 * self.WASfake[i.WI].w

    def updateWAS(self,p):
        if p in self.Inter:
            for i in self.Inter[p].WI:
                self.WAS[i].w = - self.WAS[i].w

#    def updateWAS(self,p):
#        """ fake version """
#        if p in self.Inter:
#            for i in self.Inter[p].WI:
#                self.WASfake[i].w = - self.WASfake[i].w

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
#                phi = np.zeros(self.model.k+1)
#                if p in self.Inter:
#                    for i in self.Inter[p].WI:
#                        order = len(self.WAS[i].arr)
#                        phi[order-1] = phi[order-1] + self.WAS[i].w
#
#                Z[i] = self.sumArr[p]
#                for i in range(1, self.model.k+1):
#                    if phi[i] != 0:
#                        Z[i] = Z[i] + i * phi[i]
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

    def updateImprSC(self, p, minimize):
        self.improveA.remove(p)
        if p in self.Inter:
            for i in self.Inter[p].arr:
                if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i]<self.threshold):
                    if i not in self.improveA:
                        self.improveA.append(i)
                elif i in self.improveA:
                    self.improveA.remove(i)

    def updatePertImprSC(self, p, minimize):
        if p in self.Inter:
            for i in self.Inter[p].arr:
                if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i]<self.threshold):
                    if i not in self.improveA:
                        self.improveA.append(i)
                elif i in self.improveA:
                    self.improveA.remove(i)
        if (minimize == True and self.SC[p] > self.threshold) or (minimize == False and self.SC[p]<self.threshold):
            if p not in self.improveA:
                self.improveA.append(p)
        elif p in self.improveA:
            self.improveA.remove(p)

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
        self.sumArrFake = copy.deepcopy(self.sumArr)
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
