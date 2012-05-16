# generate NK-landscapes instances
## the encoding is in right to left fashion
## 00-0, 01-1, 10-2, 11-3 import WalshAnalysis as wal import random import numpy as np import math

import random
import WalshAnalysis as wal
import numpy as np
import tool as tl
import math
from sets import Set
import pdb


class NKLandscape:
    """ NK-landscape class """
    def __init__(self,inN,inK, fileName = None):
        self.n = inN
        self.k = inK

        # for run experiments
        if fileName == None:
            self.genNeigh()
            self.genFunc()
        else:
            self.readFile(fileName)

#        # for generating benchmarks
#        self.genNeigh()
#        self.genFunc()
#        self.exportToFile(fileName)
#
        self.Kbits = tl.genSeqBits(self.k+1)

    def exportToFile(self, fileName):
        f = open(fileName, 'w')
        for i in range(self.n): 
            for j in range(len(self.neighs[i])):
                print >>f, self.neighs[i][j], '\t',
            print >>f
        for i in range(self.n): 
            for j in range(len(self.func[i])):
                print >>f, self.func[i][j], '\t',
            print >>f

    def readFile(self, fName):
        self.neighs = np.genfromtxt(fName, delimiter="\t", dtype='int', skip_footer=self.n, autostrip=True, usecols = range(self.k)).tolist()
        self.func = np.genfromtxt(fName, delimiter="\t", skip_header=self.n, autostrip=True, usecols = range(int(math.pow(2,self.k+1)))).tolist()
        
    """ generate neighborhood """
    def genNeigh(self):
        self.neighs = []
        for i in range(self.n):
            oneNeigh = random.sample(range(self.n), self.k)
            while i in oneNeigh:
                oneNeigh = random.sample(range(self.n), self.k)
            self.neighs.append(oneNeigh)
    def getNeigh(self):
        return self.neighs

    """ generate function value """
    def genFunc(self):
        self.func = []
        for i in range(self.n):
            oneFunc = []
            for j in range(int(math.pow(2,self.k+1))):
                oneFunc.append(random.random())
            self.func.append(oneFunc)
    def getFunc(self):
        return self.func
    def getN(self):
        return self.n
    def genK(self):
        return self.k

    """ compute the fitness value"""
    def compFit(self, bitStr): 
        sum = 0
        for i in range(self.n):
            """ compose interacting bits """
            if len(self.neighs) > 0:
                interBit = self.neighs[i][:]
            else:
                interBit = []
            interBit.append(i)
            interBit.sort()
            """ extract corresponding bits """
            bits = [ bitStr[int(j)] for j in interBit ]
            interStr = ''.join(bits)
            """ sum up the sub-function values """ 
            sum = sum + self.func[i][int(interStr,2)]
        return sum/float(self.n)

    def WalCof(self):
        """ compute the Walsh coefficients """
        subW = [] # subW is a N*2^K matrix
        for i in range(self.n):
            """ 1. Compute coefficients for each sub-functions """
            subWone = wal.computeW(self.Kbits, self.func[i])
            subW.append(subWone)
        w = np.zeros(math.pow(2,self.n))
        for i in range(int(math.pow(2,self.n))): # for every candidate solution
            iStr = bin(i)
            iStr = iStr[2:]
            if len(iStr) < self.n:
                iStr = (self.n - len(iStr))*'0' + iStr
            for j in range(self.n): # for every sub-function
                maskJ = self.neighs[j][:]
                maskJ.append(j)
                maskJ.sort()
                # pack iStr to a (k+1) length one
                maskStr = [iStr[k] for k in maskJ]
                maskStr = ''.join(maskStr)
                occurOneBit = self.indexOneBit(iStr)
                if self.checkInclude(occurOneBit, maskJ) == True :
                    extractBit = maskStr
                    w[i] = w[i] + subW[j][int(extractBit, 2)]
        return w

    def WalCofLinear(self):
        """ compute the Walsh coefficients in a linear time """
        subW = [] # subW is a N*2^K matrix
        for i in range(self.n):
            """ Compute coefficients for each sub-functions """
            subWone = wal.computeW(self.Kbits, self.func[i])
            subW.append(subWone)
        print 'len', math.pow(2,self.n)
        w = np.zeros(math.pow(2,self.n))
        for i in range(self.n): # i: index of subfunction
            interBits = self.neighs[i][:]
            interBits.append(i)
            interBits.sort()
            for j in range(int(math.pow(2, self.k+1))): # j: index of substrings
                indexW = self.composeFullStr(i, j, interBits, self.n)
                w[indexW] = w[indexW] + subW[i][j]

        return w/float(self.n)

    def WalshCofLinearLinklist(self):
        """ compute the Walsh Coefficients in a liner time with linear space """
        subW = [] # subW is a N * 2^K matrix
        for i in range(self.n):
            """ Compute coefficients for each sub-functions """
            subWone = wal.computeW(self.Kbits, self.func[i])
            subW.append(subWone)
        """ use dict to represent all non-zero Walsh Coefficients"""
        w = dict()
        for i in range(self.n): # i: index of sub-function
            if len(self.neighs)!=0:
                interBits = self.neighs[i][:]
            else:
                interBits = []
            interBits.append(i)
            interBits.sort()
            for j in range(int(math.pow(2, self.k+1))): # j: index of substrings
                indexW = self.composeFullBitStr(i, j, interBits, self.n)
                if w.has_key(indexW):
                    w[indexW] = w[indexW] + subW[i][j]
                else:
                    w[indexW] = subW[i][j]
        for k in w.keys():
            w[k] = w[k]/float(self.n)
        self.w = w
        return w

    """ 
    for Walsh Local Search
    """
    def transWal(self):
        """
        translate bitstring represented Walsh terms into arrays of bits that they touches
        """
        self.WA = [] # array representing Walsh terms
        for k in self.w.keys(): 
            if self.w[k] != 0:
                a = []
                for i in range(self.n):
                    if k[i] == '1':
                        a.append(i)
                self.WA.append( Struct(arr = a, w = self.w[k]) )

    def genHyperVote(self):
        """
        using the voting strategy where only best hyperplane have the chance to vote
        """
        self.transWal()
#        bit,fit = tl.compFit(self)
#        a = sorted(zip(bit,fit), key=lambda a_entry: a_entry[1]) 
#        optBit = a[0][0]
#        optFit = a[0][1]
#        print 'opti\n',optBit, optFit

        #for i in range(len(a)): 
#        for i in range(10): 
#            print '%s\t%.3f' %(a[i][0],a[i][1])

        # initialize sumFitA 
        self.sumFitA = []
        evalSubFunc = []
        for i in range(self.n):
            self.sumFitA.append(Struct(one=0,zero=0))
        
        for i in range(self.n):
            subBit = self.neighs[i][:]
            subBit.append(i)
            subBit.sort()

            if subBit not in evalSubFunc:
                evalSubFunc.append(subBit)

                # check every template that matches the subfunction
                seqBits = tl.genSeqBits(len(subBit))
                schFitArr = []
                walTouch = []

                # compute schema fitness
                for k in self.WA:
                    subset = True
                    for l in k.arr:
                        if l not in subBit:
                            subset = False
                            break
                    if subset == True:
                        walTouch.append(k)

                for j in seqBits:
                    schFit = 0

                    # convert bit string to array representation
                    schTpl = []
                    for k in range(len(j)):
                        if j[k] == '1':
                            schTpl.append(subBit[k])

                        for k in walTouch:
                            schFit = schFit + int(math.pow(-1, self.binCountArr(k.arr, schTpl))) * k.w

                    schFitArr.append(Struct(fit=schFit,arr=schTpl))
#                    print subBit, j, schFit
#                print 

                schFitArrSort = sorted(schFitArr, key = lambda i: i.fit)

                # perform voting from the best hyperplane associated with the subfunction
                #for k in range(self.k+1):
                for k in range(1):
                #for k in range(self.k*2):
                    for j in subBit:
                        if j in schFitArrSort[k].arr:
                            #self.sumFitA[j].one = self.sumFitA[j].one + schFitArrSort[k].fit
                            self.sumFitA[j].one = self.sumFitA[j].one + 1
                        else:
                            #self.sumFitA[j].zero = self.sumFitA[j].zero + schFitArrSort[k].fit
                            self.sumFitA[j].zero = self.sumFitA[j].zero + 1


#        for i in range(self.n):
#            print '%d\tOne: %.2f\tZero: %.2f' %(i, self.sumFitA[i].one, self.sumFitA[i].zero)

#            hamDist = 0
#            # compute the hamming distance
#            for i in range(self.n):
#                if sol[i] != optBit[i]:
#                    hamDist = hamDist + 1
#            print 'Hyper solution\t', sol, self.func(sol), hamDist
#
#        randSol = self.initIndiv(self.n)
#        hamDistRand = 0
#        for i in range(self.n):
#            if randSol.bit[i] != optBit[i]:
#                hamDistRand = hamDistRand + 1
#        print 'Random Solution\t', self.func(randSol.bit), hamDistRand
#        return {'nEvals': 0, 'sol': self.func(sol), 'bit': hamDist, 'init': self.func(randSol.bit), 'update': hamDistRand}

    def genHyperSqVote(self):
        """
        using the voting strategy where only best hyperplane have the chance to vote
        compose the template on the bases of union of two subfunction, in this way each variable can have more than one vote
        """
        self.transWal()
#        print 'genHyperSqVote'
#
#        bit,fit = tl.compFit(self)
#        a = sorted(zip(bit,fit), key=lambda a_entry: a_entry[1]) 
#        optBit = a[0][0]
#        optFit = a[0][1]
#        print 'opti\n',optBit, optFit
#
#        for i in range(len(a)): 
##        for i in range(10): 
#            print '%s\t%.3f' %(a[i][0],a[i][1])
        # initialize sumFitA 
    
        self.sumFitA = []
        for i in range(self.n):
            self.sumFitA.append(Struct(one=0,zero=0))

#        scan = 0
#        reuse = 0
        
        evalOuterFunc = []
        mergeFunc = []
        for i in range(self.n):
            subBitOut = self.neighs[i][:]
            subBitOut.append(i)
            subBitOut.sort()
#            print 'subBitOut', subBitOut

            if subBitOut not in evalOuterFunc:
                evalOuterFunc.append(subBitOut)

                evalInnerFunc = []

                for ii in range(i+1,self.n):
                    subBitIn = self.neighs[ii][:]
                    subBitIn.append(ii)
                    subBitIn.sort()
#                    print '\tsubBitIn', subBitIn

                    if subBitIn != subBitOut and subBitIn not in evalInnerFunc:
                        evalInnerFunc.append(subBitIn)
                        subBitIn = tl.listMerge(subBitOut,subBitIn)
                        subBitIn.sort()
                        
                        if subBitIn not in mergeFunc:
                            mergeFunc.append(subBitIn)
#                            print '\t\tsubMerge', subBitIn
                            # check every template that matches the subfunction
                            seqBits = tl.genSeqBits(len(subBitIn))
                            schFitArr = []
                            walTouch = []
                            init = False

                            for j in seqBits:
                                schFit = 0

                                # convert bit string to array representation
                                schTpl = []
                                for k in range(len(j)):
                                    if j[k] == '1':
                                        schTpl.append(subBitIn[k])

                                if init == False: 
                                    # compute schema fitness from scan over all wal cof
                                    for k in self.WA:
                                        subset = True
                                        for l in k.arr:
                                            if l not in subBitIn:
                                                subset = False
                                                break
                                        if subset == True:
                                            schFit = schFit + int(math.pow(-1, self.binCountArr(k.arr, schTpl))) * k.w
                                            walTouch.append(k)
                                    init = True
#                                    scan = scan + 1
                                else:
                                    for k in walTouch:
                                        schFit = schFit + int(math.pow(-1, self.binCountArr(k.arr, schTpl))) * k.w
#                                    reuse = reuse + 1 

                                schFitArr.append(Struct(fit=schFit,arr=schTpl))
                                #print subBitIn, j, schFit
#                            print 

                            schFitArrSort = sorted(schFitArr, key = lambda i: i.fit)

                            # perform voting from the best hyperplane associated with the subfunction
                            #for k in range(self.k+1):
                            for k in range(1):
                            #for k in range(self.k*2):
                                for j in subBitIn:
                                    if j in schFitArrSort[k].arr:
                                        #self.sumFitA[j].one = self.sumFitA[j].one + schFitArrSort[k].fit
                                        self.sumFitA[j].one = self.sumFitA[j].one + 1
                                    else:
                                        #self.sumFitA[j].zero = self.sumFitA[j].zero + schFitArrSort[k].fit
                                        self.sumFitA[j].zero = self.sumFitA[j].zero + 1

    def genHyperWalVote(self):
        """
        using the voting strategy where only best hyperplane have the chance to vote
        selecting hyperplane template on the basis of nonzero Walsh coefficients
        """
        self.transWal()

        # initialize sumFitA 
        self.sumFitA = []
        evalSubFunc = []
        for i in range(self.n):
            self.sumFitA.append(Struct(one=0,zero=0))
        
        for i in self.WA:
            subBit = i.arr

            if subBit not in evalSubFunc and i.arr:
                evalSubFunc.append(subBit)

                # check every template that matches the subfunction
                seqBits = tl.genSeqBits(len(subBit))
                schFitArr = []
                for j in seqBits:
                    schFit = 0

                    # convert bit string to array representation
                    schTpl = []
                    for k in range(len(j)):
                        if j[k] == '1':
                            schTpl.append(subBit[k])

                    # compute schema fitness
                    for k in self.WA:
                        subset = True
                        for l in k.arr:
                            if l not in subBit:
                                subset = False
                                break
                        if subset == True:
                            schFit = schFit + int(math.pow(-1, self.binCountArr(k.arr, schTpl))) * k.w

                    schFitArr.append(Struct(fit=schFit,arr=schTpl))
#                    print subBit, j, schFit
#                print 

                schFitArrSort = sorted(schFitArr, key = lambda i: i.fit)

                # perform voting from the best hyperplane associated with the subfunction
                #for k in range(self.k+1):
                for k in range(1):
                #for k in range(self.k*2):
                    for j in subBit:
                        if j in schFitArrSort[k].arr:
                            #self.sumFitA[j].one = self.sumFitA[j].one + schFitArrSort[k].fit
                            self.sumFitA[j].one = self.sumFitA[j].one + 1
                        else:
                            #self.sumFitA[j].zero = self.sumFitA[j].zero + schFitArrSort[k].fit
                            self.sumFitA[j].zero = self.sumFitA[j].zero + 1


    def composeFullStr(self, i, j, interBits, n):
        """ return the integer representation of Full String """
        subStr = bin(j)
        subStr = subStr[2:]
        if len(subStr) < self.k+1:
            subStr = '0'*(self.k+1-len(subStr)) + subStr
        indexSubOneBit = self.indexOneBit(subStr)

    def composeFullStr(self, i, j, interBits, n):
        """ return the integer representation of Full String """
        subStr = bin(j)
        subStr = subStr[2:]
        if len(subStr) < self.k+1:
            subStr = '0'*(self.k+1-len(subStr)) + subStr
        indexSubOneBit = self.indexOneBit(subStr)
        iStr = ['0']*n
        for k in range(len(indexSubOneBit)):
            iStr[int(interBits[indexSubOneBit[k]])] = subStr[indexSubOneBit[k]]
        iStr = ''.join(iStr)
        return int(iStr, 2)

    def composeFullBitStr(self, i, j, interBits, n):
        """ return the original full string """
        subStr = bin(j)
        subStr = subStr[2:]
        if len(subStr) < self.k+1:
            subStr = '0'*(self.k+1-len(subStr)) + subStr
        indexSubOneBit = self.indexOneBit(subStr)
        iStr = ['0']*n
        for k in range(len(indexSubOneBit)):
            iStr[int(interBits[indexSubOneBit[k]])] = subStr[indexSubOneBit[k]]
        iStr = ''.join(iStr)
        return iStr

    def checkInclude(self, occurOneBit, mask):
        for i in range(len(occurOneBit)):
            if occurOneBit[i] not in mask:
                return False
        return True

    def indexOneBit(self, iStr):
        range1 = range(len(iStr))
        return [ i for i in range1 if iStr[i] == '1']
            
    def dispNK(self):
        print self.n, self.k

    def binCountArr(self, a1, a2):
        """
        count the number of one bits appearing in both a1 and a2
        """
        s = 0
        for i in a1:
            if i in a2:
                s = s + 1

        return s
    def genSolProp(self, sumFitA):
        sol = []
        for i in range(self.n):
            if random.random() < sumFitA[i].zero / (sumFitA[i].one + sumFitA[i].zero + 0.0):
                sol.append('0')
            else:
                sol.append('1')
        return sol

    def countInterBits(self):
        """
        count the number of subfunctions that touch a particular bit
        """
        self.interBit = dict()

        for i in range(self.n):
            sub = self.neighs[i][:]
            sub.append(i)
            for j in range(len(sub)):
                for k in [h for h in range(len(sub)) if h > j]:
                    if sub[j] not in self.interBit:
                        self.interBit[sub[j]] = [sub[j],sub[k]]
                    elif sub[k] not in self.interBit[sub[j]]:
                        self.interBit[sub[j]].append(sub[k])

                    if sub[k] not in self.interBit:
                        self.interBit[sub[k]] = [sub[j],sub[k]]
                    elif sub[j] not in self.interBit[sub[k]]:
                        self.interBit[sub[k]].append(sub[j])

    def initInter(self):
        self.Inter = dict()
        for i in range(len(self.WA)):
#            W = int(math.pow(-1, self.binCount(self.WA[i].arr, self.bit))) * self.WA[i].w
#            self.WAS[i] = Struct(arr = self.WA[i].arr, w = W)
#            comb = self.genComb(len(self.WA[i].arr))

            for j in self.WA[i].arr:
#                self.sumArr[j] = self.sumArr[j] + W
                if len(self.WA[i].arr)>1: # for at least order Walsh terms
                    if j not in self.Inter: # the entry of i doesn't exist yet
                        self.Inter[j] = Struct(arr=Set(), WI=Set())

                    for k in self.WA[i].arr:
                        if k != j:
                            self.Inter[j].arr.add(k)
                    self.Inter[j].WI.add(i)

class Struct:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
