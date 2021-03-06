# generate NK-landscapes instances
## the encoding is in right to left fashion
## 00-0, 01-1, 10-2, 11-3 import WalshAnalysis as wal import random import numpy as np import math

import random
import time 
import sys
# sys.path.append("../fht/build/lib.linux-x86_64-2.7/fht/")

import fht

import numpy as np
import WalshAnalysis as wal
cimport numpy as np
import tool as tl
import math
from sets import Set
from libc.stdlib cimport malloc, free
from libcpp.set cimport set
import cython 
from cython.operator cimport dereference as deref, preincrement as inc
cimport cython

import pdb

ctypedef struct ComArr:
    int** arr
    int size

ctypedef struct InTer:
    set[int]* arr

    
cdef class NKLandscape:
    cdef ComArr** lookup
    # cdef InTer** Inter
    cdef public list Inter # the list of variables that interact with
                                          # ith variable
    cdef public int n                            # number of variables
    cdef public int k                            
    cdef public int c                            # number of clauses
    cdef public list neighs
    cdef public list f
    cdef float** func
    cdef list Kbits
    cdef public dict w
    cdef public list WA
    # cdef public list U
    cdef public dict U
    cdef public list listSubFunc
    # cdef double* subFit
    
    """ NK-landscape class """
    def __init__(self,inN, inK, inC, fileName = None):
        self.n = inN
        self.k = inK
        self.c = inC
        
        # self.m = self.n * self.m
        # print 'init NK'
        
        if fileName == None:
            print 'generating problem'
            self.genNeigh()
            self.genFunc()
            print 'done'
        else:
            # print fileName
            self.readFile(fileName)
            
        self.Kbits = tl.genSeqBits(self.k+1)
        self.lookup = <ComArr**> malloc(sizeof(ComArr*)*self.n)
        for i in xrange(self.n):
            self.lookup[i] = NULL

        # self.subFit = <double*> malloc(sizeof(double)*self.c) 

    # cpdef double getSubFitArr(self, i):
    #     return self.subFit[i]
        
    def exportToFile(self, fileName):
        # print 'fileName', fileName
        # print self.neighs
        # print len(self.neighs)
        f = open(fileName, 'w')
        for i in xrange(self.c):
            # print i
            for j in xrange(len(self.neighs[i])):
                print >>f, self.neighs[i][j], '\t',
            print >>f
        for i in xrange(self.c):
            for j in xrange(len(self.f[i])):
                print >>f, self.f[i][j], '\t',
            print >>f

    def readFile(self, fName):
        cdef int i
        cdef list f
        
        self.neighs = np.genfromtxt(fName, delimiter="\t", dtype='int', skip_footer=self.c, autostrip=True, usecols = range(self.k+1)).tolist()
        # sort the neighs, and overwrite self.neighs 
        for i in xrange(self.c):
            interbit = self.neighs[i][:]
            interbit.sort()
            self.neighs[i][:] = interbit
        
        self.f = np.genfromtxt(fName, delimiter="\t", skip_header=self.c, autostrip=True, usecols = range(int(math.pow(2,self.k+1)))).tolist()
        self.convertFunc()


    cdef convertFunc(self):
        """
        convert list of lists f to float matrix self.func
        """
        # allocate space for the matrix 
        self.func = <float **>malloc(sizeof(float *) * len(self.f))
        for i in xrange(len(self.f)):
            self.func[i] = <float *> malloc(sizeof(float) * len(self.f[i]) )
            for j in xrange(len(self.f[i])):
                self.func[i][j] = self.f[i][j]
                # print self.func[i][j]

    cpdef float getFuncVal(self,i, j):
        return self.func[i][j]
    
                

    def convertFuncDict(self):
        """
        the second dimension of self.func is a dict with the indices for subfunction input

        convert self.func into its dict() representation
        """
        print 'readFileDict'
        
    
    def genNeigh(self):
        """ 
        generate neighborhood of K+1 randomly picked positions out of N ones 
        self.neighs include the ith variable itself
        """
        self.neighs = []
        print 'genNeigh'
        if self.c == self.n:              # strictly NK-landscapses
             # print self.c, self.n
            for i in xrange(self.c):
                oneNeigh = random.sample(range(self.n-1),self.k)
                # print 'before', oneNeigh
                for j in xrange(len(oneNeigh)):
                    if oneNeigh[j] >= i:
                        oneNeigh[j] = oneNeigh[j] + 1
                oneNeigh.append(i)
                # print 'after', oneNeigh
                oneNeigh.sort()
                self.neighs.append(oneNeigh)
        else:                             # MaxSAT instances
            for i in range(self.c):
                oneNeigh = random.sample(range(self.n), self.k+1)
                oneNeigh.sort()
                self.neighs.append(oneNeigh)
        
    def getNeigh(self):
        return self.neighs

    def genFunc(self):
        """ generate function value """
        # print 'genFunc'
        self.f = []
        for i in xrange(self.c):
            oneFunc = []
            for j in xrange(int(math.pow(2,self.k+1))):
                oneFunc.append(random.random())
            self.f.append(oneFunc)
        self.convertFunc()

    def getFunc(self):
        return self.f

    def getN(self):
        return self.n

    def genK(self):
        return self.k

    cpdef float compFit(self, bitStr):
        """ compute the fitness value"""
        #       print bitStr
        cdef int  i
        cdef float s = 0.0
        
        # self.subFit = np.empty(self.c)
        
        for i in range(self.c):
            """ compose interacting bits """
            interBit = self.neighs[i][:]
            # don't sort the list every time
            # interBit.sort() 
            """ extract corresponding bits """
            bits = [ bitStr[int(j)] for j in interBit ]
            interStr = ''.join(bits)
            """ sum up the sub-function values """
            #            print self.func[i][int(interStr,2)]
            # self.subFit[i] = self.func[i][int(interStr,2)]
            s = s + self.func[i][int(interStr,2)]
        # print sum/float(self.c)
        # print
        # return sum/float(self.c)
        # print('sum')
        return s

    def destFunc(self):
        for i in xrange(self.c):
            free(self.func[i])
        free(self.func)
        

    
    cpdef float compSubFit(self, bitStr, i):
        """
        compute i^th subfunction fitness value according to bitStr
        """ 
        
        # print 'bitStr', bitStr
        # time.sleep(2)
        
        """ compose interacting bits """
        interBit = self.neighs[i][:]
        # print 'interbit'
        # don't sort the list every time
        # interBit.sort()
        """ extract corresponding bits """
        bits = [ bitStr[int(j)] for j in interBit ]
        # print bits
        interStr = ''.join(bits)
        """ sum up the sub-function values """
        # print 'interStr',interStr
        
        return self.func[i][int(interStr,2)] 
        
    # cpdef float sumTerm(self, list bitStr, int i, int p, int q):
    #     """
    #     compute the inner summation term 
        
    #     f_i(x^{qp}) - f_i(x^{q}) - f_i(x^{p}) + f_i(x)
    #     """
    #     cdef float s = 0.0
    #     cdef int j
    #     cdef list bits
    #     # cdef str pb, qb

    #     """ extract corresponding bits """
    #     # x
        
    #     # O(k), according to %timeit
    #     bits = [ bitStr[j] for j in self.neighs[i][:] ]
    #     # print bits
    #     # O(k), %timeit int(''.join(['1']*10000),2)
    #     s = s + self.func[i][int(''.join(bits),2)]   
        
    #     # p
    #     if bitStr[p] == '0':
    #         bitStr[p] = '1'
    #         pb = '0'
    #     else:
    #         bitStr[p] = '0'
    #         pb = '1'
    #     bits = [ bitStr[j] for j in self.neighs[i][:] ]
    #     s = s - self.func[i][int(''.join(bits),2)] 
        
    #     # pq
    #     if bitStr[q] == '0':
    #         bitStr[q] = '1'
    #         qb = '0'
    #     else:
    #         bitStr[q] = '0'
    #         qb = '1'
    #     bits = [ bitStr[j] for j in self.neighs[i][:] ]
    #     s = s + self.func[i][int(''.join(bits),2)] 

    #     # q
    #     bitStr[p] = pb
    #     # if bitStr[p] == '0':
    #     #     bitStr[p] = '1'
    #     # else:
    #     #     bitStr[p] = '0'
    #     bits = [ bitStr[j] for j in self.neighs[i][:] ]
    #     s = s - self.func[i][int(''.join(bits),2)] 
        
    #     # back to x
    #     bitStr[q] = qb
    #     # if bitStr[q] == '0':
    #     #     bitStr[q] = '1'
    #     # else:
    #     #     bitStr[q] = '0'
        
    #     return s
    
    
        
    def WalCof(self):
        """ compute the Walsh coefficients """
        subW = [] # subW is a N*2^K matrix
        for i in range(self.n):
            """ 1. Compute coefficients for each sub-functions """
            subWone = wal.computeW(self.Kbits, self.f[i])
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
            subWone = wal.computeW(self.Kbits, self.f[i])
            subW.append(subWone)
        # print 'len', math.pow(2,self.n)
        w = np.zeros(math.pow(2,self.n))
        for i in range(self.n): # i: index of subfunction
            interBits = self.neighs[i][:]
            interBits.append(i)
            interBits.sort()
            for j in range(int(math.pow(2, self.k+1))): # j: index of substrings
                indexW = self.composeFullStr(i, j, interBits, self.n)
                w[indexW] = w[indexW] + subW[i][j]

        return w/float(self.n)

    # @cython.boundscheck(False)
    def WalshCofLinearLinklist(self):
        """ compute the Walsh Coefficients in a liner time with linear space using dict() """

        # print('WalshCofLinearLinklist')
        # cdef list subW, subWone
        # cdef dict w
        # cdef tuple indexW
        # cdef np.ndarray masks = np.array([False]*(self.k+1), dtype=bool)
        # cdef int i, j
        # t0 = 0.0
        # t1 = 0.0
        # t2 = 0.0
        # t3 = 0.0
        
        subW = [] # subW is a N * 2^K matrix

        # start = time.time() 
        """ Compute coefficients for each sub-functions """
        # print self.k+1
        # print math.pow(2,(self.k+1)/2.0)
        for i in xrange(self.c):
            # subWoneF = wal.computeW(self.Kbits, self.f[i])
            # subWone = fht.fht(np.asarray(self.f[i]))/(2*math.sqrt(2))
            subWone = fht.fht(np.asarray(self.f[i]))/math.pow(2,(self.k+1)/2.0)
            # print subWone
            # print subWoneF
            # print
            subW.append(subWone)
        # t0 = t0 + time.time() - start
        
        """ use dict to represent all non-zero Walsh Coefficients"""
        w = dict()
        for i in xrange(self.c): # i: index of sub-function
            interBits = self.neighs[i]
            interBits = np.array(interBits)
            # print(type(interBits))
            
            masks = np.array([False]*(self.k+1), dtype=bool)   # initialize mask array 
            
            # print
            # print interBits
            for j in xrange(int(math.pow(2, self.k+1))): # j: index of substrings
                # indexW = interBits[ii] for ii in masks
                # start = time.time()
                indexW =tuple(interBits[masks])
                # print indexW
                # t1 = t1 + time.time() - start

                # start = time.time()
                # print 'indexW', indexW, type(indexW)
                # indexW = self.composeFullBitStr(i, j, interBits, self.n)
                if w.has_key(indexW):
                    w[indexW] = w[indexW] + subW[i][j]
                else:
                    w[indexW] = subW[i][j]
                # t2 = t2 + time.time() - start
                    
                # print(masks)
                # start = time.time()
                self.addOne(masks)
                # t3 = t3 + time.time() - start
                # print(masks), '\n'
            
            # exit()
        # for k in w.keys():
        #     w[k] = w[k]/float(self.n)
        # print t0, t1, t2, t3
        
        self.w = w
        return w

    cdef addOne(self, a):
        """
        add one to the boolean numpy array, a
        """
        cdef int i
        for i in xrange(len(a)-1,-1,-1):
             # range(6,0,-1)
            if a[i]==False:
                a[i] = True
                break
            else:
                a[i] = False
        
    """
    for Walsh Local Search
    """
    def transWal(self):
        """
        translate bitstring represented Walsh terms into arrays of bits that they touches
        """
        # t=time.time()
        self.WA = [] # array representing Walsh terms
        # print time.time()-t

        # print 'WA len', len(self.WA)
        # t=time.time()
        for k in self.w.keys():
            if self.w[k] != 0:
                self.WA.append( Struct(arr = list(k), w = self.w[k]))
                
        # print time.time()-t
        # print 'dict len', len(self.w), 'WA len', len(self.WA)
        # # density of Non-Zero Walsh terms
        # print len(self.WA), self.n * 2**(self.k+1), (len(self.WA)+0.)/(self.n * 2**(self.k+1))

        
    def clearWalModel(self):
        del self.WA[:]
        
    def genHyperVote(self):
        """
        using the voting strategy where only best hyperplane have the chance to vote
        """
        self.transWal()
#        bit,fit = tl.compFit(self)
#        a = sorted(zip(bit,fit), key=lambda a_entry: a_entry[1])
#        optBit = a[0][0]
#        optFit = a[0][1]
#        print 'opti\n',optBit, optFitee

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
        for k in xrange(len(indexSubOneBit)):
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

        for i in range(self.c):
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

    # def initInter(self):
    #     """ 
    #     old implementation used by walsh search
    #     """
    #     self.Inter = dict()
    #     for i in range(len(self.WA)):
    #         for j in self.WA[i].arr:
    #             if len(self.WA[i].arr)>1: # for at least order Walsh terms
    #                 if j not in self.Inter: # the entry of i doesn't exist yet
    #                     self.Inter[j] = Struct(arr=Set(), WI=Set())
    #                 for k in self.WA[i].arr:
    #                     if k != j:
    #                         self.Inter[j].arr.add(k)
    #                 self.Inter[j].WI.add(i)
    
    cpdef genInter(self):
        """ 
        - initialization of interaction information
        - generate self.Inter
        - the list of variables interacting with i if any, other self.Inter[i]==NULL
        """
        cdef list i
        cdef int j
        cdef int j0, j1
        cdef InTer* inter
        cdef ComArr* comb
        cdef InTer** CInter
        cdef set[int].iterator it

        # print 'init'
        Cinter = < InTer** > malloc(sizeof(void *)*self.n)
        for j in xrange(self.n):
            Cinter[j] = NULL

        self.lookup = <ComArr**> malloc(sizeof(ComArr*)*self.n)
        for j in xrange(self.n):
            self.lookup[j] = NULL
            
        # print 'gen CInter'
        # merely the function itself, check every pair of 
        for i in self.neighs:
            # generete all possible pairs
            # print 'genComb', len(i)
            comb = self.genComb(len(i))
            # print 'done genComb', len(i)
                        
            for j in xrange(comb.size):
                # print 'j'
                j0 = i[comb.arr[j][0]]
                j1 = i[comb.arr[j][1]]
                
                if Cinter[j0] == NULL:
                    inter = <InTer*> malloc(sizeof(InTer))
                    inter[0].arr = new set[int]()
                    Cinter[j0] = inter

                if Cinter[j1] == NULL:
                    inter = <InTer*> malloc(sizeof(InTer))
                    inter[0].arr = new set[int]()
                    Cinter[j1] = inter

                Cinter[j0].arr.insert(j1)
                Cinter[j1].arr.insert(j0)

        # print 'copy'
        # copy data to python variable for future use in other moduler
        # initialize
        self.Inter = [None]*self.n
        for j in xrange(self.n):
            if Cinter[j] != NULL:
                self.Inter[j] = []
                it = Cinter[j].arr.begin()
                while it!=Cinter[j].arr.end():
                    self.Inter[j].append(deref(it))
                    inc(it)

        # print 'free'
        # free space associated with Inter
        for j in xrange(self.n):
            if Cinter[j] != NULL:
                free(Cinter[j][0].arr)
                free(Cinter[j])
                Cinter[j] = NULL
        free(Cinter)

    def printInter(self):
        print 'full Inter'
        for i in xrange(self.n):
            print self.Inter[i]
        print
        
    # cdef getInter(self, i):
    #     """
    #     return the vector that contains all subfunction 9that contain ith variable
    #     """
    #     return self.Inter[i]

    def genListSubFunc(self):
        """
        - generate the list of subfunctions that interact with ith variable
        - self.listSubFunc as output
        """
        self.listSubFunc = [None] * self.n
        # go through all subfunction indices
        for i in xrange(len(self.neighs)):
            indices = self.neighs[i]
            for j in indices:
                if self.listSubFunc[j] == None:
                    self.listSubFunc[j] = []
                self.listSubFunc[j].append(i)

    def printListSubFunc(self):
        print 'listSubFunc'
        for i in self.listSubFunc:
            print i

    cpdef genUdict(self):
        """
        U array will be sparse, initializing the whole U is a huge
        waste. The indices of subfunction that contains both i and j.
        Assume i<=j.
        """
        cdef int i,j,j0,j1
        
        cdef int c=0

        # t0 = 0
        # t1 = 0 
        # t2 = 0
        # t3 = 0
        
        self.U = dict()
        for i in xrange(len(self.neighs)):   # for each subfunction $O(n)$
            comb = self.genComb(len(self.neighs[i]))
            for j in xrange(comb.size):
                
                # start = time.time()
                c = c + 1
                j0 = self.neighs[i][comb.arr[j][0]] 
                j1 = self.neighs[i][comb.arr[j][1]] 
                # t0 = t0 + time.time() - start

                # start = time.time()
                (j0,j1) not in self.U
                # t1 = t1 + time.time() - start
                
                if (j0,j1) not in self.U:
                    # start = time.time()
                    self.U[(j0,j1)] = [i]
                    # t2 = t2 + time.time() - start
                else:
                    # start = time.time()
                    self.U[(j0,j1)].append(i)
                    # t3 = t3 + time.time() - start
                    
        # print 'c', c, 't0', t0, 't1', t1, 't2', t2, 't3', t3

        
    # cpdef genU(self):
    #     """
    #     construct a matrix U for updating purposes. where only i<=j
    #     entries are non-empty. Each non-empty entry is a list of indices
    #     of subfunctions that includes both i and j. The list can be
    #     constructed on the basis of 'neighs' data structur. Using
    #     one-dimensional representation to avoid the odds of generating
    #     matrix.
        
    #     Note: need an analogy implementation for Walsh-based search, the
    #     initialization is stupid

    #     should probably use (i,j) as dictionary key
    #     """
    #     cdef int i, j, l, j0, j1, pos
    #     cdef list neigh
    #     cdef ComArr* comb
    #     cdef InTer* inter
    #     cdef InTer** u
    #     cdef set[int].iterator it
        
    #     # initialize U matrix (in 1-D representation)
    #     l = (self.c * (self.c-1))/2
    #     # print 'l', l, 'c', self.c
    #     u = < InTer** > malloc( sizeof(void *) * l )
        
    #     for i in xrange(l):
    #         u[i] = NULL

    #     # print 'init'
    #     # generate entries for U matrix based on neighs
    #     for i in xrange(len(self.neighs)):
    #         comb = self.genComb(len(self.neighs[i]))
    #         for j in xrange(comb.size):
    #             j0 = self.neighs[i][comb.arr[j][0]] 
    #             j1 = self.neighs[i][comb.arr[j][1]] 
    #             # calculate the position in U
    #             if j0<=j1:
    #                 pos = j1*(j1 - 1)/2 + j0
    #             else:
    #                 pos = j0*(j0 - 1)/2 + j1
    #             # print 'pos', pos
                
    #             if u[pos] == NULL:
    #                 inter = <InTer*> malloc(sizeof(InTer))
    #                 inter[0].arr = new set[int]()
    #                 u[pos] = inter

    #             #  variable j0 and j1 appear in ith subfunction
    #             u[pos].arr.insert(i)
                
    #     # print 'copy'        
    #     # copy data C space to Python space        
    #     self.U = [None]*l
    #     for i in xrange(l):
    #         self.U[i] = []
    #         if u[i] != NULL:
    #             it = u[i].arr.begin()  
    #             while it != u[i].arr.end():
    #                 self.U[i].append(deref(it))  
    #                 inc(it)

    #     # print 'free'
    #     # free memory
    #     for i in xrange(l):
    #         if u[i] != NULL:
    #             free(u[i].arr)
    #             free(u[i])
    #     free(u)
        
    cpdef list getU(self, int i, int j):
        """
        return the list of clause indices that touches both i and j
        """
        cdef int pos
        if i<=j:
            return self.U[(i,j)]            
        else:
            return self.U[(j,i)]            
        
    # def printU(self):
    #     """
    #      print self.U for verifying purposes
    #     """
    #     cdef set[int].iterator it
    #     cdef int ii, i, j, pos

    #     print 'printU'
    #     for i in xrange(self.n):
    #         for j in xrange(i):
    #             pos = i*(i - 1)/2 + j
    #             print '*********'
    #             print 'i', i, 'j', j, 'pos', pos

    #             if self.U[pos] != NULL:
    #                 it = self.U[pos].arr.begin()
    #                 while it != self.U[pos].arr.end():
    #                     ii = deref(it)
    #                     print ii
    #                     inc(it)
    #             else:
    #                 print i, j, 'not'
            

    # cdef ComArr* genComb(self,int N) nogil:
    cdef ComArr* genComb(self,int N):
        """
        Generate C_N^2 sequence, index are stored, because they are more general, Implemented in an *incremental* fashion.
        """
        cdef int c, j, i, counter
        cdef ComArr* ptr

        if self.lookup[N] != NULL: # the key exists before
            # print 'exist'
            return self.lookup[N]
        else : # the key not found
            # print 'not found', 'N',  N
            
            # c = biomial(N, 2)
            c = N*(N-1)/2
            # print 'biomial', c
            counter = 0

            # print 'allocate'
            ptr = <ComArr*> malloc(sizeof(ComArr))
            ptr.arr = <int**> malloc(sizeof(int*)*c)
            ptr.size = c
            
            # print 'init'
            for i in xrange(c):
                ptr.arr[i] = <int*> malloc(sizeof(int)*2)

            # print 'assign'
            for i in xrange(N):
                for j in xrange(i+1,N):
                    ptr.arr[counter][0] = i
                    ptr.arr[counter][1] = j
                    counter = counter + 1
                    
            # print 'ptr'
            self.lookup[N] = ptr
            return ptr

            
    def countFreqInFunc(self):
        """ count the frequency of a variable i appearing in a function """
        freq = np.zeros(self.n)
        for i in self.neighs:
            for j in i:
                freq[j] = freq[j] + 1

        return freq


    def __del__(self):
        """ 
        destructor to be called automatically when instance is about to be destroyed
        """
        cdef int i,j

        # free(self.subFit)
        
        for i in xrange(self.n):
            free(self.lookup[i])
        free(self.lookup)

        # for j in xrange(self.n):
        #     if self.Inter[j] != NULL:
        #         free(self.Inter[j][0].arr)
        #         free(self.Inter[j])
        #         self.Inter[j] = NULL
        # free(self.Inter)

        
        # print 'del nklandscape' 


# @cython.cdivision(True)
# # cdef int biomial(int N, int K) nogil:
# cdef int biomial(long N, long K):
#     """ compute the combination of N choose K """
#     print N
#     print factorial(N)
#     print factorial(K)
#     print factorial(N-K)
#     return factorial(N)/( factorial(K) * factorial(N-K) )


# cdef int factorial(long N) nogil:
#     """ compute N! """
#     cdef int c, fact = 1
#     for c in xrange(1,N+1):
#         fact = fact*c
#     return fact

        
class Struct:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

cdef class NonNKLandscape(NKLandscape):
    """ Non-uniform random NK-Landscapes """
    def __init__(self, inN, inK, inC, fileName = None):
        self.n = inN
        self.k = inK
        self.c = inC
        NKLandscape.__init__(self, inN, inK, inC, fileName)

        if fileName == None:
            self.genNonNeigh() # clear from parental class, re-generate
        else:
            self.readFile(fileName)

           
        self.Kbits = tl.genSeqBits(self.k+1)
        # print 'after init', self.neighs
        
    def genNonNeigh(self):
        self.neighs = []
        # enforce n=c
        if self.c == self.n:
            # print 'NonNK', 'N', self.n, 'C', self.c
            
            # generate permutation
            pi = range(self.n)
            random.shuffle(pi)
            # print pi
            for i in range(self.n):
                # print 'i', i
                oneNeigh = [i]
                while len(oneNeigh)<self.k+1:
                    index = np.random.binomial(len(pi), 0.5)
                    # print 'index', index
                    rand = pi[index]
                    if rand not in oneNeigh:
                        oneNeigh.append(rand)
                        # print 'rand', rand
                self.neighs.append(oneNeigh)
            # print 'neighs', self.neighs
        else:
            print 'Non-uniform NK(q) landscape instances require c==n'
