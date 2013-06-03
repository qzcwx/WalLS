# carry the Walsh analysis information with individual
import numpy as np
import math
import random
import copy
import time
from sets import Set

# import for Cython
cimport cython
from libc.math cimport pow
from cpython cimport bool
from cython.parallel import prange, parallel, threadid
from libcpp.vector cimport vector
from libcpp.set cimport set
from libc.stdlib cimport malloc, free
from libc.string cimport strcpy
from cython.operator cimport dereference as deref, preincrement as inc

class orderSum:
    def __init_(self, k=0):
        self.totalSum = 0
        self.subSum = np.zeros(k)

class Struct:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

ctypedef struct InfBit:
    vector[int]* arr
    int WI

# ctypedef struct Uelem:                      # each element in the U array/matrix is a list
#     vector[int]* arr                      # indices of WAS terms
    
ctypedef struct InTer:
    set[int]* arr                         # list of variables that are interacting with p
    set[int]* WI                          # list of the addresses of Walsh terms that include p bit

ctypedef struct Was:
    int* arr
    float w

ctypedef struct ComArr:
    int** arr
    int size

cdef class Individual:
    cdef InTer** Inter                    # the list of variables that interact with the ith variable
    cdef vector[InfBit*]** infectBit
    cdef Was* WAS
    cdef float* SC
    cdef float* Z
    cdef float** C
    cdef float** orderC
    # cdef Uelem** U
    cdef dict U
    # cdef public list improveA
    # cdef set[int] improveA
    cdef public object improveA
    cdef public list improveSC
    cdef object func
    cdef object model
    cdef int MaxFit
    cdef public float threshold
    cdef int fitEval
    cdef ComArr** lookup 
    cdef public float fit
    cdef public float fitG
    cdef public list bit
    cdef float* sumArr
    cdef public int dim
    cdef public int addC
    cdef public int addWAS
        
    
    def __init__( self, n=0, neigh=False, oldIndiv=False ):
        # self.bit = NULL
        self.bit = []
        self.fit = 0
        self.dim = n
        if neigh == True:
            self.fitG = 0
        if oldIndiv != False:
            self.fit = oldIndiv.fit
            self.dim = oldIndiv.dim
            self.fitG = oldIndiv.fitG
            # print 'oldindiv.bit', oldIndiv.bit
            # time.sleep(2)
            # self.copyBit(oldIndiv.bit)
            self.bit = oldIndiv.bit[:]
        self.threshold = 1e-15
        self.addC = 0
        self.addWAS = 0
        
    # def __del__(self):
    #     # free(self.bit)
        
    def init(self):
        # self.initIndiv()
        self.initIndivPy()

    # def copyBit(self,  bit):
    #     # print  bit
    #     # time.sleep(2)
        
    #     # self.bit = new char[self.dim+1]()
        
    #     # self.bit = <char*>malloc( (self.dim+1) * sizeof(char))
    #     # self.bit=bit[:]
    #     # print bit
    #     # time.sleep(2)

    #     # print bit
    #     # for j in xrange(self.dim):
    #     #     randBitStr[j] = bit[j]
    #     for j in xrange(self.dim):
    #             self.bit[j] = bit[j]
    #     self.bit[self.dim] = '\0'
        
    #     # print 'randBitStr',randBitStr
    #     print bit
        
        # # time.sleep(2)
        # exit()
        
        # # print 'copybit', randBitStr
        
        # print 'self.bit',self.bit
        # time.sleep(2)

    def initWal(self, model):
        self.model = model

        cdef int i,j,k
        cdef float W
        # cdef InTer* inter
        cdef vector[InfBit*]* vectPtr
        cdef InfBit* strPtr
        cdef ComArr* comb
        cdef Was* was

        self.sumArr = <float*>malloc(self.dim * sizeof(float))
        for i in xrange(self.dim):
            self.sumArr[i] = 0

        self.infectBit = < vector[InfBit*]** > malloc(sizeof(void *) * self.dim)
        for i in xrange(self.dim):
            vectPtr = new vector[InfBit*]()
            self.infectBit[i] = vectPtr

        self.C = <float **>malloc(sizeof(float *) * self.dim)
        for i in xrange(self.dim) :
            self.C[i] = <float *> malloc(sizeof(float) * self.dim)
            for j in xrange(self.dim):
                self.C[i][j] = 0

        self.WAS = <Was* > malloc(sizeof(Was)* len(self.model.w.keys()))
        
        self.lookup = <ComArr**> malloc(sizeof(ComArr*)*self.dim)
        for i in xrange(self.dim):
            self.lookup[i] = NULL
            
        self.Inter = < InTer** > malloc(sizeof(void *)*self.dim)
        for i in xrange(self.dim):
            self.Inter[i] = NULL


        for i in xrange(len(self.model.WA)):
            W = int(math.pow(-1, binCount(self.model.WA[i].arr, self.bit))) * self.model.WA[i].w

            was = <Was *>malloc(sizeof(Was))
            was[0].arr = <int *>malloc(sizeof(int)*len(self.model.WA[i].arr))
            for j in xrange(len(self.model.WA[i].arr)):
                was[0].arr[j] = self.model.WA[i].arr[j]
            was[0].w = W
            self.WAS[i] = was[0]

            comb = self.genComb(len(self.model.WA[i].arr))

            for j in self.model.WA[i].arr:
                self.sumArr[j] = self.sumArr[j] + W
                if len(self.model.WA[i].arr)>1: # for at least order Walsh terms
                    #if not self.Inter[j]: # the entry of i doesn't exist yet
                    if self.Inter[j] == NULL:
                        inter = <InTer*> malloc(sizeof(InTer))
                        inter[0].arr = new set[int]()
                        inter[0].WI = new set[int]()
                        self.Inter[j] = inter

                    for k in self.model.WA[i].arr:
                        if k != j :
                            self.Inter[j].arr.insert(k)
                    self.Inter[j].WI.insert(i)

                # add list of order >= 3 Walsh terms for the purpose of updating C matrix
                if len(self.model.WA[i].arr) >= 3:
                    strPtr = <InfBit *> malloc(sizeof(InfBit))
                    strPtr.WI = i
                    strPtr.arr = new vector[int]()

                    for k in self.model.WA[i].arr:
                        strPtr.arr[0].push_back(k)
                    #self.infectBit[j][0].push_back(strPtr[0])
                    self.infectBit[j][0].push_back(strPtr)

            for l in xrange(comb.size):
                j0 = self.model.WA[i].arr[comb.arr[l][0]]
                j1 = self.model.WA[i].arr[comb.arr[l][1]]
                self.C[j0][j1] = self.C[j0][j1] + W

    def initWalU(self, model):
        
        self.model = model

        cdef int i,j,k,l,j0,j1,pos
        cdef float W
        cdef vector[InfBit*]* vectPtr
        cdef InfBit* strPtr
        cdef ComArr* comb
        cdef Was* was
        # cdef Uelem* uelem

        self.sumArr = <float*>malloc(self.dim * sizeof(float))
        for i in xrange(self.dim):
            self.sumArr[i] = 0
            
        self.WAS = <Was* > malloc(sizeof(Was)* len(self.model.w.keys()))
        self.lookup = <ComArr**> malloc(sizeof(ComArr*)*self.dim)
        for i in xrange(self.dim):
            self.lookup[i] = NULL

        self.Inter = < InTer** > malloc(sizeof(void *)*self.dim)
        for i in xrange(self.dim):
            self.Inter[i] = NULL
            
        # initialize U array
        
        # l = (self.model.c * (self.model.c-1))/2
        # self.U = < Uelem** > malloc( sizeof(void *) * l )
        # for i in xrange(l):
        #     self.U[i] = NULL
        self.U = dict()
        
        for i in xrange(len(self.model.WA)):
            # |Q|*n operation should be at least n^2
            W = int(math.pow(-1, binCount(self.model.WA[i].arr, self.bit))) * self.model.WA[i].w
            
            was = <Was *>malloc(sizeof(Was))
            was[0].arr = <int *>malloc(sizeof(int)*len(self.model.WA[i].arr))
            for j in xrange(len(self.model.WA[i].arr)):
                was[0].arr[j] = self.model.WA[i].arr[j]
            was[0].w = W
            self.WAS[i] = was[0]

            for j in self.model.WA[i].arr:
                self.sumArr[j] = self.sumArr[j] + W

            if len(self.model.WA[i].arr)>1: # for at least order Walsh terms
                for j in self.model.WA[i].arr:   # 
                    #if not self.Inter[j]: # the entry of i doesn't exist yet
                    if self.Inter[j] == NULL: 
                        inter = <InTer*> malloc(sizeof(InTer))
                        inter[0].arr = new set[int]()
                        inter[0].WI = new set[int]()
                        self.Inter[j] = inter
                        
                    for k in self.model.WA[i].arr:   # k^2 part
                        if k != j :
                            self.Inter[j].arr.insert(k)
                    self.Inter[j].WI.insert(i)
                    
                # add entries in U matrix, k^2 part 
                comb = self.genComb(len(self.model.WA[i].arr)) 
                # print self.model.WA[i].arr   
                for j in xrange(comb.size):
                    j0 = self.model.WA[i].arr[comb.arr[j][0]] 
                    j1 = self.model.WA[i].arr[comb.arr[j][1]] 

                    # calculate the position in U
                    if j0<=j1:
                        if (j0,j1) not in self.U:
                            self.U[(j0,j1)] = [i]
                        else:
                            self.U[(j0,j1)].append(i)
                    else:                 # j0>j1
                        if (j1,j0) not in self.U:
                            self.U[(j1,j0)] = [i]
                        else:
                            self.U[(j1,j0)].append(i)
                        
                        # print 'pos', pos
                    # print '(j0, j1)', j0, j1
                    # print 'pos',pos
                
        # print 'finish initWalU'
    
            #     # add list of order >= 3 Walsh terms for the purpose of updating C matrix
            #     if len(self.model.WA[i].arr) >= 3:
            #         strPtr = <InfBit *> malloc(sizeof(InfBit))
            #         strPtr.WI = i
            #         strPtr.arr = new vector[int]()

            #         for k in self.model.WA[i].arr:
            #             strPtr.arr[0].push_back(k)
            #         #self.infectBit[j][0].push_back(strPtr[0])
            #         self.infectBit[j][0].push_back(strPtr)

            # for l in xrange(comb.size):
            #     j0 = self.model.WA[i].arr[comb.arr[l][0]]
            #     j1 = self.model.WA[i].arr[comb.arr[l][1]]
            #     self.C[j0][j1] = self.C[j0][j1] + W


    # def printWalU(self):
    #     """
    #     print the U matrix associated with indices of Walsh terms
    #     """
    #     cdef vector[int].iterator it
    #     print 'printWalU'
        
    #     for i in xrange(self.dim):
    #         for j in xrange(i):
    #             pos = i*(i - 1)/2 + j
    #             print '*********'
    #             print 'i', i, 'j', j, 'pos', pos
                
    #             if self.U[pos] != NULL:
    #                 # for k in self.U[pos].arr:
    #                 #     print k
    #                 it = self.U[pos].arr.begin()
    #                 while it != self.U[pos].arr.end():
    #                     ii = deref(it)
    #                     print ii
    #                     inc(it)
    #             else:
    #                 print i, j, 'not'
            
    def printWAS(self):
        """
        print WAS array, all Walsh terms
        """
        print 'WAS'
        # print len(self.model.WA)
        for i in xrange(len(self.model.WA)):
            print 'i',i
            print 'arr',
            for j in xrange(len(self.model.WA[i].arr)):
                print self.WAS[i].arr[j],
            print
            print self.WAS[i].w
            print '***************'
        
    cpdef initBfUpdate(self, old, bool minimize, object model):
        """
        initialize data structures for performing partial update with S vector maintained
        """
        cdef int i,j
        cdef str oldI
        self.model = model
         

        # start = time.time() 
        # initialize S vector, first derivative
        self.sumArr = <float*>malloc(self.dim * sizeof(float))
        for i in xrange(self.dim):
            self.sumArr[i] = 0
        # print time.time() - start
            
        # start = time.time()
        for i in xrange(self.dim):
            # indiv = Individual(oldIndiv=old)
            # indiv.flip(i)
            
            for j in self.model.listSubFunc[i]:
                # get ride of compSubFit completely

                # subtract the old value
                bits = [old.bit[k] for k in self.model.neighs[j][:]]
                self.sumArr[i] = self.sumArr[i] - self.model.getFuncVal(j, int(''.join(bits),2))

                # flip i 
                if old.bit[i] == '0':
                    old.bit[i] = '1'
                    oldI = '0'
                else:
                    old.bit[i] = '0'
                    oldI = '1'
                    
                bits = [ old.bit[k] for k in self.model.neighs[j][:] ]
                self.sumArr[i] = self.sumArr[i] + self.model.getFuncVal(j, int(''.join(bits),2))
                old.bit[i] = oldI
                
                # self.sumArr[i] = self.sumArr[i] + (self.model.compSubFit(indiv.bit, j) - self.model.compSubFit(old.bit, j))
        # print time.time() - start

        # start = time.time() 
        self.genImproveSpartialUpdate(minimize)
        # print time.time() - start
        
        
            
    def checkWalshSum(self):
        """
        compute the sum of Walsh terms grouped by the order (1 to k)

        * compute the sum of order seperately for each possible bit-flip (k*n)
        """
        self.walSumArr = [orderSum(self.model.k) for i in range(self.dim)]
        

    def initSC(self):
        # compute the SC array
        # self.SC = np.zeros(self.dim)
        # self.Z = np.zeros(self.dim)
        # self.orderC = np.zeros((self.dim,self.dim))

        self.SC = <float *>malloc(self.dim* sizeof(float))
        self.Z = <float *>malloc(self.dim* sizeof(float))

        self.orderC = <float **>malloc(sizeof(float *) * self.dim)
        for i in xrange(self.dim) :
            self.orderC[i] = <float *> malloc(sizeof(float) * self.dim)

        for p in range(self.dim):
            phi = np.zeros(self.model.k+1)
            if p in self.model.Inter:
                for i in self.model.Inter[p].WI:
                    #order = len(self.WAS[i].arr)
                    order = len(self.model.WA[i].arr)
                    phi[order-1] = phi[order-1] + self.WAS[i].w

            self.Z[p] = self.sumArr[p]
            for i in range(1, self.model.k+1):
                if phi[i] != 0:
                    self.Z[p] = self.Z[p] + i * phi[i]

            self.SC[p] = self.sumArr[p] - 2/float(self.dim) * self.Z[p]

        for i in range(len(self.model.WA)):
            lenArr = len(self.model.WA[i].arr)
            comb = self.genComb(lenArr)
            for j in xrange(comb.size):
                j0 = self.model.WA[i].arr[comb.arr[j][0]]
                j1 = self.model.WA[i].arr[comb.arr[j][1]]
                self.orderC[j0][j1] = self.orderC[j0][j1] + lenArr * self.WAS[i].w

    ## def compPSum(self,bitStr):
    ##     """
    ##     use Elementary Landscape Analysis to obtain the average of neighs of given
    ##     individual
    ##     """
    ##     p = np.zeros(self.model.k+1)
    ##     for k in self.model.w.keys():
    ##         oneC = k.count('1')
    ##         if  oneC !=0 :
    ##             p[oneC-1] = p[oneC-1] + self.model.w[k] * math.pow(-1,wal.bc(k,bitStr))
    ## #            else :
    ## #                p[0] = p[0] + self.model.w[k]

    ##     s = 0
    ##     for i in range(self.model.k+1):
    ##         s = s + (i+1)*p[i]
    ##     return s


    def compPhisum(self,p):
        """
        \varphi_{p,i}^{\prime}(x) = \Sigma_{order j terms, that touches bit p}
        """
        phi = np.zeros(self.model.k+1)
        if p in self.model.Inter:
            for i in self.model.Inter[p].WI:
                order = len(self.mode.WA[i].arr)
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
            s = s + self.C[i][p]

        for i in range(p+1, self.dim):
            s = s + self.C[p][i]
        s = s + self.sumArr[p]

        return s

        
    cpdef update(self,int p):
        """
        By keeping track of coincidence matrix,
        Cij stands for S_i(y_j) = S_i(x) - C_ij
        partially update the Sum Array and self.WAS, given the bit which is changed
        """
        cdef int i,ii, k0, k1, k
        cdef int len1
        cdef set[int].iterator it
        cdef vector[int].iterator itt
        cdef vector[int] arr
        cdef InfBit I
        cdef ComArr* comb

        # print 'update'
        self.sumArr[p] = - self.sumArr[p]
        if self.Inter[p]!=NULL:
            """ iterative over self.Inter[p].arr """
            #print 'Not NULL', p
            #for i in prange(len(self.Inter[p].arr), nogil= True):
            it = self.Inter[p].arr.begin()
            #print 'begin'
            while it != self.Inter[p].arr.end():
                ii = deref(it)
                #print ii
                if ii < p:
                    self.sumArr[ii] = self.sumArr[ii] - 2*self.C[ii][p]
                    self.C[ii][p] = - self.C[ii][p]
                else:
                    self.sumArr[ii] = self.sumArr[ii] - 2*self.C[p][ii]
                    self.C[p][ii] = - self.C[p][ii]
                inc(it)

        # update the rest of elements in C matrix
        if self.infectBit[p].size() != 0:
            # for i in prange(self.infectBit[p].size(), nogil=True):
            for i in xrange(self.infectBit[p].size()):
                I = self.infectBit[p][0][i][0]
                arr = I.arr[0]
                itt = arr.begin()
                while itt != arr.end():
                    if deref(itt) == p:
                        arr.erase(itt)
                        break
                    inc(itt)
                comb = self.genComb(arr.size())
                for k in xrange(comb.size):
                    k0 = arr[int(comb.arr[k][0])]
                    k1 = arr[int(comb.arr[k][1])]
                    self.C[k0][k1] = self.C[k0][k1] - 2 * self.WAS[I.WI].w


                    
    cpdef updateU(self,int p):
        """
        partially update the sum array directly by refering to U _dict_,
        without C matrix

        complexity of ck(k-1)2^{k-2}
        """
        cdef int i,ii
        cdef float s
        cdef set[int].iterator it
        cdef list wList

        self.sumArr[p] = - self.sumArr[p]
        if self.Inter[p]!=NULL:
            """ iterative over self.Inter[p].arr """
            #for i in prange(len(self.Inter[p].arr), nogil= True):
            it = self.Inter[p].arr.begin()
            while it != self.Inter[p].arr.end():
                ii = deref(it)
                # compute the sum of walsh terms that touches both ii and p
                if ii<=p:
                    wList = self.U[(ii,p)]
                else:
                    wList = self.U[(p,ii)]
                    
                s = 0
                for i in wList:
                    s = s + self.WAS[i].w
                    self.addC = self.addC + 1
                    
                self.sumArr[ii] = self.sumArr[ii] - 2 * s
                inc(it)
                # if ii < p:
                #     self.sumArr[ii] = self.sumArr[ii] - 2*self.C[ii][p]
                #     self.C[ii][p] = - self.C[ii][p]
                # else:
                #     self.sumArr[ii] = self.sumArr[ii] - 2*self.C[p][ii]
                #     self.C[p][ii] = - self.C[p][ii]
                # inc(it)

                
    cpdef updateImprS(self, p, minimize):
        cdef int i
        cdef set[int].iterator it
        
        if p in self.improveA:
            self.improveA.remove(p)
            
            # while it!=self.model.Inter[p].arr.end():
        # if p in self.model.Inter:
        # if self.model.Inter[p]:
            # for i in self.model.Inter[p]:
        if self.Inter[p]!=NULL:
            it = self.Inter[p].arr.begin()
            while it != self.Inter[p].arr.end():
                i = deref(it)
                """ equal moves """
                # if (minimize == True and self.sumArr[i] > - self.threshold) or (minimize == False and self.sumArr[i] < self.threshold ):
                """ NOT equal moves """
                if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i]< - self.threshold ):
                    if i not in self.improveA:
                    # if self.improveA.empty() == False:
                        # self.improveA.append(i)
                        self.improveA.add(i)
                elif i in self.improveA:
                # elif i in 
                    self.improveA.remove(i)
                inc(it)
                
    cpdef updateImprSpartialUpdate(self, int p, minimize):
        cdef int i
        
        # cdef set[int].iterator it

        if p in self.improveA:
            self.improveA.remove(p)
            
            # while it!=self.model.Inter[p].arr.end():
        # if p in self.model.Inter:
        if self.model.Inter[p]!=None:
            for i in self.model.Inter[p]:
        # if self.Inter[p]!=NULL:
            # it = self.Inter[p].arr.begin()
            # while it != self.Inter[p].arr.end():
                # i = deref(it)
                """ equal moves """
                # if (minimize == True and self.sumArr[i] > - self.threshold) or (minimize == False and self.sumArr[i] < self.threshold ):
                """ NOT equal moves """
                if (minimize == False and self.sumArr[i] > self.threshold) or (minimize == True and self.sumArr[i]< - self.threshold ):
                    if i not in self.improveA:
                        # self.improveA.append(i)
                        self.improveA.add(i)
                elif i in self.improveA:
                    self.improveA.remove(i)
                # inc(it)

    
    cpdef updatePertImprS(self, p, minimize):
        cdef int i
        cdef set[int].iterator it

        # print 'updatePertImprS'
        if self.Inter[p]!= NULL:
        # if self.model.Inter[p]:
        # if p in self.model.Inter:
            # for i in self.model.Inter[p]:
            it = self.Inter[p].arr.begin()
            while it!=self.Inter[p].arr.end():
                i = deref(it)
                """ equal moves """
                # if (minimize == True and self.sumArr[i] > - self.threshold) or (minimize == False and self.sumArr[i]< self.threshold ):
                """ NOT equal moves """
                if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i]< - self.threshold ):
                    if i not in self.improveA:
                        self.improveA.add(i)
                elif i in self.improveA:
                    self.improveA.remove(i)
                inc(it)

        """ equal move """
        #if (minimize == True and self.sumArr[p] > - self.threshold) or (minimize == False and self.sumArr[p] < self.threshold ):
        """ NOT equal move """
        if (minimize == True and self.sumArr[p] > self.threshold) or (minimize == False and self.sumArr[p] < - self.threshold ):
            if p not in self.improveA:
                self.improveA.add(p)
        elif p in self.improveA:
            self.improveA.remove(p)

    cpdef updatePertImprSpartialUpdate(self, p, minimize):
        cdef int i
        # cdef set[int].iterator it

        # print 'updatePertImprS'
        # if self.Inter[p]!= NULL:
        if self.model.Inter[p] != None:
        # if p in self.model.Inter:
            for i in self.model.Inter[p]:
            # it = self.Inter[p].arr.begin()
            # while it!=self.Inter[p].arr.end():
                # i = deref(it)
                """ equal moves """
                # if (minimize == True and self.sumArr[i] > - self.threshold) or (minimize == False and self.sumArr[i]< self.threshold ):
                """ NOT equal moves """
                if (minimize == False and self.sumArr[i] > self.threshold) or (minimize == True and self.sumArr[i]< - self.threshold ):
                    if i not in self.improveA:
                        self.improveA.append(i)
                elif i in self.improveA:
                    self.improveA.remove(i)
                # inc(it)

        """ equal move """
        #if (minimize == True and self.sumArr[p] > - self.threshold) or (minimize == False and self.sumArr[p] < self.threshold ):
        """ NOT equal move """
        if (minimize == False and self.sumArr[p] > self.threshold) or (minimize == True and self.sumArr[p] < - self.threshold ):
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

        if p in self.model.Inter:
            for i in self.model.Inter[p].arr:
                if i < p:
                    self.sumArrFake[i] = self.sumArrFake[i] - 2*self.Cfake[i,p]
                else:
                    self.sumArrFake[i] = self.sumArrFake[i] - 2*self.Cfake[p,i]
                    
    cpdef updateWAS(self,int p):
        # complexity of k*2^{k-1}
        cdef int ii

        cdef set[int].iterator it
        
        # if p in self.model.Inter:
        if self.Inter[p]!=NULL:
            # for i in self.Inter[p].WI:
            it = self.Inter[p].WI.begin()
            while it != self.Inter[p].WI.end():
                ii = deref(it)
                self.WAS[ii].w = - self.WAS[ii].w
                self.addWAS = self.addWAS + 1
                inc(it)
                

    def updateSC(self, p):
        cdef int i,ii, k0, k1
        cdef int len1
        cdef set[int].iterator it
        cdef vector[int].iterator itt
        cdef vector[int] arr
        cdef InfBit I
        cdef ComArr* comb

        self.SC[p] = - self.SC[p]
        self.Z[p] = - self.Z[p]

        #update Z array
        if self.Inter[p]!=NULL:
            it = self.Inter[p].arr.begin()
            while it != self.Inter[p].arr.end():
                ii = deref(it)
                if ii < p :
                    self.Z[ii] = self.Z[ii]  - 2* self.orderC[ii][p]
                    self.orderC[ii][p] = - self.orderC[ii][p]
                else :
                    self.Z[ii] = self.Z[ii]  - 2* self.orderC[p][ii]
                    self.orderC[p][ii] = - self.orderC[p][ii]
                self.SC[ii] = self.sumArr[ii] - 2/float(self.dim) * self.Z[ii]
                inc(it)

        # update the rest of elements in C matrix
        if self.infectBit[p].size() != 0:
            for i in range(self.infectBit[p].size()):
                I = self.infectBit[p][0][i][0]
                arr = I.arr[0]
                itt = arr.begin()
                while itt != arr.end():
                    if deref(itt) == p:
                        arr.erase(itt)
                        break
                    inc(itt)
                comb = self.genComb(arr.size())
                for k in xrange(comb.size):
                    k0 = arr[int(comb.arr[k][0])]
                    k1 = arr[int(comb.arr[k][1])]
                    self.orderC[k0][k1] = self.orderC[k0][k1] - 2 * (arr.size()+ 1)* self.WAS[I.WI].w

    def updateImprSC(self, p, minimize):
        if p in self.improveSC:
            self.improveSC.remove(p)
        if p in self.model.Inter:
            for i in self.model.Inter[p].arr:
                """ equal move """
                #if (minimize == True and self.SC[i] > - self.threshold) or (minimize == False and self.SC[i] < self.threshold):
                """ NOT equal move """
                if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i] < - self.threshold):
                    if i not in self.improveSC:
                        self.improveSC.append(i)
                elif i in self.improveSC:
                    self.improveSC.remove(i)

    def updatePertImprSC(self, p, minimize):
        if p in self.model.Inter:
            for i in self.model.Inter[p].arr:
                """ equal move """
                #if (minimize == True and self.SC[i] > - self.threshold) or (minimize == False and self.SC[i] < self.threshold):
                """ NOT equal move """
                if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i] < - self.threshold):

                    if i not in self.improveSC:
                        self.improveSC.append(i)
                elif i in self.improveSC:
                    self.improveSC.remove(i)
        """ equal move """
        #if (minimize == True and self.SC[p] > - self.threshold) or (minimize == False and self.SC[p] < self.threshold):
        """ not equal move """
        if (minimize == True and self.SC[p] > self.threshold) or (minimize == False and self.SC[p] < - self.threshold):
            if p not in self.improveSC:
                self.improveSC.append(p)
        elif p in self.improveSC:
            self.improveSC.remove(p)

    def genImproveS(self,minimize):
        """
        generate the index of best neigh according to sumArr only (surrogate of fitness)
        """
        # check improving move
        # self.improveA = []
        self.improveA = Set()

        for i in range(self.dim):
            if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i]< - self.threshold):            # NOT equal moves
            # if (minimize == True and self.sumArr[i] > - self.threshold) or (minimize == False and self.sumArr[i] < self.threshold):              equal move
                self.improveA.add(i)
                
    def genImproveSpartialUpdate(self,minimize):
        """
        generate the index of best neigh according to sumArr only (surrogate of fitness)
        """
        # check improving move
        # self.improveA = new set[int]()
        # self.improveA = []
        self.improveA = Set()
        
        for i in range(self.dim):
            # print self.sumArr[i]
            if (minimize == True and self.sumArr[i] < - self.threshold) or (minimize == False and self.sumArr[i] >  self.threshold):            # NOT equal moves
            # if (minimize == True and self.sumArr[i] > - self.threshold) or (minimize == False and self.sumArr[i] < self.threshold):              equal move
                # self.improveA.append(i)
                self.improveA.add(i)
                # print 'add', i

        # print 'set', self.improveA
        
##     def genFitNext(self,minimize):
##         """
##         generate the index of next improving neigh according to sumArr only (surrogate of fitness)
##         """
##         # check improving move
##         improve = False
##         self.improveA = []
##         for i in range(self.dim):
##             #if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i] < - self.threshold):
##             if (minimize == True and self.sumArr[i] > - self.threshold) or (minimize == False and self.sumArr[i] < self.threshold):
##                 self.improveA.append(i)
##                 improve = True

##         if improve == False:
##             return False, None

##         # randomly pick an improving move, which takes only constant time
##         bestI = random.choice(self.improveA)

##         return True, bestI

##     def genFitBest2(self,minimize):
##         """
##         generate the index of best distance 2 neighborhoods according to sumArr only (surrogate of fitness)

##         return: 1) whether there is an improving move in distance 2 neigh
##                 2) the index of best distance 1 neigh for taking the next move
##                 3) the number of evaluations consumed by this step
##         """
##         improve = False
##         self.improveA = []
##         neighImprove = []

##         # checking the distance 1 neigh
##         for i in range(self.dim):
##             #if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i] < -self.threshold):
##             if (minimize == True and self.sumArr[i] > - self.threshold) or (minimize == False and self.sumArr[i] < self.threshold ):
##                 self.improveA.append(i)
##                 neighImprove.append(Struct(index =[i], val = self.sumArr[i])) # add dist 1 neigh into consideration as well
##                 improve = True

## #        print self.sumArr
## #        print self.improveA

##         # checking the distance 2 neigh, remember to preserve context
##         for i in range(self.dim):
##             self.mirrorParam() # everything is reset, pretending nothing happened
##             self.updateFake(i)
##             #self.updateWASfake(i)
##             for j in [k for k in range(self.dim) if k!=i]:
##                 self.sumArrFake[j] = self.sumArrFake[j]+self.sumArr[i]
##                 if (minimize == True and self.sumArrFake[j] > - self.threshold) or (minimize == False and self.sumArrFake[j]< self.threshold):
##                 #if (minimize == True and self.sumArrFake[j] > self.threshold) or (minimize == False and self.sumArrFake[j]< -self.threshold):
##                     neighImprove.append(Struct(index =[i,j], val = self.sumArrFake[j]))
##                     improve = True

## #        for i in range(len(neighImprove)):
## #            print neighImprove[i].index, neighImprove[i].val

## #        for i in neighImprove:
## #            print i.index, i.val

##         if improve == False:
##             #return False, None, self.dim*self.dim
##             return False, None, self.dim

##         for i in range(len(neighImprove)):
##             if i == 0:
##                 best = neighImprove[i].val
##                 bestI = neighImprove[i].index
##             # if ( best<neighImprove[i].val - self.threshold and minimize == True) or ( best>neighImprove[i].val + self.threshold and minimize == False ): # seek for max S
##             elif ( best<neighImprove[i].val + self.threshold and minimize == True) or ( best>neighImprove[i].val - self.threshold and minimize == False ): # seek for max S
##                 best = neighImprove[i].val
##                 bestI = neighImprove[i].index

##         bestIlist = []
##         for i in range(len(neighImprove)):
##             #if abs(best - neighImprove[i].val) < self.threshold:
##             if abs(best - neighImprove[i].val) < self.threshold:
##                 candI = neighImprove[i].index
##                 if candI not in bestIlist:
##                     bestIlist.append(candI)

##         #print 'bestIlist',bestIlist
##         bestI = random.choice(bestIlist)
## #        print 'bestList', bestIlist
## #        print 'bestI', bestI
##         if type(bestI) is int:
##             # make a consistent interface
##             bestI = [bestI]

##         #return True, bestI, self.dim*self.dim
##         return True, bestI, self.dim

##     def genFitBestsm(self,minimize):
##         """
##         generate the index of best neigh according to sumArr only (surrogate of fitness)
##         """
##         #TODO need to update the threshold
##         # check improving move
##         improve = False
##         self.Buffer = []
##         for i in range(self.dim):
##             if (minimize == True and self.sumArr[i] > 0) or (minimize == False and self.sumArr[i]<0):
##                 self.Buffer.append(i)
##                 improve = True

##         if improve == False:
##             return False, None

##         for i in self.Buffer:
##             if i == self.Buffer[0]:
##                 best = self.sumArr[i]
##                 bestI = i
##             elif (best<self.sumArr[i] and minimize == True) or (best>self.sumArr[i] and minimize == False): # seek for max S
##                 best = self.sumArr[i]
##                 bestI = i

##         self.P = [bestI]

##         # produce buffer list (the independent improving set)
##         for i in [ j for j in self.Buffer if j != bestI ]:
##             if i not in self.Inter:
##                 self.P.append(i)
##             else :
##                 inter = False
##                 for j in self.P:
##                     if j in self.Inter[i].arr:
##                         inter = True
##                 if inter == False:
##                     self.P.append(i)

##         return True, bestI


    ## def updateFitBest(self, p, minimize):
    ##     self.improveA.remove(p)

    ##     if p in self.Inter:
    ##         for i in self.Inter[p].arr:
    ##             # if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i] < - self.threshold ):
    ##             if (minimize == True and self.sumArr[i] > - self.threshold) or (minimize == False and self.sumArr[i]< self.threshold ):
    ##                 if i not in self.improveA:
    ##                     self.improveA.append(i)
    ##             elif i in self.improveA:
    ##                 self.improveA.remove(i)

    ##     if not self.improveA:
    ##         return False, None

    ##     random.shuffle(self.improveA)

    ##     for i in self.improveA:
    ##         if i == self.improveA[0]:
    ##             best = self.sumArr[i]
    ##             bestI = i
    ##         elif (best<self.sumArr[i] and minimize == True) or (best>self.sumArr[i] and minimize == False):
    ##             best = self.sumArr[i]
    ##             bestI = i

    ##     return True, bestI

    def steepFitDesc(self, minimize):
        # pick take the move that yields the most improvement in evaluation
        if not self.improveA:
            return False, None
        bestList = []
        
        # random.shuffle(self.improveA)

        init = False
        # find the best evaluation
        for i in self.improveA:
            if init == False:
                best = self.sumArr[i] 
                init = True
            elif (best<self.sumArr[i] + self.threshold and minimize == True) or (best > self.sumArr[i] - self.threshold and minimize == False):
                best = self.sumArr[i]
            
            # if i == self.improveA[0]:
            #     best = self.sumArr[i]
            #     # bestI = i
            # #elif (best<self.sumArr[i] - self.threshold and minimize == True) or (best>self.sumArr[i] + self.threshold and minimize == False):
            # elif (best<self.sumArr[i] + self.threshold and minimize == True) or (best > self.sumArr[i] - self.threshold and minimize == False):
            #     best = self.sumArr[i]
            #     # bestI = i

        # locate equally good best-moves
        for i in self.improveA:
            # if (abs(best - self.sumArr[i])<self.threshold):
            if (abs(best - self.sumArr[i])<self.threshold):
                bestList.append(i)
                
        # bestList.sort()
        # print 'bestList', bestList
        return True, random.choice(bestList)

    def steepFitDescPartialUpdate(self, minimize):
        if not self.improveA:
            return False, None
        bestList = []

        # random.shuffle(self.improveA)

        for i in self.improveA:
            if i == self.improveA[0]:
                best = self.sumArr[i]
                # bestI = i
            #elif (best<self.sumArr[i] - self.threshold and minimize == True) or (best>self.sumArr[i] + self.threshold and minimize == False):
            elif (best<self.sumArr[i] + self.threshold and minimize == False) or (best > self.sumArr[i] - self.threshold and minimize == True):
                best = self.sumArr[i]
                # bestI = i

        for i in self.improveA:
            if (abs(best - self.sumArr[i])<self.threshold):
                bestList.append(i)
        # bestList.sort()
        
        return True, random.choice(bestList)

    def nextFitDescPartialUpdate(self):
        if not self.improveA:
            return False, None
        
        return True, self.improveA.pop()

        
    def steepMeanDesc(self, minimize):
        if not self.improveSC:
            return False, None

        random.shuffle(self.improveSC)

        # find the best value
        for i in self.improveSC:
            if i == self.improveSC[0]:
                best = self.SC[i]
                bestI = i
            # elif ( best<self.SC[i] - self.threshold and minimize == True ) or ( best>self.SC[i] + self.threshold and minimize == False ): # seek for max S
            elif ( best<self.SC[i] + self.threshold and minimize == True ) or ( best>self.SC[i] - self.threshold and minimize == False ): # seek for max S
                best = self.SC[i]
                bestI = i
        return True, bestI

    ## def updateFitNext(self, p, minimize):
    ##     """
    ##     find the next improving move by the similar update trick
    ##     """
    ##     self.improveA.remove(p)
    ##     if p in self.Inter:
    ##         for i in self.Inter[p].arr:
    ##             #if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i]< - self.threshold ):
    ##             if (minimize == True and self.sumArr[i] > - self.threshold) or (minimize == False and self.sumArr[i]< self.threshold ):
    ##                 if i not in self.improveA:
    ##                     self.improveA.append(i)
    ##             elif i in self.improveA:
    ##                 self.improveA.remove(i)

    ##     if not self.improveA:
    ##         return False, None, evalCount

    ##     # randomly pick an improving move, which takes only constant time
    ##     bestI = random.choice(self.improveA)

    ##     return True, bestI, evalCount

    def nextDesc(self):
        """
        find the next improving move by the similar update trick
        """
        if not self.improveA:
            return False, None

        # randomly pick an improving move, which takes only constant time
        # bestI = random.choice(self.improveA)
        bestI = self.improveA.pop()
        
        # for i in range(self.dim):
        #     print self.sumArr[i]

        return True, bestI

##     def updateFitBest2(self, P, minimize):
##         """
##         generate the index of best distance 2 neighborhoods according to sumArr only (surrogate of fitness), by performing partial updates
##         """
##         neighImprove = []
##         evalCount = 0

##         for p in P:
##             if p in self.improveA:
##                 self.improveA.remove(p)
##             if p in self.Inter:
##                 for i in self.Inter[p].arr:
##                     evalCount = evalCount + 1
##                     # if (minimize == True and self.sumArr[i] > self.threshold) or (minimize == False and self.sumArr[i]<-self.threshold):
##                     if (minimize == True and self.sumArr[i] > - self.threshold) or (minimize == False and self.sumArr[i]< self.threshold):
##                         if i not in self.improveA:
##                             self.improveA.append(i)
##                     elif i in self.improveA:
##                         self.improveA.remove(i)
##         p = P[-1]

##         for i in self.improveA:
##             """ add distance 1 neigh under consideration """
##             neighImprove.append(Struct(index=[i], val = self.sumArr[i]))

##         # checking the distance 2 neigh, remember to preserve context
##         for i in [k for k in range(self.dim) if k!=p]:
##             self.mirrorParam() # everything is reset, pretending nothing happened
##             self.updateFake(i)
##             for j in [k for k in range(self.dim) if k!=i]:
##                 self.sumArrFake[j] = self.sumArrFake[j]+self.sumArr[i]
##                 # if (minimize == True and self.sumArrFake[j] > self.threshold) or (minimize == False and self.sumArrFake[j]<-self.threshold):
##                 if (minimize == True and self.sumArrFake[j] > - self.threshold) or (minimize == False and self.sumArrFake[j] < self.threshold):
##                     neighImprove.append(Struct(index =[i,j], val = self.sumArrFake[j]))

## #        for i in range(len(neighImprove)):
## #            print neighImprove[i].index, neighImprove[i].val

##         if not neighImprove:
##             #return False, None, evalCount + self.dim * (self.dim-1)
##             return False, None, evalCount

##         for i in range(len(neighImprove)):
##             if i == 0:
##                 best = neighImprove[i].val
##                 bestI = neighImprove[i].index[0]
##             # elif (best<neighImprove[i].val - self.threshold and minimize == True) or (best>neighImprove[i].val + self.threshold and minimize == False): # seek for max S
##             elif (best<neighImprove[i].val + self.threshold and minimize == True) or (best>neighImprove[i].val - self.threshold and minimize == False): # seek for max S
##                 best = neighImprove[i].val
##                 bestI = neighImprove[i].index[0]

##         bestIlist = []
##         for i in range(len(neighImprove)):
##             if abs(best - neighImprove[i].val) < self.threshold:
##                 candI = neighImprove[i].index
##                 if candI not in bestIlist:
##                     bestIlist.append(candI)

##         bestI = random.choice(bestIlist)
##         if type(bestI) is int:
##             # make a consistent interface
##             bestI = [bestI]

##         #return True, bestI, evalCount + self.dim * (self.dim-1)
##         return True, bestI, evalCount

    def genImproveSC(self,minimize):
        """
        generate the index of best neigh according to {S_p(X)-2/N \Sigma_{i=1}^{N}C_{ip}(X)} only (surrogate of fitness)
        """
        # check improving move
        self.improveSC = []
        for i in range(self.dim):
            # if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i]< - self.threshold):
            if (minimize == True and self.SC[i] > - self.threshold) or (minimize == False and self.SC[i] < self.threshold):
                self.improveSC.append(i)

    cpdef flip(self, i):
        if self.bit[i]=='0':
            self.bit[i]='1'
        else:
            self.bit[i]='0'
        # print self.bit
            
    cpdef destructorWal(self,fitName):
        """ 
        free memory to avoid memory leaks, espcially for executing multiple runs
        """
        cdef int i,j
        # free(self.bit)

        free(self.sumArr)

        for i in xrange(self.dim):
            free(self.C[i])
        free(self.C)

        for i in xrange(len(self.model.WA)):
            free(self.WAS[i].arr)
            free(self.WAS+i)
        free(self.WAS)

        # for i in xrange(self.dim):
        #     free(self.lookup[i])
        # free(self.lookup)

        if self.lookup != NULL:
            for i in xrange(self.dim):
                if self.lookup[i] != NULL:
                    free(self.lookup[i])
            free(self.lookup)

        
        # for i in xrange(len(self.model.WA)):
        #     for j in self.model.WA[i].arr:
        #         if len(self.model.WA[i].arr)>1: # for at least order Walsh terms
        #             if self.Inter[j] != NULL:
        #                 free(self.Inter[j][0].arr)
        #                 free(self.Inter[j][0].WI)
        #                 free(self.Inter[j])
        #                 self.Inter[j] = NULL
        # free(self.Inter)

        # for j in xrange(self.dim):
        #     if self.Inter[j] != NULL:
        #         free(self.Inter[j][0].arr)
        #         free(self.Inter[j])
        #         self.Inter[j] = NULL
        # free(self.Inter)
        
        for i in xrange(self.dim):
            for j in xrange(self.infectBit[i][0].size()):
                #del self.infectBit[i][0][j].arr
                free( self.infectBit[i][0][j].arr )
                free( self.infectBit[i][0][j] )
            #del self.infectBit[i]
            free( self.infectBit[i] )
        free(self.infectBit)

        if fitName == 'mean':
            """ release extra memory for performing local search using mean """
            free(self.Z)
            free(self.SC)
            for i in xrange(self.dim):
                free(self.orderC[i])
            free(self.orderC)

    def destructorWalU(self,fitName):
        """ 
        free memory to avoid memory leaks, espcially for executing multiple runs
        """
        cdef int i,j

        free(self.sumArr)

        for i in xrange(len(self.model.WA)):
            free(self.WAS[i].arr)
        free(self.WAS)
        
        for i in xrange(self.dim):

            if self.Inter[i] != NULL:
                free(self.Inter[i].arr)
                free(self.Inter[i].WI)
                free(self.Inter[i])

        if self.lookup != NULL:
            for i in xrange(self.dim):
                if self.lookup[i] != NULL:
                    free(self.lookup[i])
            free(self.lookup)


        if fitName == 'mean':
            """ release extra memory for performing local search using mean """
            free(self.Z)
            free(self.SC)
            for i in xrange(self.dim):
                free(self.orderC[i])
            free(self.orderC)


    def destructorBfUpdate(self):
        """
        destructor for Bfupdate implementation
        """
        # free(self.bit)
        free(self.sumArr)

        if self.lookup != NULL:
            for i in xrange(self.dim):
                if self.lookup[i] != NULL:
                    free(self.lookup[i])
            free(self.lookup)
        
        
    ## def genMeanBest(self,minimize):
    ##     """
    ##     generate the index of best neigh according to {S_p(X)-2/N \Sigma_{i=1}^{N}C_{ip}(X)} only (surrogate of fitness)
    ##     """
    ##     # check improving move
    ##     improve = False
    ##     self.improveA = []
    ##     for i in range(self.dim):
    ##         # if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i]< - self.threshold):
    ##         if (minimize == True and self.SC[i] > - self.threshold) or (minimize == False and self.SC[i] < self.threshold):
    ##             self.improveA.append(i)
    ##             improve = True

    ##     if improve == False:
    ##         return False, None

    ##     random.shuffle(self.improveA)

    ##     # find the best value
    ##     for i in self.improveA:
    ##         if i == self.improveA[0]:
    ##             best = self.SC[i]
    ##             bestI = i
    ##         # elif ( best<self.SC[i] - self.threshold and minimize == True ) or ( best>self.SC[i] + self.threshold and minimize == False ): # seek for max S
    ##         elif ( best<self.SC[i] + self.threshold and minimize == True ) or ( best>self.SC[i] - self.threshold and minimize == False ): # seek for max S
    ##             best = self.SC[i]
    ##             bestI = i
    ##     return True, bestI

    ## def genMeanNext(self,minimize):
    ##     """
    ##     generate the index of next improving neigh according to {S_p(X)-2/N \Sigma_{i=1}^{N}C_{ip}(X)} only (surrogate of fitness)
    ##     """
    ##     # check improving move
    ##     improve = False
    ##     self.improveA = []
    ##     for i in range(self.dim):
    ##         # if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i] < - self.threshold):
    ##         if (minimize == True and self.SC[i] > - self.threshold) or (minimize == False and self.SC[i] < self.threshold):
    ##             self.improveA.append(i)
    ##             improve = True

    ##     if improve == False:
    ##         return False, None

    ##     bestI = random.choice(self.improveA)
    ##     return True, bestI

    ## def updateMeanBest(self, p, minimize):
    ##     self.improveA.remove(p)
    ##     if p in self.Inter:
    ##         for i in self.Inter[p].arr:
    ##             # if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i]< - self.threshold):
    ##             if (minimize == True and self.SC[i] > - self.threshold) or (minimize == False and self.SC[i] < self.threshold):
    ##                 if i not in self.improveA:
    ##                     self.improveA.append(i)
    ##             elif i in self.improveA:
    ##                 self.improveA.remove(i)

    ##     if not self.improveA:
    ##         return False, None

    ##     random.shuffle(self.improveA)

    ##     for i in self.improveA:
    ##         if i == self.improveA[0]:
    ##             best = self.SC[i]
    ##             bestI = i
    ##         # elif (best<self.SC[i] - self.threshold and minimize == True) or (best>self.SC[i]+self.threshold and minimize == False): # seek for max S
    ##         elif (best<self.SC[i] + self.threshold and minimize == True) or (best>self.SC[i] - self.threshold and minimize == False): # seek for max S
    ##             best = self.SC[i]
    ##             bestI = i

    ##     return True, bestI

    ## def updateMeanNext(self, p, minimize):
    ##     self.improveA.remove(p)
    ##     if p in self.Inter:
    ##         for i in self.Inter[p].arr:
    ##             # if (minimize == True and self.SC[i] > self.threshold) or (minimize == False and self.SC[i]< - self.threshold):
    ##             if (minimize == True and self.SC[i] > - self.threshold) or (minimize == False and self.SC[i]< self.threshold):
    ##                 if i not in self.improveA:
    ##                     self.improveA.append(i)
    ##             elif i in self.improveA:
    ##                 self.improveA.remove(i)

    ##     if not self.improveA:
    ##         return False, None

    ##     bestI = random.choice(self.improveA)
    ##     return True, bestI

    ## def binCount(self, arr, bit):
    ##     """
    ##     count the one bit of union self.model.WA[i].arr and bit
    ##     """
    ##     s = 0
    ##     for i in arr:
    ##         if bit[i] == '1':
    ##             s = s + 1
    ##     return s

    def initIndivPy(self):
        """ initial the search inidividual with random bit string """
        # self.bit = <char*>malloc( (self.dim+1) * sizeof(char) )
        # self.bit = new char[self.dim+1]()
        self.bit=[]
        for j in xrange(self.dim):
            if random.random()<0.5:
                self.bit.append('0')
            else:
                # self.bit[j] = '1'
                self.bit.append('1')

    # cdef initIndiv(self):
    #     """ initial the search inidividual with random bit string """
    #     self.bit = <char*>malloc( (self.dim+1) * sizeof(char) )
    #     # self.bit = new char[self.dim+1]()
    #     for j in xrange(self.dim):
    #         if random.random()<0.5:
    #             self.bit[j] = '0'
    #         else:
    #             self.bit[j] = '1'
    #     self.bit[self.dim] = '\0'
    #     self.bit[self.dim-1] = 'a'

    cpdef updateEval(self, I):
        """
        update the evaluation function according to sumArr
        """
        self.fit = self.fit - 2*self.sumArr[I]

    cpdef updateEvalPartialUpdate(self, I):
        self.fit = self.fit + self.sumArr[I]
        
    cpdef updateSumArr(self, int q, list old):
        """
        for the partial update implementation:
        update the first derivative according to second derivative

        only consider the p setting that possibly interacts with q
        """
        cdef int i, p
        
        if self.model.Inter[q]!=None:
            for p in self.model.Inter[q]:
                # print 'p', p
                # calculate the sum of second derivative
                # manipulate on the length-k extracted string
                for i in self.model.getU(p,q):
                    # print 'update subfunc i', i
                    # print 'before', old.bit
                    self.sumArr[p] = self.sumArr[p] + self.sumTerm(old, i, p, q)
                    
        self.sumArr[q] = - self.sumArr[q]
        
        
    cdef float sumTerm(self, list bitStr, int i, int p, int q):
        """
        compute the inner summation term, flip p and q bits for i
        subfunction
        
        f_i(x^{qp}) - f_i(x^{q}) - f_i(x^{p}) + f_i(x)
        """
        cdef float s, saveDiffp,
        cdef int j, realI, pI, qI
        cdef list bits

        """ extract corresponding bits """
        bits = []
        for j in xrange(self.model.k+1):
            realI = self.model.neighs[i][j]
            bits.append(bitStr[realI])
            if  realI == p:
                pI = j
            elif realI == q:
                qI = j
        I = int(''.join(bits),2)
        s = self.model.getFuncVal(i, I)

        # print 'p'
        saveDiffp = (1-2*int(bitStr[p])) * pow(2,self.model.k-pI)
        I = int(I + saveDiffp)
        s = s - self.model.getFuncVal(i,I)   
        
        # pq, flip q
        I = int(I + (1-2*int(bitStr[q])) * pow(2,self.model.k-qI))
        s = s + self.model.getFuncVal(i,I) 
        
        # q
        I = int(I - saveDiffp)
        s = s - self.model.getFuncVal(i,I) 
        
        return s

                
    cpdef printSumArr(self):
        for i in range(self.dim):
            print i, self.sumArr[i]
        print

    cpdef compFitG(self):
        self.fitG = self.fit - 2/ float(self.dim) * (sumC(self.sumArr, self.dim))

        
    # cdef ComArr* genComb(self,int N) nogil:
    cdef ComArr* genComb(self,int N):
        """
        Generate C_N^2 sequence, index are stored, because they are more general, Implemented in an *incremental* fashion.
        """
        cdef int c, j, i, counter
        cdef ComArr* ptr

        if self.lookup[N] != NULL: # the key exists before
            return self.lookup[N]
        else : # the key not found
            
            # c = biomial(N, 2)
            c = N*(N-1)/2
            counter = 0

            ptr = <ComArr*> malloc(sizeof(ComArr))
            ptr.arr = <int**> malloc(sizeof(int*)*c)
            ptr.size = c

            for i in xrange(c):
                ptr.arr[i] = <int*> malloc(sizeof(int)*2)

            for i in xrange(N):
                for j in xrange(i+1,N):
                    ptr.arr[counter][0] = i
                    ptr.arr[counter][1] = j
                    counter = counter + 1
            self.lookup[N] = ptr
            return ptr



cdef int binCount(list arr, list bit):
    """
    count the one bit of union self.model.WA[i].arr and bit
    """
    cdef int s = 0
    cdef int i
    for i in arr:
        if bit[i] == '1':
            s = s + 1
    return s

cdef float sumC(float * a, int d):
    cdef float s = 0

    for i in xrange(d):
        s = s+a[i]

    return s

# @cython.cdivision(True)
# cdef int biomial(int N, int K) nogil:
#     """ compute the combination of N choose K """
#     return factorial(N)/( factorial(K) * factorial(N-K) )

# cdef int factorial(int N) nogil:
#     """ compute N! """
#     cdef int c, fact = 1
#     for c in xrange(1,N+1):
#         fact = fact*c
#     return fact
