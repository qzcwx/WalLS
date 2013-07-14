import tool as tl
import random
import numpy as np
import math
import pdb
import cython
import nkLandscape as nk

class NKQLandscape(nk.NKLandscape):

    # cdef public int n                            # number of variables
    # cdef public int k                            
    # cdef public int c                            # number of clauses
    # cdef public int q
    # cdef public int t
    
    def __init__(self, inN, inK, inC, inQ, inT=-1, fileName = None):
        self.n = inN
        self.k = inK
        self.c = inC
        self.q = inQ
        self.t = inT

#        self.genNeigh()
#        self.genFunc()
#        self.genFuncQ()
#        self.exportToFile(fileName)

        nk.NKLandscape.__init__(self, inN, inK, inC, fileName)
        if fileName == None:
            if self.t == -1:
                self.genFuncQ()
            else:
                self.genFuncQT()
        else:
            self.readFile(fileName)
        self.Kbits = tl.genSeqBits(self.k+1)
        
    def genFuncQ(self):
        self.f = []
        for i in range(self.c):
            oneFunc = []
            for j in range(int(math.pow(2,self.k+1))):
                oneFunc.append(np.random.randint(self.q))
            self.f.append(oneFunc)

    def genFuncQT(self):
        self.f = []
        for i in range(self.c):
            oneFunc = []
            tAss = random.sample(range(int(math.pow(2,self.k+1))),self.t)
            for j in range(int(math.pow(2,self.k+1))):
                if j in tAss:
                    oneFunc.append(1)
                else:
                    oneFunc.append(0)
            self.f.append(oneFunc)


class NonNKQLandscape(NKQLandscape):
    """ Non-uniform random NKQ-landscapes """
    def __init__(self, inN, inK, inC, inQ, inT=-1, fileName = None):
        self.n = inN
        self.k = inK
        self.c = inC
        self.q = inQ
        self.t = inT
        NKQLandscape.__init__(self, inN, inK, inC, inQ, inT, fileName)
        if fileName == None:
            # regenerate neighs
            NKQLandscape.genNonNeigh()
            # if self.t == -1:
            #     self.genFuncQ()
            # else:
            #     self.genFuncQT()
        else:
            self.readFile(fileName)

