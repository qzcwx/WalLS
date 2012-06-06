import nkLandscape as nk
import tool as tl
import random
import numpy as np
import math
import pdb

class NKQLandcape(nk.NKLandscape):
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
        self.func = []
        for i in range(self.c):
            oneFunc = []
            for j in range(int(math.pow(2,self.k+1))):
                oneFunc.append(np.random.randint(self.q))
            self.func.append(oneFunc)

    def genFuncQT(self):
        self.func = []
        for i in range(self.c):
            oneFunc = []
            tAss = random.sample(range(int(math.pow(2,self.k+1))),self.t)
            for j in range(int(math.pow(2,self.k+1))):
                if j in tAss:
                    oneFunc.append(1)
                else:
                    oneFunc.append(0)
            self.func.append(oneFunc)
