import nkLandscape as nk
import tool as tl
import random
import numpy as np
import math
import pdb

class NKQLandcape(nk.NKLandscape):
    def __init__(self, inN, inK, inQ, fileName = None):
        self.q = inQ
        self.n = inN
        self.k = inK

#        self.genNeigh()
#        self.genFunc()
#        self.genFuncQ()
#        self.exportToFile(fileName)

        nk.NKLandscape.__init__(self, inN, inK, fileName)
        if fileName == None:
            self.genFuncQ()
        else:
            print 'read', fileName
            self.readFile(fileName)
        self.Kbits = tl.genSeqBits(self.k+1)

    def genFuncQ(self):
        self.func = []
        for i in range(self.n):
            oneFunc = []
            for j in range(int(math.pow(2,self.k+1))):
                oneFunc.append(np.random.randint(self.q))
            self.func.append(oneFunc)
