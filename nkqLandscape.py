import nkLandscape as nk
import random
import numpy as np
import math
import pdb

class NKQLandcape(nk.NKLandscape):
    def __init__(self, inN, inK, inQ, fileName = None):
        nk.NKLandscape.__init__(self, inN, inK, fileName)
        self.q = inQ
        if fileName == None:
            self.genFuncQ()
        #self.exportToFile(fileName)
        else:
            self.readFile(fileName)

    def genFuncQ(self):
        self.func = []
        for i in range(self.n):
            oneFunc = []
            for j in range(int(math.pow(2,self.k+1))):
                oneFunc.append(np.random.randint(self.q))
            self.func.append(oneFunc)
