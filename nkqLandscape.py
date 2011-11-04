import nkLandscape as nk
import random
import numpy as np
import math
import pdb

class NKQLandcape(nk.NKLandscape):
    def __init__(self,inN,inK,inQ):
        nk.NKLandscape.__init__(self, inN, inK)
        self.q = inQ
        self.genFuncQ()

    def genFuncQ(self):
        self.func = []
        for i in range(self.n):
            oneFunc = []
            for j in range(int(math.pow(2,self.k+1))):
                oneFunc.append(random.randint(0,self.q))
            self.func.append(oneFunc)
