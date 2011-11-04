import numpy as np
from StringIO import StringIO

class MAXSAT:
    def __init__(self):
        pass

    def setInstance(self, i):
        self.readTestsuite(i)

    def readTestsuite(self,i):
        file = open('./benchmark/uf100-0'+str(i)+'.cnf', 'r')
        dataStr = file.readlines()
        dataStr = dataStr[8:-3]
        dataStr = StringIO("".join(dataStr))
        self.data =  np.genfromtxt(dataStr, usecols = range(3))

    def compFit(self, s):
        sum = 0
        for i in range(self.data.shape[0]):
            neg = self.data[i] > 0
            sum += np.logical_not( np.logical_xor( int(s[abs(int(self.data[i,0]))-1]), neg[0])  ) or np.logical_not( np.logical_xor( int(s[abs(int(self.data[i,1]))-1]), neg[1])  ) or np.logical_not( np.logical_xor( int(s[abs(int(self.data[i,2]))-1]), neg[2])  )
        return - sum

#maxsat = MAXSAT()
#maxsat.setInstance(1)
#print maxsat.compFit(100*'1')
