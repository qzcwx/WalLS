# spin glasses class
import sys
import numpy as np
import math
class Struct:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


cdef class SpinGlass:
    cdef public int n              # number of variables
    cdef public list dim
    cdef public int m              # number of subfunctions
    cdef public int k              # number of variables per subfunction
    cdef public list WA
    
    def __init__(self, inN, inSize):
        self.n = inN
        self.dim = inSize
        # print self.dim
        self.readFile()
        
    cdef readFile(self):
        # read walsh terms directly from files
        self.WA = []
        # print 'benchmark/spin/spin2g-%dX%d-I0' %(self.dim[0],self.dim[1])
        fName = 'benchmark/spin/spin2g-%dX%d-I0' %(self.dim[0],self.dim[1])
        with open(fName, 'r') as f:
            fl = f.readline()
            self.n = int(fl.split()[0])
            self.m = int(fl.split()[1])

        v = np.genfromtxt(fName, skip_header=1, usecols=(0,1),dtype='int')
        c = np.genfromtxt(fName, skip_header=1, usecols=(2))
        # print self.n, self.m
        # print v
        # print c
        for i in xrange(len(c)):
            if c[i] != 0:
                t = list(v[i]-1)
                t.sort()
                self.WA.append(Struct(arr=t, w = c[i]))
                # print self.WA[i].arr
                # print self.WA[i].w
        
        # print self.WA

    cpdef double compFit(self, bitStr):
        """ compute fitness  """
        """ 1 -> -1; 0 -> +1 """
        cdef int i
        cdef double s = 0.0
        cdef int minusCount
        # print 'len',len(self.WA)
        # print 'str len', len(bitStr)
        for i in xrange(len(self.WA)):
            # print 'i',i
            print self.WA[i].arr
            minusCount = 0
            for j in self.WA[i].arr:
                # print 'j',j
                if bitStr[j] == '1':
                    minusCount = minusCount + 1
            print math.pow(-1,minusCount) * self.WA[i].w 
            s = s + math.pow(-1,minusCount) * self.WA[i].w 
        return -s
