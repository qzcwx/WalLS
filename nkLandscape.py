# generate NK-landscapes instances
## the encoding is in right to left fashion
## 00-0, 01-1, 10-2, 11-3 import WalshAnalysis as wal import random import numpy as np import math

import random
import WalshAnalysis as wal
import numpy as np
import math
import pdb

class NKLandscape:
    """ NK-landscape class """
    def __init__(self,inN,inK, fileName = None):
        self.n = inN
        self.k = inK

#        # for run experiments
#        if fileName == None:
#            self.genNeigh()
#            self.genFunc()
#        else:
#            self.readFile(fileName)

        # for generating benchmarks
        self.genNeigh()
        self.genFunc()
        self.exportToFile(fileName)

        self.Kbits = genSeqBits(self.k+1)

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
        self.neighs = np.genfromtxt(fName, delimiter="\t", skip_footer=self.n, autostrip=True, usecols = range(self.k)).tolist()
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

def genSeqBits(n):
    bitStr = []
    for i in range(int(math.pow(2,n))):
       bit = bin(i)
       bit = bit[2:]
       if len(bit) < n:
           bit = (n - len(bit))*'0' + bit
       bitStr.append(bit)
    return bitStr
