# this function counts the number of local optimum over the entire search space
import random
import numpy as np
import math

def localOpt(bitStr, f):
    n = len(bitStr[0])
    num = 0
    compCount = 0
    loMark = len(bitStr)*[True]
    loMark = np.array(loMark)

    """ generate the random solution under consideration for local optimum """
    while True:
        indice = np.where(loMark==True)[0]
        if len(indice) == 0:
            break

        """ consider all the neighborhoods """
        solIndex = random.choice(indice)
        sol =  bin(solIndex)[2:]
        if len(sol) < n:
            sol = (n - len(sol))*'0' + sol 

        for i in range(n):
            neigh = sol[0:i] + flipBit(sol[i]) + sol[i+1:]
            neighIndex = int(neigh, 2)
            compCount = compCount + 1
            if f[solIndex]>f[neighIndex]:
                loMark[solIndex] = False
                break
            elif f[solIndex]<f[neighIndex]:
                loMark[neighIndex] = False
        if loMark[solIndex]==True:
            num = num + 1
            loMark[solIndex] = False
 #       print loMark
    return num

def flipBit(bit):
    if bit == '0':
        return '1'
    else:
        return '0'
