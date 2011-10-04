# Compute Auto-Correlation with respect to Walsh coefficients w_i
import math

def autoCorr(s, w):
    allSolNum = len(w)
    n = int(math.log(allSolNum, 2))
    print 'numOfSol', allSolNum, 'n', n, 'rangeI', [j for j in range(allSolNum) if j != 0]
    sum1 = 0
    sum2 = 0
    for i in [j for j in range(allSolNum) if j != 0]:
        print 'i', i
        squaredW = w[i] * w[i]
        print 'squaredW', squaredW
        sum1 = sum1 + math.pow(lamd(i, n),s) * squaredW
        sum2 = sum2 + squaredW
    return sum1/sum2

def lamd(i, n):
    oneBit = bin(i)
    oneBit = oneBit[2:]
    numOfOne = oneBit.count('1')
    print 'i', i, 'str', oneBit, 'numOfOne', numOfOne
    print 'lamd', (1 -  2.0*numOfOne/n)
    return (1 -  2.0*numOfOne/n)
