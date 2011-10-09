# Compute Auto-Correlation with respect to Walsh coefficients w_i
import math

def autoCorr(s, w):
    allSolNum = len(w)
    n = int(math.log(allSolNum, 2))
   # print 'numOfSol', allSolNum, 'n', n, 'rangeI', [j for j in range(allSolNum) if j != 0]
    sum1 = 0
    sum2 = 0
    for i in [j for j in range(allSolNum) if j != 0]:
    #    print 'i', i
        squaredW = w[i] * w[i]
    #    print 'squaredW', squaredW
        sum1 = sum1 + math.pow(lamd(i, n),s) * squaredW
        sum2 = sum2 + squaredW
    return sum1/sum2

def autoCorrLinerSpace(s, w):
    sum1 = 0
    sum2 = 0
    n = len(list(w.viewkeys())[0])
    for bitstr in [ i for i in w.viewkeys() if i != n*'0']:
        squaredW = w[bitstr] * w[bitstr]
        sum1 = sum1 + math.pow(lamdLinear(bitstr, n), s) * squaredW
        sum2 = sum2 + squaredW
    return sum1/sum2

def lamdLinear(i, n):
    numOfOne = i.count('1')
    return (1 -  2.0*numOfOne/n)

def lamd(i, n):
    oneBit = bin(i)
    oneBit = oneBit[2:]
    numOfOne = oneBit.count('1')
    return (1 -  2.0*numOfOne/n)
