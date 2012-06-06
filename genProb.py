"""
Generates 10 instances for each of the following configuraion (10 combinations)
        NK Landcapes: K = 0, 2, 4, 8, 16; N = 20, 50, 100. (5*3=15 combinations)
        NK-Q Landscapes: K = 0, 2, 4, 8, 16; N = 20, 50, 100; q = 2, 4, 8, 16. (5*3*4=60 combinations)
"""
import nkLandscape as nk
import nkqLandscape as nkq
import random as random


iRange = range(10)
rseed = 0

#nRange = [4]
nRange = [15]
kRange = [2]
cRange = [15, 65, 120]

#cRange = [4,17,32]
tRange = [4,5,6,7]
qRange = [2]


prob = 'NKQ'

prefixNK = './benchmark/NK/'
prefixNKQ = './benchmark/NKQ/'

for n in nRange:
    for k in kRange:
        for c in cRange:
            for t in tRange:
                for i in iRange:
                    random.seed(rseed+i)
                    model = nk.NKLandscape(n,k,c)
                    if prob == 'NK':
                        # NK landscapes:
                        model.exportToFile(prefixNK+'NK-N'+str(n)+'-K'+str(k)+'-C'+str(c)+'-I'+str(i))
                    elif prob == 'NKQ':
                        # NKq landscapes:
                        for q in qRange:
                            model = nkq.NKQLandcape(n, k, c, q, t)
                            model.exportToFile(prefixNKQ+'NKQ-N'+str(n)+'-K'+str(k)+'-C'+str(c)+'-I'+str(i)+'-Q'+str(q)+'-T'+str(t))
