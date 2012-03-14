"""
Generates 10 instances for each of the following configuraion (10 combinations)
        NK Landcapes: K = 0, 2, 4, 8, 16; N = 20, 50, 100. (5*3=15 combinations)
        NK-Q Landscapes: K = 0, 2, 4, 8, 16; N = 20, 50, 100; q = 2, 4, 8, 16. (5*3*4=60 combinations)
"""
import nkLandscape as nk
import nkqLandscape as nkq
import pdb
import random as random

numOfInstances = 10
prefixNK = './benchmark/NK/'
prefixNKQ = './benchmark/NKQ/'

for n in [1000]:
    for k in [2,4,8]:
        for i in range(numOfInstances):
            # NK landscapes:
            model = nk.NKLandscape(n,k)
            model.exportToFile(prefixNK+'NK-N'+str(n)+'-K'+str(k)+'-I'+str(i))
        for q in [2]:
            for i in range(numOfInstances):
                model = nkq.NKQLandcape(n, k, q)
                model.exportToFile(prefixNKQ+'NKQ-N'+str(n)+'-K'+str(k)+'-I'+str(i)+'-Q'+str(q))
