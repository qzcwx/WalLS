import matplotlib.pyplot as plt
import numpy as np

nameOfDir = './result/'
pRange = ['NK','NKQ']
algoName = 'rLS'
fitName = 'fit'
compMeth = 'walWalk'
#compMeth = 'walRest'
cRange = ['walWalk','walRest']
compName = ['WalLS-Walk','WalLS-Rest']
inst = 0
color = ['red','green']
s = 30
n = 100
k = 4
q = 2
i = 0

plt.ion()
for probName in pRange:
    plt.subplot(1,2,pRange.index(probName)+1)
    for compMeth in cRange:
        if probName == 'NKQ':
            nameOfF = nameOfDir+'Trace-'+probName+'-'+algoName+'-F'+fitName+'-C'+compMeth+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'-K'+str(k)+'-Q'+str(q)+'-R'+str(i)+'.txt'
        elif probName == 'NK':
            nameOfF = nameOfDir+'Trace-'+probName+'-'+algoName+'-F'+fitName+'-C'+compMeth+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'-K'+str(k)+'-R'+str(i)+'.txt'

        trace = np.genfromtxt(nameOfF)
        plt.plot(trace[:,0],trace[:,1],'.',label=compName[cRange.index(compMeth)],c=color[cRange.index(compMeth)])
        plt.legend()
        if probName == 'NKQ':
            plt.title("NKq-Landscape (q=2)")
            plt.xlabel("The indices of sampling points")
            plt.ylabel("Function Value")
            #plt.title("Sampling points on NKq-Landscape (q=2), N=100, K=4")
        elif probName == 'NK':
            plt.title("NK-Landscape")
            plt.xlabel("The indices of sampling points")
            plt.ylabel("Function Value")
