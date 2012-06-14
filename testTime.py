"""
Test the runtime for rLS with mean of neighborhoods as fitness, 
the time for computing walsh coefficients is stored as well.
"""
import os

def resultExist(probName,algoName,fitName,inst,s,w,m,n,k,q,t,c):
    """ check whether the results have been produced or not """
    nameOfDir = './result/'
    if probName == 'NKQ':
        nameOfF = nameOfDir+probName+'-'+algoName+'-F'+fitName+'-M'+m+'-I'+str(inst)+'-S'+str(s)+'-W'+str(w)+'-N'+str(n)+'-K'+str(k)+'-C'+str(c)+'-Q'+str(q)+'-T'+str(t)+'.txt'
    elif opt.probName == 'NK':
        nameOfF = nameOfDir+probName+'-'+algoName+'-F'+fitName+'-M'+m+'-I'+str(inst)+'-S'+str(s)+'-W'+str(w)+'-N'+str(n)+'-C'+str(c)+'-K'+str(k)+'.txt'

    # if probName == 'NKQ':
    #     nameOfF = './result/'+probName+'-'+algoName+'-F'+fitName+'-C'+c+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'-K'+str(k)+'-Q'+str(q)+'.txt'
    # elif probName == 'NK':
    #     nameOfF = './result/'+probName+'-'+algoName+'-F'+fitName+'-C'+c+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'-K'+str(k)+'.txt'

    if os.path.isfile(nameOfF)==True:
        print nameOfF, 'exists!!!'
    return os.path.isfile(nameOfF)

nRange = [15]
kRange = [2]
tRange = [4,5,6,7]
cRange = [15, 65, 120]
iRange = range(1)
w = 0
aRange = ['rLS']
fRange = ['fit','mean']
mRange = ['walWalk']
pRange = ['NK','NKQ']


# python run.py -m walWalk -p NKQ -a rLS -f fit -i 0 -s 1 -n 4 -k 2 -q 2 -t 7 -c 4

for k in kRange:
    for n in nRange :
        for i in iRange :
            for a in aRange :
                for m in mRange:
                    for f in fRange :
                        for t in tRange:
                            for c in cRange :
                                for p in pRange :
                                    if p == 'NKQ':
                                        for q in [2]:
                                            s = 1
                                            if resultExist(p,a,f,i,s,w,m,n,k,q,t,c) == False:
                                                os.system('python run.py -m'+m+' -p'+p+' -a'+a+' -f'+f+' -i'+str(i)+' -s'+str(s)+' -n'+str(n)+' -k'+str(k)+' -q'+str(q)+' -t'+str(t)+' -c'+str(c))
                                            elif p == 'NK': # for NK problem
                                                q = 0
                                                s = 1
                                                if resultExist(p,a,f,i,s,w,m,n,k,q,t,c) == False:
                                                    os.system('python run.py -m'+m+' -p'+p+' -a'+a+' -f'+f+' -i'+str(i)+' -s'+str(s)+' -n'+str(n)+' -k'+str(k)+' -q'+str(q) + ' -c'+str(c) )
                                                
