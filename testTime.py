"""
Test the runtime for rLS with mean of neighborhoods as fitness, 
the time for computing walsh coefficients is stored as well.
"""
import os

def resultExist(probName,algoName,fitName,inst,s,c,n,k,q):
    """ check whether the results have been produced or not """
    if probName == 'NKQ':
        nameOfF = './result/'+probName+'-'+algoName+'-F'+fitName+'-C'+c+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'-K'+str(k)+'-Q'+str(q)+'.txt'
    elif probName == 'NK':
        nameOfF = './result/'+probName+'-'+algoName+'-F'+fitName+'-C'+c+'-I'+str(inst)+'-S'+str(s)+'-N'+str(n)+'-K'+str(k)+'.txt'

    if os.path.isfile(nameOfF)==True:
        print nameOfF, 'exists!!!'
    return os.path.isfile(nameOfF)

nRange = [20]
kRange = [2]
#kRange = [2,4,8]
#nRange = [20,50,100]
iRange = [0]
aRange = ['rLS']
fRange = ['mean']
#cRange = ['supm', 'wal']
cRange = ['wal']
pRange = ['NK']

for k in kRange:
    for n in nRange :
        for i in iRange :
            for a in aRange :
                for f in fRange :
                    for c in cRange :
                        for p in pRange :
                            if p == 'NKQ':
                                for q in [2]:
                                    s = 1
                                    print 'python demo.py '+c+' '+p+' '+a+' '+f+' '+str(i)+' '+str(s)+' '+str(n)+' '+str(k)+' '+str(q)
                                    if resultExist(p,a,f,i,s,c,n,k,q) == False:
                                        os.system('python demo.py '+c+' '+p+' '+a+' '+f+' '+str(i)+' '+str(s)+' '+str(n)+' '+str(k)+' '+str(q))
                            elif p == 'NK': # for NK problem
                                q = 0
                                s = 1
                                print 'python demo.py '+c+' '+p+' '+a+' '+f+' '+str(i)+' '+str(s)+' '+str(n)+' '+str(k)+' '+str(q)
                                if resultExist(p,a,f,i,s,c,n,k,q) == False:
                                    os.system('python demo.py '+c+' '+p+' '+a+' '+f+' '+str(i)+' '+str(s)+' '+str(n)+' '+str(k)+' '+str(q) )
