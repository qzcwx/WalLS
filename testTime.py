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

for k in [2, 4, 8]:
    for n in [200, 500]:
        for i in [0]:
            for a in ['rLS']:
                for f in ['fit','mean']:
                    for c in ['wal']:
                        for p in ['NK','NKQ']:
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
                                    os.system('python demo.py '+c+' '+p+' '+a+' '+f+' '+str(i)+' '+str(s)+' '+str(n)+' '+str(k)+' '+str(q) + ' >' + p+'-'+a+'-F'+f+'-C'+c+'-I'+str(i)+'-S'+str(s)+'-N'+str(n)+'-K'+str(k)+'-Q'+str(q)+'&' )
