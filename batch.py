import os

#[NameOfProblem] [NameOfAlgorithm] [fit/mean/std] [I] [N] [K] [Q]

i = 0
for n in [20, 50, 100] :
    for k in [0, 2, 4, 8, 16]:
        for a in ['LS', 'GA', 'CHC']:
            if a == 'LS':
                s = 1
                for f in ['fit', 'mean', 'std']:
                    for p in ['NKQ','NK']:
                        if p == 'NKQ':
                            for q in [2, 4, 8, 16]:
                                i = i+1
                                print 'qsub sat.pbs -v P='+p+',A='+a+',F='+f+',I=1,S='+str(s)+',N='+str(n)+',K='+str(k)+',Q='+str(q)
                                os.system( 'qsub sat.pbs -v P='+p+',A='+a+',F='+f+',I=1,S='+str(s)+',N='+str(n)+',K='+str(k)+',Q='+str(q) )
                        else:
                            q = 0
                            i = i+1
                            print 'qsub sat.pbs -v P='+p+',A='+a+',F='+f+',I=1,S='+str(s)+',N='+str(n)+',K='+str(k)+',Q='+str(q)
                            os.system( 'qsub sat.pbs -v P='+p+',A='+a+',F='+f+',I=1,S='+str(s)+',N='+str(n)+',K='+str(k)+',Q='+str(q) )
            else:
                for s in [30, 100]:
                    for f in ['fit', 'mean', 'std']:
                        for p in ['NKQ','NK']:
                            if p == 'NKQ':
                                for q in [2, 4, 8, 16]:
                                    i = i+1
                                    print 'qsub sat.pbs -v P='+p+',A='+a+',F='+f+',I=1,S='+str(s)+',N='+str(n)+',K='+str(k)+',Q='+str(q)
                                    os.system( 'qsub sat.pbs -v P='+p+',A='+a+',F='+f+',I=1,S='+str(s)+',N='+str(n)+',K='+str(k)+',Q='+str(q) )
                            else:
                                q = 0
                                i = i+1
                                print 'qsub sat.pbs -v P='+p+',A='+a+',F='+f+',I=1,S='+str(s)+',N='+str(n)+',K='+str(k)+',Q='+str(q)
                                os.system( 'qsub sat.pbs -v P='+p+',A='+a+',F='+f+',I=1,S='+str(s)+',N='+str(n)+',K='+str(k)+',Q='+str(q) )

for a in ['LS', 'GA', 'CHC']:
    if a == 'LS':
        s =1
        for f in ['fit', 'mean', 'std']:
            i = i+1
            print 'qsub sat.pbs -v P='+p+',A='+a+',F='+f+',I=1,S='+str(s)+',N='+str(n)+',K='+str(k)+',Q='+str(q)
            os.system( 'qsub sat.pbs -v P='+p+',A='+a+',F='+f+',I=1,S='+str(s)+',N='+str(n)+',K='+str(k)+',Q='+str(q) )
    else:
        for s in [30,100]:
            for f in ['fit', 'mean', 'std']:
                i = i+1
                print 'qsub sat.pbs -v P='+p+',A='+a+',F='+f+',I=1,S='+str(s)+',N='+str(n)+',K='+str(k)+',Q='+str(q)
                os.system( 'qsub sat.pbs -v P='+p+',A='+a+',F='+f+',I=1,S='+str(s)+',N='+str(n)+',K='+str(k)+',Q='+str(q) )

print i, 'jobs submitted'
