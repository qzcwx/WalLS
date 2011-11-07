import os

i = 0
for n in [20, 50, 100] :
    for k in [0, 2, 4, 8, 16]:
        for a in ['SATGA', 'GA', 'LS', 'SATLS']:
            for p in ['NKQ','NK']:
                if p == 'NKQ':
                    for q in [2, 4, 8, 16]:
                        i = i+1
                        print 'qsub sat.pbs -v P='+p+',A='+a+',N='+str(n)+',K='+str(k)+',Q='+str(q)
                        os.system( 'qsub sat.pbs -v P='+p+',A='+a+',N='+str(n)+',K='+str(k)+',Q='+str(q) )
                else:
                    q = 0
                    i = i+1
                    print 'qsub sat.pbs -v P='+p+',A='+a+',N='+str(n)+',K='+str(k)+',Q='+str(q)
                    os.system( 'qsub sat.pbs -v P='+p+',A='+a+',N='+str(n)+',K='+str(k)+',Q='+str(q) )

print i, 'jobs submitted'
