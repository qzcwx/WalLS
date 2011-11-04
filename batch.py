import os

q = 2
for n in [20]:
    for k in [0, 2, 4, 8, 16]:
        if k+1<=n:
            for a in ['SATGA', 'SATLS', 'GA', 'LS']:
                for p in ['NKQ','NK']:
                    print 'qsub sat.pbs -v P='+p+',A='+a+',N='+str(n)+',K='+str(k)+',Q='+str(q)
                    os.system( 'qsub sat.pbs -v P='+p+',A='+a+',N='+str(n)+',K='+str(k)+',Q='+str(q) )
