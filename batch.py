import os

for n in [20, 50, 100] :
    for k in [0, 2, 4, 8, 16]:
        for a in ['SATGA', 'GA']:
            for q in [2, 4, 8, 16]:
                for p in ['NKQ','NK']:
                    print 'qsub sat.pbs -v P='+p+',A='+a+',N='+str(n)+',K='+str(k)+',Q='+str(q)
                    os.system( 'qsub sat.pbs -v P='+p+',A='+a+',N='+str(n)+',K='+str(k)+',Q='+str(q) )
