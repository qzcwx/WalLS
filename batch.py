import os

for n in [100, 50, 2]:
    for k in [0, 2, 4, 8, 16]:
        if k+1<=n:
            for a in ['SATGA', 'SATLS', 'GA', 'LS']:
                print 'algo', a, 'n', n, 'k', k 
                os.system('qsub sat.pbs -v A='+a+',N='+str(n)+',K='+str(k) )
