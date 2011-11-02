import os

for n in [100, 50, 2]:
    for k in [0, 2, 4, 8, 16]:
        if k+1<=n:
            print 'n', n, 'k', k 
            print 'qsub sat.pbs -v N='+str(n)+',K='+str(k)
            os.system('qsub sat.pbs -v N='+str(n)+',K='+str(k) )
