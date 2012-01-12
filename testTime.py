"""
Test the runtime for rLS with mean of neighborhoods as fitness, 
the time for computing walsh coefficients is stored as well.
"""
import os

for k in [2, 4, 8]:
    for n in [200, 500]:
        for i in [0]:
            for a in ['rLS']:
                for f in ['fit']:
                    for c in ['bf']:
                        for p in ['NK','NKQ']:
                            if p == 'NKQ':
                                for q in [2]:
                                    s = 1
                                    print 'python demo.py '+c+' '+p+' '+a+' '+f+' '+str(i)+' '+str(s)+' '+str(n)+' '+str(k)+' '+str(q)
                                    os.system('python demo.py '+c+' '+p+' '+a+' '+f+' '+str(i)+' '+str(s)+' '+str(n)+' '+str(k)+' '+str(q))
                            elif p == 'NK': # for NK problem
                                q = 0
                                s = 1
                                print 'python demo.py '+c+' '+p+' '+a+' '+f+' '+str(i)+' '+str(s)+' '+str(n)+' '+str(k)+' '+str(q)
                                os.system('python demo.py '+c+' '+p+' '+a+' '+f+' '+str(i)+' '+str(s)+' '+str(n)+' '+str(k)+' '+str(q))
