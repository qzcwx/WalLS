# this function counts the number of local optimum over the entire search space
import random
from libcpp.vector cimport vector
from libc.math cimport pow
from cython.operator cimport dereference as deref, preincrement as inc
from libcpp cimport bool
from libcpp.set cimport set
#from libcpp cimport advance

cpdef int localOpt(bitStr, f):
    cdef int i, j, solIndex, neigh, bitStrLen = len(bitStr), weigh, n, numLocalOpt
    cdef vector[int] loMark
    cdef bool localOpt

    n = len(bitStr[0])
    numLocalOpt = 0

    """ initialize all markers to true """
    for i in xrange( bitStrLen ):
        loMark.push_back(i)
        
    """ generate the random solution under consideration for local optimum """
    while loMark.empty()==False :
        # randomly pick a solution
        solIndex = loMark[random.randint(0,loMark.size())]
        localOpt = True
        # print 'solIndex', solIndex, 'loMark size', loMark.size() # 
        """ considering all neighbors of solIndex """
        for i in xrange(n):
            weigh = <int> pow(2, i)
            # print 'weigh', weigh
            # print 'bit at position', i, 'is', solIndex & weigh
            if solIndex & weigh != 0:
                # the bit sitting at position i is 1, set it to 0
                neigh = solIndex - weigh
            else :
                # the bit sitting at position i is 0, set it to 1
                neigh = solIndex + weigh
            if f[solIndex]>f[neigh]:
                # print 'remove sol', solIndex
                remove(loMark, solIndex)
                # print 'greater'
                localOpt = False
                break
            elif f[solIndex] < f[neigh]:
                remove(loMark, neigh)
        if localOpt == True:
            numLocalOpt = numLocalOpt + 1
            remove(loMark, solIndex)

    return numLocalOpt

cpdef plateauCount(bitStr, f, opt):
    """ count the number plateau by Enumeration """
    cdef int i, bitStrLen = len(bitStr), solIndex, n
    cdef vector[int] mark
#    cdef set[int] realLocOpt
    cdef vector[set[int]] platSet, exitSet
    n = len(bitStr[0])

    """ initialize all markers to true """
    for i in xrange( bitStrLen ):
        mark.push_back(i)

    while mark.empty() == False :
        # randomly pick a solution
        solIndex = mark[random.randint(0,mark.size())]
        # nonplat = True
        # consider all neighs of solIndex 
        for i in xrange(n):
            weigh = <int> pow(2, i)
            if solIndex & weigh != 0:
                # the bit sitting at position i is 1, set it to 0
                neigh = solIndex - weigh
            else :
                # the bit sitting at position i is 0, set it to 1
                neigh = solIndex + weigh
               # print 'solIndex', solIndex, 'neigh', neigh
            if f[solIndex] == f[neigh]:
                # if nonplat == True:
                #     remove(mark, solIndex)
                #     nonplat = False
                # both solIndex and neigh have been inspected
#                remove(mark, neigh)

                # merge two set of plateaus relative to solIndex and neigh,
                # if not exists yet, create a plateau containing two points
                merge(platSet, solIndex, neigh) 
 #               print 'merge'
            
        remove(mark, solIndex)
    # for i in xrange(platSet.size()):
    #      print "%g\t%g" %(platSet[i].size(), f[deref(platSet[i].begin())])

    # print 'exitConut'
    exitSet = exitCount(platSet, f, n)

    writeToFile(platSet, exitSet, realLocalOptimaCount(platSet, f, n, bitStrLen), f, opt)

cdef vector[set[int]] exitCount(vector[set[int]] plat, object f, int n):
    """ count the number of exists out of all plateaus

    exist are those plateaus that have an improving neighborhoods,
    think about improving move as the the one with higher fitness,
    just to be consistent with MaxSAT problem
    """
    # make a copy of platSet
    cdef vector[set[int]] exitSet
    cdef int i, weigh, solIndex, neigh
    cdef set[int].iterator it
    cdef set[int] tempS
    cdef bool impr
    
    for i in xrange(plat.size()):
        it = plat[i].begin()
        tempS.clear()
        while it != plat[i].end():
            tempS.insert(deref(it))
            inc(it)
        exitSet.push_back(tempS)

    # for i in xrange(exitSet.size()):
    #     print "%g\t%g" %(exitSet[i].size(), f[deref(exitSet[i].

    # check the neighborhood for every points in the plateauos
    for i in xrange(exitSet.size()):
        it = exitSet[i].begin()
        while it != exitSet[i].end():
            solIndex = deref(it)
    #        print 'sol', solIndex
            impr = False
            for j in xrange(n):
                weigh = <int> pow(2,j)
                if solIndex & weigh != 0:
                    # the bit sitting at position i is 1, set it to 0
                    neigh = solIndex - weigh
                else :
                    # the bit sitting at position i is 0, set it to 1
                    neigh = solIndex + weigh
   #             print 'neigh', neigh
                if f[solIndex] < f[neigh]:
                    impr = True
  #                  print 'impr', neigh
                    break
            if impr == False:
 #               print 'erase', solIndex
                exitSet[i].erase(it)
            inc(it)
#            print

    return exitSet
  
    # for i in xrange(exitSet.size()):
    #     print "%g\t%g" %(exitSet[i].size(), f[deref(exitSet[i].begin())])

cdef set[int] realLocalOptimaCount(vector[set[int]] plat, object f, int n, int bitStrLen):
    """ count the number of real optima, those points that have no equal moves nor improving

    one can guarantee that non of real optima will be on plateau, thus we only need to consider
    the rest of elements other than those on the plateau
    """
    cdef set[int] mark                    # possible real local optima that are not on plateau
    cdef int i, c=0, solIndex
    cdef set[int].iterator it

    for i in xrange( bitStrLen ):
        mark.insert(i)
    
    for i in xrange(plat.size()):
        it  = plat[i].begin()
        while it != plat[i].end():
            c = c + 1
            mark.erase(mark.find(deref(it)))
            inc(it)

    it = mark.begin()
    while it != mark.end():
        solIndex = deref(it)
        for j in xrange(n):
            weigh = <int> pow(2,j)
            if solIndex & weigh != 0:
                # the bit sitting at position i is 1, set it to 0
                neigh = solIndex - weigh
            else :
                # the bit sitting at position i is 0, set it to 1
                neigh = solIndex + weigh
   #             print 'neigh', neigh
                if f[solIndex] < f[neigh]: # find an improving move, no longer a potential real local optima
                    mark.erase(it)
                    break
        inc(it)

    
cdef void writeToFile(vector[set[int]] plat, vector[set[int]] exitSet, set[int] realLocOpt, object fit, object opt):
    cdef str dirPrefix = 'plateau/', realLocOptPrefix = 'plateau/realLocOpt-'
    cdef set[int].iterator it
    
    # if opt.probName == 'NKQ':
    #     nameOfF = dirPrefix+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-C'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-N'+str(opt.n)+'-K'+str(opt.k)+'-Q'+str(opt.q)+'-T'+str(opt.t)+'.txt'
    # elif opt.probName == 'NK':
    #     nameOfF = dirPrefix+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-C'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-N'+str(opt.n)+'-K'+str(opt.k)+'.txt'

    # store the plateau info to files
    if opt.probName == 'NKQ':           
        nameOfF = dirPrefix+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-W'+str(opt.w)+'-N'+str(opt.n)+'-K'+str(opt.k)+'-C'+str(opt.c)+'-Q'+str(opt.q)+'-T'+str(opt.t)+'.txt'
    elif opt.probName == 'NK':
        nameOfF = dirPrefix+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-W'+str(opt.w)+'-N'+str(opt.n)+'-K'+str(opt.k)+'-C'+str(opt.c)+'.txt'

    f = open(nameOfF, 'w')
    print >>f, "PlateauSize\tLocalOptimaSize\tExitSize\tFitness"
    for i in xrange(plat.size()):
        print >>f, "%g\t%g\t%g\t%g" %(plat[i].size(), plat[i].size()-exitSet[i].size(), exitSet[i].size(), fit[deref(plat[i].begin())])
    f.close()

    # store the real local optima to files
    if opt.probName == 'NKQ':
        nameOfF = realLocOptPrefix+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-W'+str(opt.w)+'-N'+str(opt.n)+'-K'+str(opt.k)+'-C'+str(opt.c)+'-Q'+str(opt.q)+'-T'+str(opt.t)+'.txt'
    elif opt.probName == 'NK':
        nameOfF = realLocOptPrefix+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-M'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-W'+str(opt.w)+'-N'+str(opt.n)+'-K'+str(opt.k)+'-C'+str(opt.c)+'.txt'
        
    f = open(nameOfF, 'w')
    print >>f, "Fitness"
    it = realLocOpt.begin()
    while it != realLocOpt.end():
        print >>f, f[deref(it)]
        inc(it)
    f.close()
    

    
    
cdef void merge(vector[set[int]] &platSet, int e1, int e2):
    """ merge two set containing e1 and e2, respectively """
    cdef int i, set1=-1, set2=-1
    cdef bool s1Exist=False, s2Exist=False
    cdef set[int].iterator it
    cdef vector[set[int]].iterator vit
    cdef set[int] tempS

    # locate the set containing e1
    for i in xrange(platSet.size()):
        if (platSet[i].find(e1) != platSet[i].end()):
            set1 = i
            s1Exist = True
            break

    # locate the set containing e2
    for i in xrange(platSet.size()):
        if (platSet[i].find(e2) != platSet[i].end()):
            set2 = i
            s2Exist = True
            break

#    print s1Exist, s2Exist
    
    if s1Exist & s2Exist == True and (set1 != set2):
        # both e1 and e2 are contained in some sets
        # merge them to the first one, and delete the second one
        it = platSet[set2].begin()
        
        while it != platSet[set2].end():
            platSet[set1].insert(deref(it))
            inc(it)
        platSet[set2].clear()
        vit = platSet.begin()
        #        advance(vit, set2)
        #        vit = vit + set2
        for i in xrange(set2):
            inc(vit)
        platSet.erase(vit)
    elif s1Exist == True and s2Exist == False:
        # e2 is a newly discovered plateau, add it to s1
        platSet[set1].insert(e2)
    elif s1Exist == False and s2Exist == True:
        # e1 is a newly discovered plateau, add it to s2
        platSet[set2].insert(e1)
    elif set1 == -1 and set2 == -1: # both points have not existed in set list, create a new set containing both points and append the created to vectorx
        tempS.insert(e1)
        tempS.insert(e2)
        platSet.push_back(tempS)
        
cdef void remove(vector[int] &vec, int elem):
    cdef vector[int].iterator it
    it = vec.begin()
    while it != vec.end():
        if deref(it) == elem:
            vec.erase(it)
            break
        inc(it)
