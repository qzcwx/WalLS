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
    """ count the number plateau using BFS """
    cdef int i, bitStrLen = len(bitStr), solIndex, n
    cdef vector[int] mark
    cdef vector[set[int]] platSet
    cdef nonplat
    n = len(bitStr[0])

    """ initialize all markers to true """
    for i in xrange( bitStrLen ):
        mark.push_back(i)

    while mark.empty() == False :
        # randomly pick a solution
        solIndex = mark[random.randint(0,mark.size())]
        nonplat = True
        # consider all neighs of solIndex 
        for i in xrange(n):
            weigh = <int> pow(2, i)
            if solIndex & weigh != 0:
                # the bit sitting at position i is 1, set it to 0
                neigh = solIndex - weigh
            else :
                # the bit sitting at position i is 0, set it to 1
                neigh = solIndex + weigh
            # print 'solIndex', solIndex, 'neigh', neigh, 'mark.size', mark.size() 
            if f[solIndex] == f[neigh]:
                if nonplat==True:
                    remove(mark, solIndex)
                    nonplat = False
                # both solIndex and neigh have been inspected
                remove(mark, neigh)
                # merge two set of plateaus relative to solIndex and neigh,
                # if not exists yet, create a plateau containing two points
                merge(platSet, solIndex, neigh)
        if nonplat == True:
            remove(mark, solIndex)

    writeToFile(platSet, f, opt)
    return platSet.size()

cdef void writeToFile(vector[set[int]] plat,object fit, object opt):
    cdef str dirPrefix = 'plateau/'
    if opt.probName == 'NKQ':
        nameOfF = dirPrefix+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-C'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-N'+str(opt.n)+'-K'+str(opt.k)+'-Q'+str(opt.q)+'-T'+str(opt.t)+'.txt'
    elif opt.probName == 'NK':
        nameOfF = dirPrefix+opt.probName+'-'+opt.algoName+'-F'+opt.fitName+'-C'+opt.compMeth+'-I'+str(opt.inst)+'-S'+str(opt.s)+'-N'+str(opt.n)+'-K'+str(opt.k)+'.txt'

    f = open(nameOfF, 'w')
    print >>f, "Plateau Size\tFitness"
    for i in xrange(plat.size()):
        print >>f, "%g\t%g" %(plat[i].size(), fit[deref(plat[i].begin())])
    f.close()
    
cdef void merge(vector[set[int]] &platSet, int e1, int e2):
    """ merge two set containing e1 and e2, respectively """
    cdef int i, set1, set2
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

    if s1Exist & s2Exist == True:
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
    else: # both points have not existed in set list, create a new set containing both points and append the created to vectorx
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
