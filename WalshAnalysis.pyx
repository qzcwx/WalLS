# Walsh Analysis Module
import numpy as np

cpdef computeW(i,f):
    """ Walsh coefficient computation 

        This function receives a tuple of bit-string and a tuple of 
        corresponding function value, return the Walsh coefficient.
    """
    cdef int L
    cdef double factor
    L = len(i[0])
    factor = (1.0 / len(f))
    # return np.array([ factor * subsum(i,i[j],f) for j in range(len(i)) ])
    return [ factor * subsum(i,i[j],f) for j in range(len(i)) ]

cdef subsum(list i,str j,list f):
    """ compute the subsum for j bit-string """
    cdef double s = 0
    cdef int index
    for index in xrange(len(i)):
       if bc(i[index],j)%2 == 0:
           s = s + f[index]
       else:
           s = s - f[index]
    return s

cdef bc(str x,str y):
    """ return the bit count of x AND y """
#    return bin( int(x,2) & int(y,2) ).count('1')
    cdef int sum, i
    sum = 0
    for i in xrange(len(x)):
        if x[i] == '1' and y[i] == '1':
           sum = sum + 1
    return sum 
