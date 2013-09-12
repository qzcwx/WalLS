# spin glasses class

cdef class SpinGlass:
    cdef public int n              # number of variables
    cdef public int dim
    cdef public int m              # number of subfunctions
    cdef public int k              # number of variables per subfunction

    def __init__(self, inSize):
        print inSize
