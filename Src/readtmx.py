#############################################################################
from scipy.io import FortranFile
import numpy as np


def readtmx(fname='triples.tmx',numTJ=40188):
    """
    Parameter:
    numTJ:  int
            number of triple junctions
    Returns:
    huge_kcell: ndarray (numTJ,3,36)
            Cell indices that have non-zero coefficient in A, data type in this
            array is integer. First dimension is the TJ index, second is the 
            boundary index, third is the equivalent point. Using notations in 
            paper, it is (J,s,?)
    huge_coeff: ndarray (numTJ,3,36,3,3)
            Coefficients in A, data type is float. Using notations in paper, it
            is (J,s,?,i,l)
    """
    inttype=np.int32
    floattype=np.float32
    #headertype=np.uint32
    headertype=np.uint64

    f = FortranFile(fname, 'r',header_dtype=headertype)

    number_of_points=numTJ

    huge_coeff=np.empty((number_of_points,3,36,3,3)).astype(floattype)
    huge_kcell=np.empty((number_of_points,3,36)).astype(inttype)


    for nn in range(number_of_points):
        first=f.read_ints(dtype=inttype)
        for ss in range(3):
            huge_kcell[nn,ss]=f.read_ints(dtype=inttype)
            for kk in range(36):
                #Tested: Fortran has the weird column first order, python has row first.
                #The transpose makes huge_coeff[nn,ss,kk,i,l] equals the xcoeffo(ss,kk,i,l) in Morawiec's code
                huge_coeff[nn,ss,kk]=f.read_record(dtype=floattype).reshape((3,3)).transpose()
        for kk in range(36):
            for ii in range(3):
                huge_coeff[nn,:,kk,ii,:]=huge_coeff[nn,:,kk,ii,:]/(np.sum(huge_coeff[nn,:,kk,ii,:]**2))**0.5 
    return huge_kcell,huge_coeff
###############################################################################


