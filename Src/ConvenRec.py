import myLOBPCG_new
from scipy.sparse import csr_matrix
from scipy.io import FortranFile
import numpy as np

prefix='../output/' # output directory
fn='triples.tmx' # input .tmx file
numTJ=60000 # number of TJs in .tmx file
xic=np.loadtxt('cap_vec.res') # file that contains the plane normal for each cell


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
                #Tested: Fortran is the column first order, python is row first.
                #The transpose makes huge_coeff[nn,ss,kk,i,l] equals the xcoeffo(ss,kk,i,l) in Morawiec's reconstruction code
                huge_coeff[nn,ss,kk]=f.read_record(dtype=floattype).reshape((3,3)).transpose()
        for kk in range(36):
            for ii in range(3):
                huge_coeff[nn,:,kk,ii,:]=huge_coeff[nn,:,kk,ii,:]/(np.sum(huge_coeff[nn,:,kk,ii,:]**2))**0.5 
    return huge_kcell,huge_coeff
###############################################################################




res=readtmx(fn,numTJ)


data=np.swapaxes(res[1],1,3).ravel()
row_ind=np.tile(np.arange(res[0].shape[0]*3),(36*3*3,1)).T.ravel()

tmp=np.swapaxes(np.swapaxes(np.tile(res[0],(3,3,1,1,1)),0,2),2,4)
tmp[:,:,:,:,1]+=(36*9**4)
tmp[:,:,:,:,2]+=(2*36*9**4)
col_ind=tmp.ravel()-1

A=csr_matrix((data,(row_ind,col_ind)),shape=(res[0].shape[0]*3,3*36*9**4))
del data
del row_ind
del col_ind
del tmp
del res



X=xic.T.ravel().reshape((-1,1))


X=xic.T.ravel().reshape((-1,1))
res=myLOBPCG_new.lobpcg(A,None,X,verbosityLevel=0,largest=False,maxiter=10000,retLambdaHistory=True)
print("number of lobpcg iterations: ",len(res2[2]))


a=res[1][:,0].reshape((-1,1))*486
tmp=np.hstack([a[:36*9**4],a[36*9**4:2*36*9**4],a[2*36*9**4:]])

energy=np.sum(xic*tmp,axis=1)
if(np.sum(energy<0)/float(len(energy))>0.5):
    energy=-energy

np.savetxt(prefix+'cap_vec.res',tmp)
np.savetxt(prefix+'energy_v.res',energy)

