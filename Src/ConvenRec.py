import sys

from readtmx import readtmx

if len(sys.argv)==1:
    fn='triples.tmx'
    numTJ=60000
elif len(sys.argv)==2:
    fn='triples.tmx'
    numTJ=int(sys.argv[1])
elif len(sys.argv)==3:
    fn=sys.argv[2]
    numTJ=int(sys.argv[1])
else:
    print('too many argments')

res=readtmx(fn,numTJ)
print res[0].shape

print res[1].shape

from scipy.sparse import csr_matrix
import numpy as np

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


xic=np.loadtxt('../Code/cap_vec.res')

X=xic.T.ravel().reshape((-1,1))


import myLOBPCG_new
X=xic.T.ravel().reshape((-1,1))
res=myLOBPCG_new.lobpcg(A,None,X,verbosityLevel=0,largest=False,maxiter=10000,retLambdaHistory=True)
print res[0]
print len(res[2])


a=res[1][:,0].reshape((-1,1))*486
tmp=np.hstack([a[:36*9**4],a[36*9**4:2*36*9**4],a[2*36*9**4:]])



np.save('MorRes.npy',tmp)

