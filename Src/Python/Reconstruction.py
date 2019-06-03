import numpy as np
from scipy.sparse import csr_matrix
import myLOBPCG_new

###########################################################################
#
#  Parameters and file path
#
###########################################################################

TJfile='../../TJdata/triples_30000.dat' 
prefix='../../output/' # directory contains rowA.binary, colA.binary, valA.binary; and will save outputs
numTJ=30000 # number of triple junctions
lamb=1000 # hyperparameter for the strength of the regularization
sym='Cubic' # Cubic or Hex, it changes the gbdat file header 
fn= prefix+'Cub.gbdat' # the name of output gbdat file

###########################################################################
#
# Define util functions
#
###########################################################################


def read_dat(datFile, numTJ):
    """
    Input: triples.dat, wrote from the fortran program Torq_trn
                size=[numTJ*8,]
                In each group, the data is [TJ directon, EA1, GB1, EA2, GB2, EA3, GB3]
    Output: TJs, direction of the triple junctions
                size = [numTJ, 3]
            EAs, the EA angles of the 3 grains at a TJ
                size = [numTJ, 3, 3]
            norms, normal direction of the 3 GB at a TJ
                size = [numTJ, 3, 3]
    """
    with open(datFile) as f:
        tmp = [line.split() for line in f if line.strip()]
    TJs = np.zeros((numTJ, 3))
    EAs = np.zeros((numTJ, 3, 3))
    norms = np.zeros((numTJ, 3, 3))
    for i in range(numTJ):

        TJs[i,:] = np.array(tmp[i*8 + 1]).astype(float)
        EAs[i,0, :] = np.array(tmp[i*8 + 2]).astype(float)
        norms[i,0, :] = np.array(tmp[i*8 + 3]).astype(float)
        EAs[i, 1, :] = np.array(tmp[i*8 + 4]).astype(float)
        norms[i, 1, :] = np.array(tmp[i*8 + 5]).astype(float)
        EAs[i, 2, :] = np.array(tmp[i*8 + 6]).astype(float)
        norms[i, 2, :] = np.array(tmp[i*8 + 7]).astype(float)

    return (TJs, EAs, norms)

def EulerZXZ2Mat(e):
    """
    Active Euler Angle (radian)  in ZXZ convention to active rotation matrix, which means newV=M*oldV
    """
    x=e[0]
    y=e[1]
    z=e[2]
    s1=np.sin(x)
    s2=np.sin(y)
    s3=np.sin(z)
    c1=np.cos(x)
    c2=np.cos(y)
    c3=np.cos(z)
    m=np.array([[c1*c3-c2*s1*s3,-c1*s3-c3*c2*s1,s1*s2],
        [s1*c3+c2*c1*s3,c1*c2*c3-s1*s3,-c1*s2],
        [s3*s2,s2*c3,c2]])
    return m

def EAtoG(EA):
    """
    Input: a set of Euler Angle
                size=[3,]
    Output: the corresponding orientation matrix g
                size = [3, 3]
    """
    g = np.zeros((3,3))
    EA = np.radians(EA)
    g=EulerZXZ2Mat(EA).T
    return g

###########################################################################
#
#  Construct and solve the minimization problem to get GB energy
#
###########################################################################



(TJs, EAs, norms) = read_dat(TJfile, numTJ)

Norm=np.empty((numTJ*3,3))
for i in range(numTJ):
    Norm[3*i]=EAtoG(EAs[i,1]).dot(norms[i,0])
    Norm[3*i+1]=EAtoG(EAs[i,2]).dot(norms[i,1])
    Norm[3*i+2]=EAtoG(EAs[i,0]).dot(norms[i,2])
for j in range(len(Norm)):
    Norm[j]=Norm[j]/(np.linalg.norm(Norm[j]))


rowB=[]
colB=[]
valB=[]

for i in range(numTJ):
    t=TJs[i]
    xls=np.array([[0,-t[2],t[1]],[t[2],0,-t[0]],[-t[1],t[0],0]])
    M=np.hstack([EAtoG(EAs[i,1]).T,EAtoG(EAs[i,2]).T,EAtoG(EAs[i,0]).T])
    M=xls.dot(M)
    rowB.extend([i*3]*9)
    rowB.extend([i*3+1]*9)
    rowB.extend([i*3+2]*9)
    colB.extend(list(range(i*9,i*9+9))*3)
    valB.extend(M.ravel())

B=csr_matrix((valB,(rowB,colB)),shape=(rowB[-1]+1,numTJ*3*3))
rowB=None
colB=None
valB=None




rowA=np.fromfile(prefix+'rowA.binary',dtype=np.uint32)
colA=np.fromfile(prefix+'colA.binary',dtype=np.uint32)
valA=np.fromfile(prefix+'valA.binary',dtype=np.float32)
A=csr_matrix((valA,(rowA,colA)),shape=(rowA[-1]+1,numTJ*3*3))
rowA=None
colA=None
valA=None



X0=Norm.ravel().reshape(-1,1)
print("initial loss of |AX|^2 (smoothness): \n", np.sum((A.dot(X0))**2))
print("initial loss of |BX|^2 (consistency with Herring's equation): \n", np.sum((B.dot(X0))**2))
res2=myLOBPCG_new.lobpcg(A,lamb*B,X0,verbosityLevel=0,largest=False,maxiter=500,retLambdaHistory=True)
print("number of lobpcg iterations: ",len(res2[2]))
print("residual loss of |AX|^2: ", np.sum((A.dot(res2[1]))**2))
print("residual loss of |BX|^2: ", np.sum((B.dot(res2[1]))**2))



resxi=np.linalg.norm(Norm.ravel())/np.linalg.norm(res2[1])*res2[1].reshape((-1,3))
resE=np.sum(resxi*Norm,axis=1)

if np.sum(resE<0)/float(len(resE))>0.5:
    resxi=-resxi
    resE=-resE
    
print("fraction of negative energy: ",np.sum(resE<0)/float(len(resE)))

np.savetxt(prefix+'resE.txt',resE)


###########################################################################
#
#  save .gbdat file for plotting
#
###########################################################################

tmp=norms
L_EA = np.zeros((3*numTJ, 3))
R_EA = np.zeros((3*numTJ, 3))
n  = np.zeros((3*numTJ, 3))
for i in range(numTJ):
        L_EA[3*i,:] = EAs[i,1,:]
        R_EA[3*i,:] = EAs[i,2,:]
        n[3*i,:] = tmp[i,0,:]
        L_EA[3*i + 1,:] = EAs[i,2,:]
        R_EA[3*i + 1,:]  = EAs[i,0,:]
        n[3*i + 1,:] =tmp[i,1,:]
        L_EA[3*i + 2,:] = EAs[i,0,:]
        R_EA[3*i + 2,:] = EAs[i,1,:]
        n[3*i + 2,:] = tmp[i,2,:]   

PHI=np.degrees(np.arccos(n[:,2])).reshape((n.shape[0],1))
THETA=np.degrees(np.arctan2(n[:,1],n[:,0])).reshape((n.shape[0],1))
THETA[THETA<0] += 360

doomy = np.ones((len(PHI),1))*1200
resE2t = resE.reshape((len(PHI),1))

toPrint = np.hstack((L_EA, R_EA, PHI, THETA, doomy, resE2t))
with open(fn,'w') as text_file:
    if sym=='Hex':
        text_file.write(
    """# This file was created by GBToolbox 1.1.20151207
# It contains boundary parameters imported from DREAM.3D output files
# EXP
# 6/mmm
# L_PHI1 L_PHI L_PHI2 R_PHI1 R_PHI R_PHI2 ZENITH AZIMUTH CORRELAT AREA\n""")
    elif sym=='Cubic':
        text_file.write(
    """# This file was created by GBToolbox 1.1.20151207
# It contains boundary parameters imported from DREAM.3D output files
# EXP
# m-3m
# L_PHI1 L_PHI L_PHI2 R_PHI1 R_PHI R_PHI2 ZENITH AZIMUTH CORRELAT AREA\n""")
    else:
        print("sym must be Cubic or Hex")
np.savetxt(open(fn,'a'), toPrint, delimiter=',', fmt='%3.4f '*8+'%4d '+'%3.7f')




