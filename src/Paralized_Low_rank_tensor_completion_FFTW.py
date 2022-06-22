# run svd in parallel
import numpy as np
import pandas as pd
import os
import scipy as sp
import math
from scipy import sparse
from numpy import random
import argparse
import time
from joblib import Parallel, delayed
import pyfftw
import ray
import logging

logging.basicConfig(
    format="%(asctime)s - %(levelname)s - %(name)s - %(message)s",
    datefmt="%m/%d/%Y %H:%M:%S",
    level=logging.INFO,
)
logger = logging.getLogger(__name__)

def Parser():
   parser=argparse.ArgumentParser('')
   parser.add_argument('-i', help = 'input file path')
   parser.add_argument('-o', help = 'output file path')
   parser.add_argument('-s', help = 'suffix')
   parser.add_argument('-t', help = 'tolerance of change of gradient',default = 1e-4,type=float)
   parser.add_argument('-max_iter', help = 'max iteration',default = 500,type=int)
   parser.add_argument('-mu', help = 'stepsize for dual variable',default = 1e-4,type=float)
   parser.add_argument('-max_mu', help = 'max stepsize for dual variable',default = 1e10,type=int)
   parser.add_argument('-rho', help = '>=1, increase mu',default = 1.1,type=float)
   parser.add_argument('-n_core',help= 'number of cores for paralization',default=10,type=int)
   return parser.parse_args()   

def tensor_fft(tensor):
    t_s = time.time()
    I = len(tensor)
    J = len(tensor[0])
    K = len(tensor[0][0])
    tmp_tensor = pyfftw.empty_aligned((I,J,K), dtype='complex128')
    fft_tensor = pyfftw.empty_aligned((I,J,K), dtype='complex128')
    fft_object = pyfftw.FFTW(tmp_tensor, fft_tensor,axes=(0,),threads=n_cores)
    tmp_tensor[:,:,:] = tensor
    fft_res = fft_object(tmp_tensor)
    t_e = time.time()
    logger.info("fft time:{}".format(t_e-t_s))
    return fft_res

def tensor_ifft(fft_tensor):
    t_s = time.time()
    I = len(fft_tensor)
    J = len(fft_tensor[0])
    K = len(fft_tensor[0][0])
    tmp_tensor = pyfftw.empty_aligned((I,J,K), dtype='complex128')
    ifft_tensor = pyfftw.empty_aligned((I,J,K), dtype='complex128')
    ifft_object = pyfftw.FFTW(tmp_tensor,ifft_tensor,axes=(0,), direction='FFTW_BACKWARD',threads=n_cores)
    tmp_tensor[:,:,:] = fft_tensor
    ifft_res = ifft_object(tmp_tensor)
    t_e = time.time()
    logger.info("ifft time:{}".format(t_e-t_s))
    return ifft_res

def prox_tnn(Y,rho):
    t_s = time.time()
    I,J,K=Y.shape
    X=np.zeros((I,J,K))
    X=X.astype('complex128')
    Y=Y.astype('complex128')
    Y=tensor_fft(Y)
    #Y=np.fft.fft(Y,axis=0)
    tnn=0
    trank=0
    U,S,V = np.linalg.svd(Y[0],hermitian=True)
    r = np.sum(S>rho)
    if r>=1:
        S=S[0:r]-rho
        X[0] = np.dot(np.dot(U[:,0:r],np.diag(S)),V[0:r,])
        tnn += sum(S)
        trank = max(trank,r)
    halfn3 = math.ceil(I/2)
    svd_list = ray_svd(Y[2:(halfn3+1)],n_cores)
    for i in range(2,halfn3+1):
        U,S,V = svd_list[i-2]
        r = np.sum(S>rho)
        if r>=1:
            S=S[0:r]-rho
            X[i] = np.dot(np.dot(U[:,0:r],np.diag(S)),V[0:r,])
            tnn += sum(S)*2
            trank = max(trank,r)
        X[I-i] = X[i].conjugate()
    
    if I%2 == 0:
        i = halfn3+1
        U,S,V = np.linalg.svd(Y[i],hermitian=True)
        r = np.sum(S>rho)
        if r>=1:
            S=S[0:r]-rho
            X[i] = np.dot(np.dot(U[:,0:r],np.diag(S)),V[0:r,])
            tnn += sum(S)
            trank = max(trank,r)
    tnn /= I
    X=tensor_ifft(X)
    #X=np.fft.ifft(X,axis=0)
    t_e = time.time()
    logger.info("prox_tnn time:{}".format(t_e-t_s))
    return X,tnn

def opt(M,omega,tol,max_iter,mu,max_mu,rho):
    I,J,K=M.shape
    X=np.zeros((I,J,K))
    X=X.astype('complex128')
    for i in omega:
        X[i] = M[i]  
    E=np.zeros((I,J,K))
    Y=E
    iter = 0
    log = []
    t_p=time.time()
    for iter in range(max_iter):
        Xk = X
        Ek = E
        X,tnnX=prox_tnn(-E+M+Y/mu,1/mu)
        E=M-X+Y/mu
        for i in omega:
            E[i] = 0
        dY = M-X-E
        D = I*J*K
        chgX = np.max( np.abs(Xk.reshape(D,)-X.reshape(D,)) )
        chgE = np.max( np.abs(Ek.reshape(D,)-E.reshape(D,)) )
        chg = np.max(  [ chgX,chgE,np.max(np.abs(dY.reshape(D,)))  ] )
        if chg < tol:
            break
        Y = Y+mu*dY
        mu = min(rho*mu,max_mu)
        #logger.info(X)
        logger.info("Iter:%0.2f" % iter)
        logger.info("Error: %0.2f" % np.linalg.norm( dY.reshape(D,) ))
        log.append([iter,np.linalg.norm( dY.reshape(D,) ),chg])
        t_e = time.time()
        logger.info(t_e-t_p)
        t_p=t_e
    obj = tnnX
    err = np.linalg.norm( dY.reshape(D,) )
    return X,log

@ray.remote
def mysvd(mat,idx):
    if idx%50 == 0:
        logger.info('Cell {} starts'.format(idx))
    #logger.info('Cell {} starts'.format(idx))
    res = np.linalg.svd(mat,hermitian=True)
    return res

def par_svd(tensor,n_core):
    res = Parallel(n_jobs=n_core)( delayed(mysvd)(i,idx) for idx,i in enumerate(tensor))
    return res

def ray_svd(tensor,n_core):
    ray.put(tensor)
    N=tensor.shape[0]
    res = ray.get([mysvd.remote(i,idx) for idx,i in enumerate(tensor) ])
    return res

args = Parser()
input_path = args.i
output_path = args.o
tol = args.t
max_iter = args.max_iter
mu = args.mu
max_mu = args.max_mu
rho = args.rho
suffix = args.s
n_cores = args.n_core


#initialize ray
ray.init(num_cpus=n_cores)
# read in data and store as a tensor
sc_tensor = []
all_file = os.listdir(input_path)
#KR_file = [i for i in all_file if 'RawCount' in i]
KR_file = all_file
KR_file.sort()
#KR_file=[]
#for i in range(1,41):
#   KR_file.append('slice_'+str(i)+'.txt')

# logger.info raw data file name
logger.info('These files are in use:\n')
[logger.info(i) for i in KR_file]
KR_file_frame = pd.DataFrame(KR_file)
KR_file_frame.to_csv(output_path+'/'+suffix+'_input_file_index.csv')
sc_tensor = pd.read_csv(input_path+"/"+KR_file[0],header=None,sep='\t')
sc_tensor = sc_tensor.values
N = sc_tensor.shape[0]
sc_tensor = np.empty(shape=[len(KR_file),N,N])

for idx,f in enumerate(KR_file):
    tmp_file = pd.read_csv(input_path+"/"+f,header=None,sep='\t')
    tmp_file = tmp_file.values
    # check if the input data has same dimensions
    if idx == 0:
        N = len(tmp_file)
    else:
        if len(tmp_file) != N:
            raise ValueError('scHi-C in this chromosome has different dimensions')
    sc_tensor[idx,:,:]=tmp_file
    logger.info('{idx} outof {n_file} done'.format(idx=idx,n_file=len(KR_file)))

sc_tensor=np.array(sc_tensor) # A[i,j,k] represent (j,k) element of single cell i
# define measurement set: non zero entries in sc_tensor
I = len(sc_tensor)
J = len(sc_tensor[0])
K = J
omega = set()
non_zero_entries=np.transpose(sc_tensor.nonzero())
for i,j,k in non_zero_entries:
    omega.add( (i,k,j) )
    
logger.info(len(omega))
X,log = opt(sc_tensor,omega,tol,max_iter,mu,max_mu,rho)
np.save(output_path+'/'+suffix, X)
log = pd.DataFrame(log)
log.to_csv(output_path+'/'+suffix+".txt",sep='\t')

