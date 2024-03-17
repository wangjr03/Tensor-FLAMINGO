import pandas as pd
import numpy as np
import argparse
import os
def Parser():
   parser=argparse.ArgumentParser('')
   parser.add_argument('-i', help = 'input npy path')
   parser.add_argument('-o', help = 'output file path')
   parser.add_argument('-t', default=1, help = 'minimum imputed IF for valid data')
   parser.add_argument('-alpha', default=0.25, help = 'alpha')
   return parser.parse_args()   
args = Parser()
input_path = args.i
output_path = args.o
dat = np.load(input_path)
threshold = args.t
alpha = args.alpha
os.mkdir(output_path)

for idx,mat in enumerate(dat):
    tmp_mat = np.real(mat)
    if_mat = tmp_mat
    if_mat[if_mat<threshold] = np.nan
    pd_mat = if_mat**(-alpha)
    np.savetxt(output_path+'/IF_Cell_'+str(idx+1)+'.txt',if_mat)
    np.savetxt(output_path+'/PD_Cell_'+str(idx+1)+'.txt',pd_mat)
