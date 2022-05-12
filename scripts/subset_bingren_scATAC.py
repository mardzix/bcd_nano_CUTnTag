import pandas as pd
import sys
import random
import math

matrix_file_path = sys.argv[1]          # Path to input matrix
ncells           = sys.argv[2]          # either ncells [>1] or fraction of dataset to return [0-1]
matrix_file_out  = sys.argv[3]          # Path to the output file

sys.stdout.write("Reading the dataset\n")
colnames          = list(pd.read_csv(matrix_file_path,sep='\t',nrows=1))
try:
    ncells = int(ncells)
except ValueError:
    ncells = float(ncells)

if ncells > 0 and ncells < 1:
    ncells = math.floor(float(ncells) * len(colnames))

sys.stdout.write("Filtering {} cells from the dataset\n".format(ncells))

colnames_selected = random.sample(colnames,ncells)
colnames_selected = [colnames[0]] + colnames_selected
del(colnames)

sys.stdout.write("Subseting the matrix \n")
matrix_small = pd.read_csv(matrix_file_path,sep='\t',usecols=colnames_selected)

sys.stdout.write("Exporting the data to csv  \n")
matrix_small.to_csv(matrix_file_out,sep='\t',index=False)