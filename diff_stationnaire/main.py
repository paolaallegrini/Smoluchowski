# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 09:50:06 2019

@author: Home
"""
from read_file import read_file
from paraview import write_file
import numpy as np

def check_mat():
    A=our_mesh.matrice_A()
    M=our_mesh.matrice_mass()
    R=our_mesh.matrice_rigidite()
    
    np.savetxt("mat_A_py.csv",A,delimiter=",")
    np.savetxt("mat_M_py.csv",M,delimiter=",")
    np.savetxt("mat_R_py.csv",R,delimiter=",")
    return
# Mesh creation from msh file:
our_mesh = read_file("maillage/square_simple.msh")

#save matrices in files
check_mat()

#print matrices in terminal
def check_print():
    print("------------------A--------------")
    A = our_mesh.matrice_A()
    print(A)
    print("------------------det(A)--------------")
    print(np.linalg.det(A))
    print("------------------M--------------")
    print(our_mesh.matrice_mass())
    print("------------------R--------------")
    print(our_mesh.matrice_rigidite())
    print("------------------U--------------")
    U=our_mesh.vector_U()
    print(U)
    print("--------------Done---------------")

#Write solution in paraview format
write_file(our_mesh)