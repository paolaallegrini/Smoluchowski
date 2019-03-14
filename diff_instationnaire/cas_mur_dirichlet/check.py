# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 14:45:14 2019

@author: Home
"""
import numpy as np

''' print in terminal  A, M, R and U  '''
def check_print(our_mesh):
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
    
''' save in csv file  matrices A, M and R '''
def check_mat(our_mesh):
    A=our_mesh.matrice_A()
    M=our_mesh.matrice_mass()
    R=our_mesh.matrice_rigidite()
    
    np.savetxt("mat_A_py.csv",A,delimiter=",")
    np.savetxt("mat_M_py.csv",M,delimiter=",")
    np.savetxt("mat_R_py.csv",R,delimiter=",")
    return

''' range with floating points'''
def frange(start, stop, step):
    i=start
    while i <stop:
        yield i
    i+=step