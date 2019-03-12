# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 09:50:06 2019

@author: Home
"""
from read_file import read_file
from paraview import write_file
#from check import check_mat,check_print
import numpy as np


coeff_d=0.01
dt=0.1
U0=20
Itf=30 # Nb iterations

# Mesh creation from msh file:

our_mesh = read_file("../../maillage/square_holes.msh")
#parameters
our_mesh.init_cond(coeff_d,dt,U0)
our_mesh.maj_matrices()



Uit=np.array([])
for it in range(Itf):
    #print('Iteration : %d'%it)
    b=our_mesh.vector_b()
    U=our_mesh.vector_U()
    Uit=np.concatenate((Uit,U))
    our_mesh.t+=dt

Uit=Uit.reshape((Itf,our_mesh.Ns))


#Write solution in paraview format
write_file(our_mesh)
print(our_mesh.Uold)
#save matrices in files
#check_mat(our_mesh)

#print matrices 
#check_print(our_mesh)

