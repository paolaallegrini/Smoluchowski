# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 09:50:06 2019

@author: Home
"""
from read_file import read_file
from paraview import write_file
from plot import plot_mesh
#from check import check_mat,check_print
import numpy as np

dt=0.1
coeff_d=0.5
U0=10.0
Itf=15 # Nb iterations

# Mesh creation from msh file:

our_mesh = read_file("../../maillage/square_4_borders.msh")
#parameters
our_mesh.init_cond(coeff_d,dt,U0)
our_mesh.maj_matrices()




''' Time loop '''
Uit=np.array([])
for it in range(Itf):
    #print('Iteration : %d'%it)
    U=our_mesh.vector_U()
    Uit=np.concatenate((Uit,U))
    our_mesh.t+=dt

Uit=Uit.reshape((Itf,our_mesh.Ns))
#np.savetxt("mat_Uit.csv",Uit,delimiter=",")


#for it in range(0,np.size(our_mesh.Uold)):
#    print('Uold({})={}'.format(it, our_mesh.Uold[it]))

''' Plot mesh '''
plot_mesh(our_mesh,Uit)

'''Write solution in paraview format'''
write_file(our_mesh)



#save matrices in files
#check_mat(our_mesh)
#print matrices 
#check_print(our_mesh)

