# -*- coding: utf-8 -*-
"""
Created on Tue May  7 10:03:13 2019

@author: Home
"""

import numpy as np 
from read_file import read_file
from paraview import write_file,erase_files

'Mesh creation from msh file'

our_mesh = read_file("C:/Users/Home/Desktop/stage_labo/Smoluchowski/maillage/square_4_borders_hole.msh")
erase_files()


'''parameters'''
dt=0.1
NB=11
coeff_d=[ 1.0/(i**(1./3.)) for i in range(1,NB+1)]
#coeff_d=[ 0.1 for i in range(1,NB+1)]
Ut0=np.zeros(NB)
Tf=3
Itf=5#int(Tf/dt) # Nb iterations
UM=[]
Uold,U=our_mesh.init_cond(coeff_d,dt,Ut0)

''' Initial situation '''
our_mesh.maj_matrices()
cl=5
for it in range(Itf):
    #print("\n Iteration :",it)

    Uold=np.array(our_mesh.Uold)  
    U=np.array(our_mesh.vector_U())
    our_mesh.t+=our_mesh.dt
    
for i in range(0,np.size(U[cl,:])):
    print('U_{}({})={},'.format(cl+1,i, U[cl,i]))