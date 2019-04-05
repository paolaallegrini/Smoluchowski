# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 16:34:02 2019

@author: Home
"""

import numpy as np 
from read_file import read_file
from paraview import write_file,erase_files

'Mesh creation from msh file'

our_mesh = read_file("C:/Users/Home/Desktop/stage_labo/Smoluchowski/maillage/square_4_borders_hole.msh")
erase_files()


'''parameters'''
dt=0.01
coeff_d=[1/i for i in range(1,17)]
Ut0=np.zeros(16)
NB=np.size(Ut0)
Itf=100# Nb iterations

Uold,U=our_mesh.init_cond(coeff_d,dt,Ut0)

''' Initial situation '''
our_mesh.maj_matrices()

for it in range(Itf):
    #print("\n Iteration :",it)

    Uold=np.array(our_mesh.Uold)
    Utot_old=[sum(col) for col in zip(*Uold)]
    
#    if((it%10==0)):
#        '''Write solution in paraview format'''
#        write_file(our_mesh,Utot_old,int(it/10))
    
    U=np.array(our_mesh.vector_U())
    Utot=[sum(col) for col in zip(*U)]
    
    our_mesh.t+=our_mesh.dt
   
    if  our_mesh.equilibrium(np.array(Utot),np.array(Utot_old),prec=1e-5):
        print("U:\n",U[NB-1,:])
        print('---Equilibrium reached---- : Iteration {} and t={}\n'.format(it,our_mesh.t))
        break;
        
print("UM:\n",U[NB-1,:])
