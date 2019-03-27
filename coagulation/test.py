# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 16:34:02 2019

@author: Home
"""

import numpy as np 
from read_file import read_file
from paraview import write_file,erase_files
import numpy as np
from main_coag import equilibrium


'Mesh creation from msh file'
our_mesh = read_file("C:/Users/Home/Desktop/stage_labo/Smoluchowski/maillage/square_4_borders_hole.msh")
erase_files()


'''parameters'''
dt=0.01
coeff_d=20
Ut0=[10.0, 0.0,0.0,0.0,0.0]
NB=np.size(Ut0)
Itf=500# Nb iterations

Uold,U=our_mesh.init_cond(coeff_d,dt,Ut0)

''' Initial situation '''
our_mesh.maj_matrices()
#print("U:",U)
for it in range(Itf):
    #print("\nIteration :",it)

    Uold=np.array(our_mesh.Uold)
#    print("\n \t Un=\n",our_mesh.U)
    U=np.array(our_mesh.vector_U())
    our_mesh.t+=dt
    #print("\n \t Un+1=\n",our_mesh.U)
    Utot=[sum(col) for col in zip(*U)]
    write_file(our_mesh,Utot,int(it))
    
    if  equilibrium(Uold[0,:],U[0,:],prec=1e-4):
        print("U:\n",U)
        print('---Equilibrium reached---- : Iteration {} and t={}\n'.format(it,our_mesh.t))
        break;
