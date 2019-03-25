# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 16:34:02 2019

@author: Home
"""

import numpy as np 
from read_file import read_file
from paraview import write_file,erase_files
import numpy as np


'Mesh creation from msh file'
our_mesh = read_file("C:/Users/Home/Desktop/stage_labo/Smoluchowski/maillage/square_4_borders.msh")
erase_files()


'''parameters'''
dt=0.1
coeff_d=2
Ut0=[10.0, 10.0,0]
NB=np.size(Ut0)
Itf=100# Nb iterations

Uold,U=our_mesh.init_cond(coeff_d,dt,Ut0)

''' Initial situation '''
our_mesh.maj_matrices()
print("U:",U)
for it in range(Itf):
    #print("\nIteration :",it)

#    Uold=np.array(our_mesh.Uold)
#    print("\n \t Un=\n",our_mesh.U)
    U=np.array(our_mesh.vector_U())
    our_mesh.t+=dt
#    print("\n \t Un+1=\n",our_mesh.U)
    
    write_file(our_mesh,U[0,:],int(it))
    print("U:\n",U)
