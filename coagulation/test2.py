# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 15:22:13 2019

@author: Home
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 16:34:02 2019

@author: Home
"""

import numpy as np 
from read_file2 import read_file
from paraview import write_file,erase_files
#import numpy as np
from main_coag import equilibrium


'Mesh creation from msh file'
our_mesh = read_file("C:/Users/Home/Desktop/stage_labo/Smoluchowski/maillage/square_holes.msh")
erase_files()


'''parameters'''
dt=0.001
coeff_d=1
Ut0=[10.0, 0.0,0.0,0.0,0.0]
NB=np.size(Ut0)
Itf=1000# Nb iterations

Uold,U=our_mesh.init_cond(coeff_d,dt,Ut0)

''' Initial situation '''
our_mesh.maj_matrices()
Z1=np.array(range(1,6))


for it in range(Itf):
    
    if (it%100==0):
        print("\nIteration :",it)

    Uold=np.array(our_mesh.Uold)
    U=np.array(our_mesh.vector_U(Rf=100))
    our_mesh.t+=dt

    'NB particles + save paraview'
    Utot=[sum(col) for col in zip(*U)]
#    Utot2=np.dot(Z1,U)
#    write_file(our_mesh,Utot2,int(it)) #MASS
     write_file(our_mesh,Utot,int(it))
    
    if  equilibrium(Uold,U,prec=1e-4):
        print("U:\n",U)
        print('---Equilibrium reached---- : Iteration {} and t={}\n'.format(it,our_mesh.t))
        break;



