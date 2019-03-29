# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 13:42:02 2019

@author: Home
"""

from read_file_source import read_file
from paraview import write_file,erase_files
import numpy as np

""" Find equilibrum
"""
def equilibrium(Uold,U,prec=1e-04):
    norm=np.linalg.norm(U-Uold,np.inf)
    if norm<=prec:
        return True
    return False

'Mesh creation from msh file'
our_mesh = read_file("C:/Users/Home/Desktop/stage_labo/Smoluchowski/maillage/square_4_borders_hole.msh")
erase_files()
'''parameters'''
dt=1
coeff_d=2
U0=10.0
Itf=10000# Nb iterations
our_mesh.init_cond(coeff_d,dt,U0)

''' Initial situation '''
our_mesh.maj_matrices()

''' Time loop '''
for it in range(Itf):

    if((it%10==0)):
        '''Write solution in paraview format'''
        write_file(our_mesh,int(it/10))
    Uold=our_mesh.Uold
    U=our_mesh.vector_U()
    our_mesh.t+=dt
    
    if equilibrium(Uold,U,prec=1e-5) :
        print('---Equilibrium reached---- : Iteration {} and t={}\n'.format(it,our_mesh.t))
        break;

'''Write Final in paraview format'''
write_file(our_mesh,int(Itf/10))

'''Print Uold'''
for it in range(0,np.size(our_mesh.Uold)):
    print('Uold({})={}'.format(it, our_mesh.Uold[it]))
