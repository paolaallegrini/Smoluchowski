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
dt=1
#coeff_d=[ i**(1./3.) for i in range(1,17)]
NB=16
coeff_d=[ 1.0/(i**(1./3.)) for i in range(1,NB+1)]
Ut0=np.zeros(NB)

Itf=1000 # Nb iterations
UM=[]
Uold,U=our_mesh.init_cond(coeff_d,dt,Ut0)

''' Initial situation '''
our_mesh.maj_matrices()

for it in range(Itf):
    #print("\n Iteration :",it)

    Uold=np.array(our_mesh.Uold)
    Utot_old=[sum(col) for col in zip(*Uold)]
    
#    if((it%100==0)):
#        '''Write solution in paraview format'''
#        write_file(our_mesh,Uold[NB-1,:],int(it/100))
    
    #write_file(our_mesh,Uold[NB-1,:],int(it))
    
    U=np.array(our_mesh.vector_U())
    Utot=[sum(col) for col in zip(*U)]
    
    our_mesh.t+=our_mesh.dt
    #print('UM({})= {}\n'.format(it,sum(U[NB-1,:])))
    UM.append(sum(U[NB-1,:]))
#    if  our_mesh.equilibrium(np.array(U[0,:]),np.array(Uold[0,:]),prec=1e-11):
#        print("U:\n",U[NB-1,:])
#        print('---Equilibrium reached---- : Iteration {} and t={}\n'.format(it,our_mesh.t))
#        break;
write_file(our_mesh,Uold[0,:],int(11))
write_file(our_mesh,Uold[4,:],int(55))
write_file(our_mesh,Uold[NB-1,:],int(NB))