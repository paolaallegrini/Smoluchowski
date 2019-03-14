# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 09:50:06 2019

@author: Home
"""
from read_file import read_file
from paraview import write_file,erase_files
import numpy as np

'Mesh creation from msh file'
our_mesh = read_file("../../../maillage/square_4_borders_hole.msh")
erase_files()
'''parameters'''
dt=2
coeff_d=1
U0=10.0
Itf=100# Nb iterations
our_mesh.init_cond(coeff_d,dt,U0)

''' Initial situation '''
our_mesh.maj_matrices()



''' Time loop '''
#Uit=np.array([])
#Uit=np.concatenate((Uit,our_mesh.U))
for it in range(Itf):
    #print('Iteration : %d'%it)
    if((it%10==0)):
        '''Write solution in paraview format'''
        write_file(our_mesh,int(it/10))

    U=our_mesh.vector_U()
    #Uit=np.concatenate((Uit,U))
    our_mesh.t+=dt

'''Write solution in paraview format'''
write_file(our_mesh,int(Itf/10))

#
#Uit=Uit.reshape((Itf+1,our_mesh.Ns))
##np.savetxt("mat_Uit.csv",Uit,delimiter=",")


'''Print Uold'''
for it in range(0,np.size(our_mesh.Uold)):
    print('Uold({})={}'.format(it, our_mesh.Uold[it]))


'''Save animation '''
#plot_animation(our_mesh,Uit)
#plt.subplot(211)
#plot_mesh(our_mesh,our_mesh.U)
#plt.subplot(212)
#plot_quadgrid(our_mesh,our_mesh.U)



