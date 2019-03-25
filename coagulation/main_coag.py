# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 09:50:06 2019

@author: Home
"""
from read_file import read_file
from paraview import write_file,erase_files
import numpy as np

""" Find equilibrum
"""
def equilibrium(Uold,U,prec=1e-04):
    norm=np.linalg.norm(U-Uold,np.inf)
    if norm<=prec:
        return True
    return False

""" Utot function 
Calculates the sum of all concentrations U
"""
def Utot(U):
    return np.array([sum(col) for col in zip(*U)])



def printU(U,Ns):
    '''Print Uold'''
    for i in range(0,Ns):
        print('U({})={}'.format(i,U[i]))


'Mesh creation from msh file'
our_mesh = read_file("C:/Users/Home/Desktop/stage_labo/Smoluchowski/maillage/square_4_borders_hole.msh")
erase_files()


'''parameters'''
dt=0.01
coeff_d=1
Ut0=[10.0, 0.0,0.0,0.0,0.0,0.0]
NB=np.size(Ut0)
Itf=5# Nb iterations

Uold,U=our_mesh.init_cond(coeff_d,dt,Ut0)


''' Initial situation '''
our_mesh.maj_matrices()
#Utotal2=Utot(our_mesh.U)

#printU(Uold[0,:],our_mesh.Ns)

''' Time loop '''
for it in range(Itf):
    #print("Iteration : ",it)
#    if((it%10==0)):
#        #Utotal1=Utotal2
#        '''Write solution in paraview format'''
#        write_file(our_mesh,U[0,:],int(it/10))
    #print(U[0,:])
    write_file(our_mesh,Utot(U),int(it))
    #print(U[0,:])
    
    Uold=np.array(our_mesh.Uold)
    U=np.array(our_mesh.vector_U())
    our_mesh.t+=dt

    if  equilibrium(Uold[0,:],U[0,:],prec=1e-5):
        print('---Equilibrium reached---- : Iteration {} and t={}\n'.format(it,our_mesh.t))
        break;


#Utotal=Utot(our_mesh.U)
'''Write Final in paraview format'''
#write_file(our_mesh,our_mesh.U[0,:],int(Itf/10))


#'''Print Utot'''
#for i in range(0,our_mesh.Ns):
#    print('Uold({})={}'.format(i,Uold[1,i]))

'''Save animation '''
#plot_animation(our_mesh,Uit)
#plt.subplot(211)
#plot_mesh(our_mesh,our_mesh.U)
#plt.subplot(212)
#plot_quadgrid(our_mesh,our_mesh.U)


