# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 09:50:06 2019

@author: Home
"""
from read_file import read_file
from paraview import write_file,erase_files
import numpy as np
from scipy.special import erf
from math import sqrt

""" Find equilibrum
"""
def equilibrium(Uold,U,prec=1e-04):
    norm=np.linalg.norm(U-Uold,np.inf)
    if norm<=prec:
        return True
    return False

""" Sol exacte
"""
def sol_exacte(x,t,coeff_d,Tl=10):
    sol=(1 -erf(x/(2*sqrt(coeff_d*t))))*Tl
    return sol
    
def vecteurX(Nodes):
    X=[]
    for e in Nodes :
        X.append(e.x)
    return np.array(X)

def U_init(Tleft,Tright,Nodes,L=100):
    U0=np.zeros(np.size(Nodes))
    for e in Nodes :
        if (e.x < 0):
            U0[e.id-1]=Tleft
        else:
            U0[e.id-1]=Tright
   
    return U0


'Mesh creation from msh file'
our_mesh = read_file("C:/Users/Home/Desktop/stage_labo/Smoluchowski/maillage/square_simple.msh")
erase_files()
'''parameters'''
L=1000
dt=1/8
Itf=3000*8 # Nb iterations
coeff_d=1
#U0=U_init(5,15,our_mesh.Nodes)
t=200
X=vecteurX(our_mesh.Nodes)
U0=sol_exacte(X,t,coeff_d)
print(U0)
our_mesh.init_cond(coeff_d,dt,U0)
our_mesh.t=t
''' Initial situation '''
our_mesh.maj_matrices()

write_file(our_mesh,"init")

''' Time loop '''
for it in range(Itf):

#    if((it%10==0)):
#        '''Write solution in paraview format'''
#        write_file(our_mesh,int(it/10))
    Uold=our_mesh.Uold
    U=our_mesh.vector_U()
    our_mesh.t+=dt
    
#    if equilibrium(Uold,U,prec=1e-5) :
#        print('---Equilibrium reached---- : Iteration {} and t={}\n'.format(it,our_mesh.t))
#        break;
#
#'''Write Final in paraview format'''
#write_file(our_mesh,int(Itf/10))
#
'''Print Uold'''
for it in range(0,np.size(our_mesh.Uold)):
    print('Uold({})={}, Uexacte={}'.format(it, our_mesh.U[it],sol_exacte(our_mesh.Nodes[it].x,our_mesh.t,coeff_d)))
    #print('Node({})={}'.format(it,our_mesh.Nodes[it].x))


'L2 error'
h=L/64
X=vecteurX(our_mesh.Nodes)
Uexact=sol_exacte(X,our_mesh.t,coeff_d)
err=sum((U-Uexact)**2)*h
print("Error =",err)


