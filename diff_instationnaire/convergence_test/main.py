# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 09:50:06 2019

@author: Home
"""
from read_file import read_file
from matrices_FE_CN import FE_method
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

'Mesh creation from msh file'
our_mesh = read_file("C:/Users/Home/Desktop/stage_labo/Smoluchowski/maillage/square_4_borders16.msh")
#erase_files()
solve=FE_method(our_mesh)


'''parameters'''
L=1000
dt=1
Tf=2000
Itf=int(Tf/dt) # Nb iterations
#dt=dt/2
#Itf=Itf*2
coeff_d=1
t=300
X=vecteurX(our_mesh.Nodes)
U0=sol_exacte(X,t,coeff_d)

''' Initial situation '''
solve.init_cond(coeff_d,dt,U0)
solve.t=t
solve.maj_matrices()

#write_file(our_mesh,"init")

''' Time loop '''
for it in range(Itf):

#    if((it%10==0)):
#        '''Write solution in paraview format'''
#        write_file(our_mesh,int(it/10))
    Uold=solve.Uold
    U=solve.vector_U()
    solve.t+=dt
    
#    if equilibrium(Uold,U,prec=1e-5) :
#        print('---Equilibrium reached---- : Iteration {} and t={}\n'.format(it,our_mesh.t))
#        break;
#
#'''Write Final in paraview format'''
#write_file(our_mesh,int(Itf/10))
#
'''Print Uold'''
#for it in range(0,np.size(solve.Uold)):
#    print('Uold({})={}, Uexacte={}'.format(it, solve.U[it],sol_exacte(our_mesh.Nodes[it].x,solve.t,coeff_d)))

'L2 error'
h=L/(15)
X=vecteurX(our_mesh.Nodes)
Uexact=sol_exacte(X,solve.t,coeff_d)
err=sqrt(sum((U-Uexact)**2)*h)
print("Error =",err)


