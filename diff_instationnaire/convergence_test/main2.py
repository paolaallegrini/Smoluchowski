# -*- coding: utf-8 -*-
"""
Created on Fri May  3 11:02:10 2019

@author: Home
"""

from read_file import read_file
#from FE_CN_check import FE_method
from matrices_FE_CN import FE_method
from paraview import write_file
import numpy as np
from scipy.special import erf
from math import sqrt,cos,sin,pi


""" Find equilibrum
"""
def equilibrium(Uold,U,prec=1e-04):
    norm=np.linalg.norm(U-Uold,np.inf)
    if norm<=prec:
        return True
    return False

""" Sol exacte
"""
#def sol_exacte(x,t,coeff_d,Tl=10):
#    sol=(1 -erf(x/(2*sqrt(coeff_d*t))))*Tl
#    return sol

def sol_exacte(X,Y,t):
    sol=[ 5*cos(10*t)*sin(2*pi*X[e])*cos(pi*Y[e]) for e in range(np.size(X))]
    return np.array(sol)
    
def vecteurX(Nodes):
    X=[]
    Y=[]
    for e in Nodes :
        X.append(e.x)
        Y.append(e.x)
    return np.array(X),np.array(Y)

'Mesh creation from msh file'
#our_mesh = read_file("C:/Users/Home/Desktop/stage_labo/Smoluchowski/maillage/mesh_laplace.msh")
our_mesh = read_file("D:/stage_labo/Smoluchowski/maillage/mesh_laplace.msh")
#erase_files()
solve=FE_method(our_mesh)


'''parameters'''
L=1
dt=1
Tf=50
Itf=int(Tf/dt) # Nb iterations
#dt=dt/2
#Itf=Itf*2
coeff_d=0.1
t=0
X,Y=vecteurX(our_mesh.Nodes)
U0=sol_exacte(X,Y,t)
U0=10*np.ones(our_mesh.Ns)
''' Initial situation '''
solve.init_cond(coeff_d,dt,U0)
#solve.t=t
solve.maj_matrices()
#print(solve.vector_b())
#write_file(our_mesh,solve.Uold,"init")

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
#write_file(our_mesh,solve.U,"Fin")
#
Uexact=sol_exacte(X,Y,solve.t)

'''Print Uold'''
#for it in range(0,np.size(solve.Uold)):
#    print('Uold({})={}, Uexacte={}'.format(it, solve.U[it],sol_exacte(our_mesh.Nodes[it].x,solve.t,coeff_d)))
#    print('Uold({})={}'.format(it, solve.Uold[it]))
    #print('Uold({})={}, Uexacte={}'.format(it, solve.U[it],Uexact[it]))
#    print('b({})={}'.format(it,solve.b[it]))

'L2 error'
#h=L/(8)
h=np.size(our_mesh.Nodes)
print(h)
#X=vecteurX(our_mesh.Nodes)
#Uexact=sol_exacte(X,solve.t,coeff_d)
Uexact=sol_exacte(X,Y,solve.t)
err=sqrt(sum((U-Uexact)**2)/h)
print("Error =",err)


