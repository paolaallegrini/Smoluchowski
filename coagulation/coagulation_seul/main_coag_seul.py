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
import sys
sys.path.append('../')
from read_file_coag_seul import read_file
from paraview import write_file,erase_files
from main_coag import equilibrium
from math import exp

'Mesh creation from msh file'
our_mesh = read_file("C:/Users/Home/Desktop/stage_labo/Smoluchowski/maillage/square_4elems.msh")
erase_files()

'Analytic sol for a=1'
def func_f(t,x,Moo=10):
    Mo=(2*Moo)/(2 + Moo*t)
    
    f=Mo*Mo*exp(-Mo*x)
    return f

def d_f(t,x,Moo=10):
    Mo=(2*Moo)/(2 + Moo*t)
    
    df=Mo*Mo*Mo*exp(-Mo*x)*(-1+0.5*x*Mo)
    return df

'''parameters'''
dt=100
coeff_d=0
Moo=10
NB=3
print(NB)
t0=1
Ut0=[func_f(t0,i) for i in range(1,NB+1)]

#NB=np.size(Ut0)


Itf=25# Nb iterations


Uold,U=our_mesh.init_cond(coeff_d,dt,Ut0)


''' Initial situation '''
our_mesh.maj_matrices()


for it in range(Itf):
    
    if (it%100==0):
        print("\nIteration :",it)
        
    Uold=np.array(our_mesh.Uold)
    Q=our_mesh.Qcalc(Uold)[0]
    print("\n t: ", our_mesh.t)
    print("\nUold:",Uold[:,0])
    print("f_0 = {}, f_1 = {},f_2 = {} ".format(func_f(our_mesh.t,1),func_f(our_mesh.t,2),func_f(our_mesh.t,3)))
    
    #print("\nQ =",Q[:,0])
    #print("df_0 = {}, df_1 = {},df_2 = {} ".format(d_f(our_mesh.t,1),d_f(our_mesh.t,2),d_f(our_mesh.t,3)))
    
    U=np.array(our_mesh.vector_U())
    #print("U:\n",U)
    our_mesh.t+=our_mesh.dt

    'Save paraview'
    Utot=[sum(col) for col in zip(*U)]
#    Utot2=np.dot(Z1,U)
#    write_file(our_mesh,Utot2,int(it)) #MASS
    write_file(our_mesh,Uold[0,:],int(it))
    
#    if  equilibrium(Uold[:,0],U[:,0],prec=1e-4):
#        print("U:\n",U)
#        print('---Equilibrium reached---- : Iteration {} and t={}\n'.format(it,our_mesh.t))
#        break;


