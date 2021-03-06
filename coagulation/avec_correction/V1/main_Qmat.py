# -*- coding: utf-8 -*-
"""
Created on Tue May  7 14:52:45 2019

@author: Home
"""
import numpy as np 
from read_file import read_file
from matrices_FE_Q import FE_method
from paraview import write_file,erase_files

'Mesh creation from msh file'

our_mesh = read_file("C:/Users/Home/Desktop/stage_labo/Smoluchowski/maillage/square_hole.msh")
erase_files()
solve=FE_method(our_mesh)


'''parameters'''
dt=0.1
NB=3
coeff_d=[ 1.0/(i**(1./3.)) for i in range(1,NB+1)]
#coeff_d=[ 0.1 for i in range(1,NB+1)]
print(coeff_d)
#coeff_d=[ 0.1 for i in range(1,NB+1)]
Ut0=np.zeros(NB)
Tf=100
Itf=4#int(Tf/dt) # Nb iterations
UM=[]
Uold,U=solve.init_cond(coeff_d,dt,Ut0)

''' Initial situation '''
solve.maj_matrices()

for it in range(Itf):
    #print("\n Iteration :",it)

    Uold=np.array(solve.Uold)
    Utot_old=[sum(col) for col in zip(*Uold)]
    
#    if((it%100==0)):
#        '''Write solution in paraview format'''
#        write_file(our_mesh,Uold[NB-1,:],int(it/100))
    
    write_file(our_mesh,Uold[NB-1,:],int(it))
    
    U=np.array(solve.vector_U())
    Utot=[sum(col) for col in zip(*U)]
    
    solve.t+=solve.dt
    #print('UM({})= {}\n'.format(it,sum(U[NB-1,:])))
    UM.append(sum(U[NB-1,:]))
#    if  our_mesh.equilibrium(np.array(U[0,:]),np.array(Uold[0,:]),prec=1e-11):
#        print("U:\n",U[NB-1,:])
#        print('---Equilibrium reached---- : Iteration {} and t={}\n'.format(it,our_mesh.t))
##        break;
#    for it in range(0,np.size(U[0,:])):
#        print('Q_{}({})={},'.format(1,it,our_mesh.Q,))


#write_file(our_mesh,Uold[NB-1,:],int(NB))

cl=0
for it in range(0,np.size(U[cl,:])):
#    print('Uold({})={}, Uexacte={}'.format(it, solve.U[it],sol_exacte(our_mesh.Nodes[it].x,solve.t,coeff_d)))
#    print('Uold({})={}'.format(it, solve.Uold[it]))
    print('U_{}({})={:.6f}'.format(cl+1,it, U[cl,it]))
    
print("Itf : {}, t: {:.2f}".format(Itf,solve.t))

print(U.min())

#'L2 error'
##X=vecteurX(our_mesh.Nodes)
##Uexact=sol_exacte(X,solve.t,coeff_d)
#Uexact=sol_exacte(X,Y,solve.t)
#err=sqrt(sum((U-Uexact)**2)*h)
#print("Error =",err)
