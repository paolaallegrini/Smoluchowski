# -*- coding: utf-8 -*-
"""
Created on Tue May  7 14:52:45 2019

@author: Home
"""
import numpy as np 
from read_file2 import read_file
from matrices_FE_Q import FE_method
from paraview import write_file,erase_files

'Mesh creation from msh file'
our_mesh = read_file("C:/Users/Studente/Documents/GitHub/Smoluchowski/maillage/holes.msh")
erase_files()
solve=FE_method(our_mesh)


'''parameters'''
dt=0.1
NB=2
coeff_d=[ 1.0/(i**(1./3.)) for i in range(1,NB+1)]
#coeff_d=[ 0.1 for i in range(1,NB+1)]
Ut0=np.zeros(NB)
Tf=0.7
Itf=4#int(Tf/dt) # Nb iterations
UM=[]
Uold,U=solve.init_cond(coeff_d,dt,Ut0)

''' Initial situation '''
solve.maj_matrices()
cl=0
it = 0
Ufin = []
while (solve.t < Tf) :
    #print("\n Iteration :",it)

    Uold=np.array(solve.Uold)
    
#    '''Write solution in paraview format'''
#    if((it%10==0)):
#        write_file(our_mesh,Uold[cl,:],int(it/10))
#        
    U=np.array(solve.vector_U(method=1))
    print('----------It {} dt = {:.3f}--------'.format(it,solve.dt))
    print('U_{}({}) = {:.4f}'.format(cl+1,0, U[cl,0]))
    Ufin.append(U[cl,0])
    print ("min value in U :{:.4f}".format(U.min()))


#    Utot=[sum(col) for col in zip(*U)]
    solve.t+=solve.dt
    it += 1

#write_file(our_mesh,Uold[cl,:],int(Itf))

print("Itf : {} Tf :{:.3f}".format(it,solve.t))
#for it in range(0,np.size(U[cl,:])):
##    print('Uold({})={}, Uexacte={}'.format(it, solve.U[it],sol_exacte(our_mesh.Nodes[it].x,solve.t,coeff_d)))
##    print('Uold({})={}'.format(it, solve.Uold[it]))
#    print('U_{}({})={}'.format(cl+1,it, U[cl,it]))

#'L2 error'
##X=vecteurX(our_mesh.Nodes)
##Uexact=sol_exacte(X,solve.t,coeff_d)
#Uexact=sol_exacte(X,Y,solve.t)
#err=sqrt(sum((U-Uexact)**2)*h)
#print("Error =",err)
