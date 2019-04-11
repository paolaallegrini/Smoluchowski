# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 17:08:31 2019

@author: Home
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 09:50:06 2019

@author: Home
"""
from read_file import read_file
from Laplace import FE_method
import numpy as np
from scipy.special import erf
from math import sqrt
from paraview import write_file
    
def vecteurX(Nodes):
    X=[]
    for e in Nodes :
        X.append(e.x)
    return np.array(X)



'Mesh creation from msh file'
our_mesh = read_file("C:/Users/Home/Desktop/stage_labo/Smoluchowski/maillage/square_simple.msh")
#erase_files()
solve=FE_method(our_mesh)


'''parameters'''
L=1000
coeff_d=1
X=vecteurX(our_mesh.Nodes)
#U0=sol_exacte(X,t,coeff_d)

''' Initial situation '''
solve.init_cond(coeff_d)
solve.maj_matrices()

U=solve.vector_U()

#'''Print U'''
#for it in range(0,np.size(solve.U)):
#    print('U({})={}, Uexacte={}'.format(it, solve.U[it],sol_exacte(our_mesh.Nodes[it].x,coeff_d)))

write_file(our_mesh,U,"INIT")
#'L2 error'
#h=L/64
#X=vecteurX(our_mesh.Nodes)
#Uexact=sol_exacte(X,coeff_d)
#err=sum((U-Uexact)**2)*h
#print("Error =",err)