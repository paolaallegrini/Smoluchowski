# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 15:58:03 2019

@author: Home
"""

#from base_FE2 import Mesh, Node, Element, Triangle,Segment
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import numpy as np


class FE_method :
    def __init__(this,mesh):
        this.mesh=mesh
        
        #default :
        this.coeff_d=1
        this.t=0
        this.dt=0
        
    """
    conditions
    """
    def init_cond(this,coeff_d, dt, U0):
        this.coeff_d =coeff_d
        this.dt = dt
        this.U0=U0
        this.Uold=this.U0  #*np.ones(this.Ns)
        
        ' Condition dirichlet '
        for id_s in this.mesh.Nodes_bords[1]:
            this.Uold[id_s-1] = 10

        for id_s in this.mesh.Nodes_bords[2]:
            this.Uold[id_s-1] = 0
    
    
    def matrice_mass(this):
        this.M = lil_matrix((this.mesh.Ns, this.mesh.Ns))#, dtype = np.complex)
    
        for p in range(0, this.mesh.Nt):
            for i in range(0, 3):
                I = this.mesh.Triangles[p].sommets[i]
                for j in range(0, 3):
    #                    J = this.Triangles[p].sommets[j]
                    if i == j : # 2
                        this.M[I-1,I-1] += this.mesh.aire_element(p+1)/6.0
                    else: # 1
                        this.M[I-1,I-1] += this.mesh.aire_element(p+1)/12.0
    
        this.M=this.M.tocsr() 
        return this.M

    def matrice_rigidite(this):
        this.D = lil_matrix((this.mesh.Ns, this.mesh.Ns))
        
        for p in range(0, this.mesh.Nt):
            B = this.mesh.matrice_B(p)
            bTb = np.dot(np.transpose(B),B)

            for i in range(0, 3):
                I = this.mesh.Triangles[p].sommets[i]
                for j in range(0, 3):
                    J = this.mesh.Triangles[p].sommets[j]
                    this.D[I-1,J-1] += (this.mesh.aire_element(p+1) ) * np.dot( np.transpose(this.mesh.grad_phi_ref[j]) ,np.dot(bTb, this.mesh.grad_phi_ref[i]))
        
        this.D=this.D.tocsr()
        return this.D
    

    """
    Matrix A's calculation: 
    """
    def matrice_A(this):
        # this.A = this.M + this.D
        this.A = lil_matrix((this.mesh.Ns, this.mesh.Ns))#, dtype = np.complex)

        this.A= this.M + this.coeff_d*this.dt*this.D
        #print ('I ={} J={} M={} D={} A={}'.format(i, j,this.M[i][j],this.D[i][j], this.A[i][j]))
        
        ' condition dirichlet bord'
        for id_s in this.mesh.Nodes_bords[1]: # Gauche
            this.A[int(id_s) -1,:] = 0
            this.A[int(id_s) -1,int(id_s) -1] = 1
            
        for id_s in this.mesh.Nodes_bords[2]: # Droite
            this.A[int(id_s) -1,:] = 0
            this.A[int(id_s) -1,int(id_s) -1] = 1

        this.A=this.A.tocsr()
        return this.A
    
    
    def vector_b(this):
        
        this.b=np.dot(this.M.toarray(),this.Uold)

        ' Condition dirichlet '
        for id_s in this.mesh.Nodes_bords[1]: #bord Gauche
            this.b[id_s-1] = 10

        for id_s in this.mesh.Nodes_bords[2]: #bord Droit
            this.b[id_s-1] = 0

        ' Condition neumann bord int fonction constante'
#        for p in range(0,this.b_int_size):
#            taille=this.aire_seg(p+1,2)
#            p1=this.Bord_exts[p].sommets[0]
#            p2=this.Bord_exts[p].sommets[1]
#            this.b[p1]+=taille*this.u_inc()
#            this.b[p2]+=taille*this.u_inc()
        
        return this.b

    def vector_U(this):
        
        this.vector_b()
        this.U = spsolve(this.A, this.b)
        this.Uold=this.U
        return this.U
    
    def maj_matrices(this):
        this.matrice_mass()
        this.matrice_rigidite()
        this.matrice_A()
        this.vector_b()
        this.vector_U()
        return