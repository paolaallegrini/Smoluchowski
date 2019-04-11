# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 17:04:18 2019

@author: Home
"""

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
    def init_cond(this,coeff_d):
        this.coeff_d =coeff_d


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
        
        #this.D=this.D.tocsr()
        return this.D
    

    """
    Matrix A's calculation: 
    """
    def matrice_A(this):
        #this.A = lil_matrix((this.mesh.Ns, this.mesh.Ns))
        this.A = this.coeff_d*this.D
        
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
        
        this.b=np.zeros(this.mesh.Ns)

        ' Condition dirichlet '
        for id_s in this.mesh.Nodes_bords[1]: #bord Gauche
            this.b[id_s-1] = 10

#        for id_s in this.mesh.Nodes_bords[2]: #bord Droit
#            this.b[id_s-1] = 0
        
        return this.b

    def vector_U(this):
        this.vector_b()
        this.U = spsolve(this.A, this.b)
        return this.U
    
    def maj_matrices(this):
        this.matrice_rigidite()
        this.matrice_A()
        this.vector_b()
        this.vector_U()
        return