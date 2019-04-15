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
    
    
    "source terme : int2d_Omega f*v"
    def func_f(this,xm,ym) :
        return 0.1*xm + 0.3*ym
    
    def U_D(this,x,y):
        return x*x*x -3*x*y*y
    
    "quadrature"
    def approx_fp(this,xm,ym,i):
        if (i==1): #psi1(1/3)=1
            fp=this.func_f(xm,ym)
        else : #psi1=psi2=1/3
            fp=1/3*this.func_f(xm,ym)
        return fp
    
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
        for id_s in this.mesh.Nodes_bords[0]:
            this.A[int(id_s) -1,:] = 0
            this.A[int(id_s) -1,int(id_s) -1] = 1
            
#        for id_s in this.mesh.Nodes_bords[2]: # Droite
#            this.A[int(id_s) -1,:] = 0
#            this.A[int(id_s) -1,int(id_s) -1] = 1

        this.A=this.A.tocsr()
        return this.A
    
    'source term in all domain'
    
    def vectorbf(this):
        mesh=this.mesh
        bf=np.zeros(mesh.Ns)
        
        for p in range(mesh.Nt) :
            xm,ym=mesh.Tk(p+1,1/3,1/3) #quadrature ordre 1, 1pt
            
            for i in range(0, 3):
                I=this.mesh.Triangles[p].sommets[i]
                bf[I]+=mesh.aire_element(p+1)*this.approx_fp(xm,ym,i)
                
        return bf         
    
    
    def vector_b(this,Ud=10):
        mesh=this.mesh
        this.b=np.zeros(mesh.Ns)

        ' Condition dirichlet '
        for id_s in mesh.Nodes_bords[0]: #bord Gauche
            x=mesh.Nodes[id_s-1].x
            y=mesh.Nodes[id_s-1].y
            this.b[id_s-1] = this.U_D(x,y)

#        for id_s in this.mesh.Nodes_bords[2]: #bord Droit
#            this.b[id_s-1] = 0
        
        'source term f in all omega'
#        this.b+=this.vector_bf()
        
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