# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 15:23:21 2019

@author: Home
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 15:58:58 2019

@author: Home
"""

# -*- coding: utf-8 -*-
"""
Mesh project script containing mesh class
Paola Allegrini
"""

import numpy as np
from math import sqrt
from scipy.sparse import lil_matrix, identity
#from scipy.sparse.linalg import spsolve
#from cmath import exp

"""
Class Mesh: Contains the essential information of a mesh: Triangles, Nodes, Exterieur Bord and interieur bord
It also aontains the necessery functions that calculate "Masse", "Rigidite" ->A matrices and the vector b 
The function vector_U calculates the solution to the given problem and puts it in the vector U
"""
class Mesh:
    def __init__(this, Ns, Nodes , Nt, Triangles, Segs,Cnt_bord):
        this.Nodes = Nodes
        this.Triangles = Triangles
        this.Ns = Ns
        this.Nt = Nt
        this.grad_phi_ref = np.array([[-1, -1], [1, 0], [0, 1]])
        
        #inter and exter borders
        this.Cnt_bord = Cnt_bord
        this.Bords = Segs
        this.Bord_1 = Segs[0] # Mur ou Ext
        this.Bord_2 = Segs[1] # Gauche ou Int
        this.Bord_3 = Segs[2] # Droite
        this.Bord_4 = Segs[3] # Cercle
        
        #find border nodes
        this.Nodes_bords=this.find_bords_nodes()
        
        # conditions by default
        this.coeff_d = 20
        this.dt = 1
        this.U0=20
        this.t=1
        this.NB=1
        #this.Uold=this.U=this.b=np.zeros([this.NB,this.Ns])
      
    """
    conditions
    """
    def init_cond(this,coeff_d, dt, U0):
        this.coeff_d =coeff_d
        this.dt = dt
        this.U0=U0
        NB=np.size(U0)
        this.NB=NB
        this.Uold=this.U=this.b=np.zeros([NB,this.Ns])
        
        ' Condition t=0 '
        for m in range(NB):
            this.Uold[m,:]=this.U0[m]*np.ones(this.Ns)
        
        this.U=np.array(this.Uold)
        
        return this.Uold, this.U

    def maj_matrices(this):
                #matrices
        this.matrice_mass()
        this.matrice_rigidite()
        this.matrice_A()
        return
            
    """
    finds bound's nodes' id_s
    """
    def find_bords_nodes(this):
        this.Nodes_bords =[[],[],[],[]]
        for i in range(0,4):
            for j in range(0, this.Cnt_bord[i]):
                for s in this.Bords[i][j].sommets:
                    if not(s in this.Bords[i]):
                        this.Nodes_bords[i].append(s)

        return this.Nodes_bords


    """
    calculates area of a triangle
    """
    def aire_element(this, id):
        p1 = this.Nodes[this.Triangles[id-1].sommets[0]- 1]
        p2 = this.Nodes[this.Triangles[id-1].sommets[1]- 1]
        p3 = this.Nodes[this.Triangles[id-1].sommets[2]- 1]

        return abs(((p2.x - p1.x)*(p3.y - p1.y) - (p3.x - p1.x)*(p2.y - p1.y)) /2.0)

    """
    calculates area of a line segment , 
    - id : id of segment to find it in segments arrays 
    - quoi = 1 indicates if this segemnts belongs to an exterieur bound / or else it belongs to an internal bound
    """
    def aire_seg(this, id, quoi):
        if quoi == 1:
            p1 = this.Nodes[this.Bord_1[id-1].sommets[0]- 1]
            p2 = this.Nodes[this.Bord_1[id-1].sommets[1]- 1]
        elif quoi==2:
            p1 = this.Nodes[this.Bord_2[id-1].sommets[0]- 1]
            p2 = this.Nodes[this.Bord_2[id-1].sommets[1]- 1]
        elif quoi==3:
            p1 = this.Nodes[this.Bord_3[id-1].sommets[0]- 1]
            p2 = this.Nodes[this.Bord_3[id-1].sommets[1]- 1]
            
        elif quoi==4:
            p1 = this.Nodes[this.Bord_4[id-1].sommets[0]- 1]
            p2 = this.Nodes[this.Bord_4[id-1].sommets[1]- 1]

        return (sqrt( (p1.x - p2.x)**2 + (p1.y - p2.y)**2 ))
    
    """
    Matrix B's calculation, this matrix is to be used in calculating matrix D
    """
    def matrice_B(this, p):

        this.B = np.zeros((2, 2))#, dtype = complex)

        p1 = this.Nodes[this.Triangles[p].sommets[0]- 1]

        p2 = this.Nodes[this.Triangles[p].sommets[1]- 1]

        p3 = this.Nodes[this.Triangles[p].sommets[2]- 1]

        jac = (p2.x - p1.x)*(p3.y - p1.y) - (p3.x - p1.x)*(p2.y - p1.y)

        this.B[0][0] = (p3.y - p1.y)*(1.0/jac)
        this.B[0][1] = (p1.y - p2.y)*(1.0/jac)
        this.B[1][0] = (p1.x - p3.x)*(1.0/jac)
        this.B[1][1] = (p2.x - p1.x)*(1.0/jac)

        return this.B
        
    def matrice_mass(this):
        this.M = lil_matrix((this.Ns, this.Ns))#, dtype = np.complex)

        for p in range(0, this.Nt):
            for i in range(0, 3):
                I = this.Triangles[p].sommets[i]
                for j in range(0, 3):
                    J = this.Triangles[p].sommets[j]
                    if i == j : # 2
                        this.M[I-1,J-1] += this.aire_element(p+1)/6.0
                    else: # 1
                        this.M[I-1,J-1] += this.aire_element(p+1)/12.0
                        
        this.M=this.M.tocsr() 
        return this.M

    def matrice_rigidite(this):

        this.D = lil_matrix((this.Ns, this.Ns))
        for p in range(0, this.Nt):
            B = this.matrice_B(p)
            bTb = np.dot(np.transpose(B),B)

            for i in range(0, 3):
                I = this.Triangles[p].sommets[i]
                for j in range(0, 3):
                    J = this.Triangles[p].sommets[j]
                    this.D[I-1,J-1] += (this.coeff_d*this.dt)*(this.aire_element(p+1) ) * np.dot( np.transpose(this.grad_phi_ref[j]) ,np.dot(bTb, this.grad_phi_ref[i]))
        
        this.D=this.D.tocsr()
        return this.D
    """
    Matrix A's calculation: 
    """
    def matrice_A(this):
        this.A=identity(this.Ns)
        return this.A
    
        """ Coagulation term : 
        - Matrix M*Ns with all the U vectors 
        - n index of the coagulation term we want ( 0->M-1)
    """
    def Qcalc(this,U):
        
        NB=this.NB #nb lignes 
        a=1*np.ones((NB,NB))
        
        Qgain=np.zeros((NB,this.Ns),dtype=np.float64)
        Qloss=np.zeros((NB,this.Ns),dtype=np.float64)
        
        for m in range(NB-1) : #no Qloss for m=NB-1
            #Qloss
            for id in range(this.Ns):
                Qloss[m,id]=np.dot(U[:,id],a[m,:])
                
                    
            Qloss[m,:]=Qloss[m,:]*U[m,:]
            
            #Qgain
            if (m!=0) : #no Qgain for m=0 
                for id in range(this.Ns):
                    for j in range(m):
                        ji=j+1 #car pb indice 0 
                        mi=m+1
                        Qgain[m,id]+=a[ji-1,mi-ji-1]*U[ji-1,id]*U[mi-ji-1,id]
            
        for id in range(this.Ns): #special case m=NB-1
            for j in range(NB-1):
                for k in range(NB-1):
                    if ((j+k)>=(NB-1)):
                        Qgain[NB-1,id]+=a[j,k]*U[j,id]*U[k,id]
        
        Q= 1/2.0*Qgain - Qloss
        return Q,Qloss
    
    
    def dtcalc(this,Qloss):
        dt=dtemp=0.1
        U=this.Uold
        NB=this.NB
        
        for m in range(NB):       
            if np.all(Qloss[m,:]!=0):
                #print("dtemp, m :\n",dtemp,m)
                dtemp=min(U[m,:]/Qloss[m,:] )           
                if (dtemp < dt):
                    #print("dtemp2, m :\n",dtemp,m)
                    dt=dtemp - 0.001
                #print("U[m], dt*Qloss :\n",min(U[m,:]),max(dt*Qloss[m,:]))           
        this.dt=dt
        return this.dt
        
        
    def vector_b(this):
        NB=this.NB
        this.b=np.zeros((NB,this.Ns))
        
        'Coagulation term'
        Q,Qloss=this.Qcalc(this.Uold)
        #print("\n Q:",Q,"\nQloss",Qloss)
        
        'Calc dt'        
        this.dt=this.dtcalc(Qloss)
        print('\n dt :',this.dt)
       
        for m in range(NB): 
            this.b[m,:]=this.dt*Q[m,:] + this.Uold[m,:]
            
        return this.b

    def vector_U(this):
        
        b=this.vector_b()
        for m in range(this.NB):
            this.U[m,:] = np.linalg.solve(this.A.toarray(),b[m,:])
            #this.U[m,this.U[m,:]<0]=0.0   # enlever val negatives
            this.Uold[m,:]=np.array(this.U[m,:])
        
        return this.U



class Node:
  def __init__(this, id, x,y, z):
    this.id = id
    this.x = x
    this.y = y
    this.z = z

class Element:
    def __init__(this, id, typeElem, tags, sommets ):
        this.id = id
        this.typeElem = typeElem
        this.tags = tags
        this.sommets = sommets

class Triangle:
    def __init__(this, id, tags, sommets ):
        this.id = id
        this.tags = tags
        this.sommets = sommets

class Segment:
    def __init__(this, id, sommets):
        this.id = id
        this.sommets = sommets