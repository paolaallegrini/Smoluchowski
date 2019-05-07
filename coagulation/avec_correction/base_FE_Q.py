# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 15:58:58 2019

@author: Home

Smoluchowski project
"""

import numpy as np
from math import sqrt,floor
from scipy.sparse import lil_matrix
#from scipy.sparse.linalg import spsolve


"""
Class Mesh: Contains the essential information of a mesh: Triangles, Nodes, Exterieur Bord and interieur bord
It also aontains the necessery functions that calculate "Masse", "Rigidite" ->A matrices and the vector b 
The function vector_U calculates the solution to the given problem and puts it in the vector U
"""
class Mesh:
    def __init__(this, Ns, Nodes , Nt, Triangles, Segs,Cnt_bord):
        this.Nodes = Nodes # list of nodes : id + coords x,y
        this.Triangles = Triangles # list of triangles id tr + ids vertices
        this.Ns = Ns # nb vertices 
        this.Nt = Nt # nb triangles
        this.grad_phi_ref = np.array([[-1, -1], [1, 0], [0, 1]]) 
        
        #inter and exter borders
        this.Cnt_bord = Cnt_bord # list nb vertices/ border
        this.Bords = Segs 
        this.Bord_1 = Segs[0] # Mur ou Ext
        this.Bord_2 = Segs[1] # Gauche ou Int
        this.Bord_3 = Segs[2] # Droite
        this.Bord_4 = Segs[3] # Cercle
        
        #find border nodes
        this.Nodes_bords=this.find_bords_nodes() # list with list of nodes for each border
        
        # conditions by default
        this.coeff_d = 20
        this.dt = 1
        this.U0=20
        this.t=0
        this.NB=1 # nb of clusters 
        #this.Uold=this.U=this.b=np.zeros([this.NB,this.Ns])
      
    """
    Initial conditions
    """
    def init_cond(this,coeff_d, dt, U0):
        
        this.coeff_d =coeff_d # diffusion coefficient
        this.dt = dt
        this.U0=U0 # initial value for each U 
        NB=np.size(U0) 
        this.NB=NB # number of clusters 
        this.Uold=this.U=this.b=np.zeros([NB,this.Ns])
        
        ' Condition t=0 '
        for m in range(NB):
            this.Uold[m,:]=this.U0[m]*np.ones(this.Ns)
        
#            ' Condition dirichlet '
#            for id_s in this.Nodes_bords[1]:
#                this.Uold[m,id_s-1] = 22
#    
#            for id_s in this.Nodes_bords[2]:
#                this.Uold[m,id_s-1] = 2
        #this.Uold[0,:]=25#this.U0[0]
        
        this.U=np.array(this.Uold) 
        
        return this.Uold, this.U
    
    """ Find equilibrum , steady state
    """
    def equilibrium(this,U,Uold,prec=1e-010):
        norm=np.linalg.norm(U-Uold,np.inf)
        if norm<=prec:
            return True
        return False
            
    """ 
    finds bound's nodes' id_s (for the borders)
    returns : list of list of ids
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
    - quoi indicates the number of the bound, if quoi=1 id belongs to the bound_1
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

    """ Calculates time step : 
        condition dt < U/Qloss
    """

    def dtcalc(this,Qloss):
        dt=dtemp=10
        U=this.Uold
        NB=this.NB
        
        for m in range(NB):       
            if np.all(Qloss[m,:]!=0):
                dtemp=min(U[m,:]/Qloss[m,:] )           
                if (dtemp < dt):
#                    print("m, Qlossm : ",(m,min(Qloss[m,:])))
#                    print("Um : ",min(U[m,:]))
                    dt=round_down(dtemp)
                    
                    
#            print("dt, U[m], dt*Qloss :\n",dt,min(U[m,:]),max(dt*Qloss[m,:]))           
#        print("dt= ",dt)
        this.dt=dt
        return this.dt
    
    """ Coagulation term :
        - Matrix NB*Ns with all the U vectors 
        
    """
    def Qcalc(this,U):
        
        NB=this.NB #Total number of clusters / max size of a cluster
        
        ' Coagulation coefficient (matrix) : a_ij=alpha/(i*j) '
        alpha=10.0
        a=np.array([[alpha/i*j for i in range(1,NB+1)] for j in range(1,NB+1)])

        ' Calculates Qgain and Qloss '
        Qgain=np.zeros((NB,this.Ns),dtype=np.float64)
        Qloss=np.zeros((NB,this.Ns),dtype=np.float64)
        
        ' Qloss for all clusters except for m=NB-1 '
        for m in range(NB-1) : #no Qloss for m=NB-1
            #Qloss
            for id in range(this.Ns):
                Qloss[m,id]=np.dot(U[:,id],a[m,:])
                
                    
            Qloss[m,:]=Qloss[m,:]*U[m,:]
            
            ' Qgain for clusters 1 ... NB-2 '
            if (m!=0) : #no Qgain for m=0 
                for id in range(this.Ns):
                    for j in range(m):
                        ji=j+1 #car pb indice 0 
                        mi=m+1
                        Qgain[m,id]+=a[ji-1,mi-ji-1]*U[ji-1,id]*U[mi-ji-1,id]
            
            
        ' Special case : Qgain for cluster NB-1 '
        for id in range(this.Ns):
            for j in range(NB-1):
                ji=j+1
                for k in range(NB-1):
                    ki=k+1
                    if ((ji+ki)>=NB):
                        Qgain[NB-1,id]+=a[j,k]*U[j,id]*U[k,id]
        ' Q total '
        Q= 1/2.0*Qgain - Qloss
        
        return Q,Qloss

    """
    Matrix A's calculation of cluster m : 
    """
    def matrice_A(this,m=0):
        # this.A = this.M + this.D
        A = lil_matrix((this.Ns, this.Ns))
        
        A = this.M + this.dt*this.coeff_d[m]*this.D
        
        
#        for i in range(0, this.Ns):
#            for j in range(0, this.Ns):
#                this.A[i,j] = this.M[i,j] + this.D[i,j]
                #print ('I ={} J={} M={} D={} A={}'.format(i, j,this.M[i][j],this.D[i][j], this.A[i][j]))
        
        ' condition dirichlet bord'
#        for id_s in this.Nodes_bords[1]: # Gauche
#            this.A[int(id_s) -1,:] = 0
#            this.A[int(id_s) -1,int(id_s) -1] = 1
#        for id_s in this.Nodes_bords[2]: # Droite
#            this.A[int(id_s) -1,:] = 0
#            this.A[int(id_s) -1,int(id_s) -1] = 1

        A=A.tocsr()
        return A
    """
    List of the NB matrices A  (coeff_d changes) : 
    """
    def matrices_Am(this):
        this.Am=[]
        for m in range(this.NB):
            this.Am.append(this.matrice_A(m))
        return this.Am
    
    
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

    """
    Mass matrix :
        - form : int2d_Omega phi_I*phi_J
    """
    def matrice_mass(this):
        this.M = lil_matrix((this.Ns, this.Ns))#, dtype = np.complex)
        
        ' Browses each triangle and adds the elementary contribution on the right node coords '
        for p in range(0, this.Nt): 
            
            'Browses each vertex of triangle p  and calculates integral int2d_p phi_I*phi_J'
            for i in range(0, 3): 
                I = this.Triangles[p].sommets[i]
                for j in range(0, 3):
#                    J = this.Triangles[p].sommets[j]
                    if i == j : # 2
#                        this.M[I-1,J-1] += this.aire_element(p+1)/6.0
                        this.M[I-1,I-1] += this.aire_element(p+1)/6.0
                    else: # 1
#                        this.M[I-1,J-1] += this.aire_element(p+1)/12.0
                        this.M[I-1,I-1] += this.aire_element(p+1)/12.0
                        
        this.M=this.M.tocsr() 
        return this.M

    """
    Rigidity matrix :
            - coeff_d(m) : diffusion coeff of cluster m
            - form :coeff_d(m) * int2d_Omega grad(phi_I)*grad(phi_J)
            
    """
    def matrice_rigidite(this,m=0):

        this.D = lil_matrix((this.Ns, this.Ns))
        ' Browses each triangle and adds the elementary contribution on the right node coords '
        for p in range(0, this.Nt):
            B = this.matrice_B(p)
            bTb = np.dot(np.transpose(B),B)
            
            'Browses each vertex of triangle p  and calculates integral coeff_d * int2d_p grad(phi_I)*grad(phi_J)'
            for i in range(0, 3):
                I = this.Triangles[p].sommets[i]
                for j in range(0, 3):
                    J = this.Triangles[p].sommets[j]
                    this.D[I-1,J-1] += this.aire_element(p+1) * np.dot( np.transpose(this.grad_phi_ref[j]) ,np.dot(bTb, this.grad_phi_ref[i]))
        
        
        this.D=this.D.tocsr()
        return this.D
    
    
    """
    Vector b  :
        b= int2d_Omega (Uold*v) + dt*int2d_Omega (Q(Uold)*v) + Neumann cond
        b= M*(Uold + dt*Q) + Neumann cond
        Neumann cond calculated with quadrature rule 
        
        Returns : Matrix NB*Ns
    """
    
    def vector_b(this,Q):
        #print("U dans b avant Q1:\n",this.Uold[0,:])
        NB=this.NB
        this.b=np.zeros((NB,this.Ns))
        psi=0.5
        
        for m in range(NB): 
            'With coag '
            this.b[m,:]=np.dot(this.M.toarray(),this.dt*Q[m,:] + this.Uold[m,:])
            
#            'Only diffusion'
#            this.b[m,:]+=np.dot(this.M.toarray(),this.Uold[m,:])
            
        'Condition neumann bord int fonction constante /uniquement pour U[0,:]=0.5' 
        for p in range(0,np.size(this.Bord_4)):
            taille=this.aire_seg(p+1,4)
            
            p1=this.Bord_4[p].sommets[0]
            p2=this.Bord_4[p].sommets[1]
            
            this.b[0,p1-1]+=(taille/2)*this.dt*psi
            this.b[0,p2-1]+=(taille/2)*this.dt*psi
            
            
        #print("\nUold=",this.Uold[0,:])
        #print("b0=",this.b[0,:])
        
        return this.b
    
    """ Vector U :
        Solves system of the form AU=b to find U
    """
    
    def vector_U(this):

        Q,Qloss=this.Qcalc(this.Uold)

        this.dt=this.dtcalc(Qloss)
        this.Am=this.matrices_Am()
        
        b=this.vector_b(Q)        
        for m in range(this.NB):

            this.U[m,:] = np.linalg.solve(this.Am[m].toarray(),b[m,:])
            #print ('U({})={},  Q({})={}'.format(m,sum(this.Uold[m,:]),m,sum(Q[m,:])))
            
        this.Uold=np.array(this.U)
        #print("Usortie",this.U[0,:])
        return this.U
        
    """ Creates the matrices M=mass, D=rigidity, A= M + dt*D
    """
    
    def maj_matrices(this):
        this.matrice_mass()
        this.matrice_rigidite()
        this.matrices_Am()
        return
    
"""
    Rounds down the number n with a certain nb of decimals 
    Used to round down the time step dt
"""
def round_down(n, decimals=4):
    multiplier = 10 ** decimals
    return floor(n * multiplier) / multiplier



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