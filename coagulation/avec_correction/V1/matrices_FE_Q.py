# -*- coding: utf-8 -*-
"""
Created on Tue May  7 12:56:42 2019

@author: Home
"""

from scipy.sparse import lil_matrix
#from scipy.sparse.linalg import spsolve
import numpy as np
from math import floor


class FE_method :
    def __init__(this,mesh):
        this.mesh=mesh
        
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
        mesh=this.mesh
        this.coeff_d =coeff_d # diffusion coefficient
        this.dt = dt
        this.U0=U0 # initial value for each U 
        NB=np.size(U0) 
        this.NB=NB # number of clusters 
        this.Uold=this.U=this.b=np.zeros([NB,mesh.Ns])
        
        ' Condition t=0 '
        for m in range(NB):
            this.Uold[m,:]=this.U0[m]*np.ones(mesh.Ns)
        
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

    """ Calculates time step : 
        condition dt < U/Qloss
    """

    def dtcalc(this,Qloss):
        dt=dtemp=1
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
        mesh=this.mesh
        NB=this.NB #Total number of clusters / max size of a cluster
        
        ' Coagulation coefficient (matrix) : a_ij=alpha/(i*j) '
        alpha=10.0
        a=np.array([[alpha/(i*j) for i in range(1,NB+1)] for j in range(1,NB+1)])
        #a=np.ones((NB,NB));
        ' Calculates Qgain and Qloss '
        Qgain=np.zeros((NB,mesh.Ns),dtype=np.float64)
        Qloss=np.zeros((NB,mesh.Ns),dtype=np.float64)
        
        ' Qloss for all clusters except for m=NB-1 '
        for m in range(NB-1) : #no Qloss for m=NB-1
            #Qloss
            for id in range(mesh.Ns):
                Qloss[m,id]=np.dot(U[:,id],a[m,:])
                
                    
            Qloss[m,:]=Qloss[m,:]*U[m,:]
            
            ' Qgain for clusters 1 ... NB-2 '
            if (m!=0) : #no Qgain for m=0 
                for id in range(mesh.Ns):
                    for j in range(m):
                        ji=j+1 #car pb indice 0 
                        mi=m+1
                        Qgain[m,id]+=a[ji-1,mi-ji-1]*U[ji-1,id]*U[mi-ji-1,id]
            
            
        ' Special case : Qgain for cluster NB-1 '
        for id in range(mesh.Ns):
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
        mesh=this.mesh
        # this.A = this.M + this.D
        A = lil_matrix((mesh.Ns, mesh.Ns))
        
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
    Mass matrix :
        - form : int2d_Omega phi_I*phi_J
    """
    def matrice_mass(this):
        mesh=this.mesh
        this.M = lil_matrix((mesh.Ns, mesh.Ns))#, dtype = np.complex)
        
        ' Browses each triangle and adds the elementary contribution on the right node coords '
        for p in range(0, mesh.Nt): 
            
            'Browses each vertex of triangle p  and calculates integral int2d_p phi_I*phi_J'
            for i in range(0, 3): 
                I = mesh.Triangles[p].sommets[i]
                for j in range(0, 3):
#                    J = this.Triangles[p].sommets[j]
                    if i == j : # 2
#                        this.M[I-1,J-1] += this.aire_element(p+1)/6.0
                        this.M[I-1,I-1] += mesh.aire_element(p+1)/6.0
                    else: # 1
#                        this.M[I-1,J-1] += this.aire_element(p+1)/12.0
                        this.M[I-1,I-1] += mesh.aire_element(p+1)/12.0
                        
        this.M=this.M.tocsr() 
        return this.M

    """
    Rigidity matrix :
            - coeff_d(m) : diffusion coeff of cluster m
            - form :coeff_d(m) * int2d_Omega grad(phi_I)*grad(phi_J)
            
    """
    def matrice_rigidite(this,m=0):
        mesh=this.mesh
        this.D = lil_matrix((mesh.Ns, mesh.Ns))
        ' Browses each triangle and adds the elementary contribution on the right node coords '
        for p in range(0, mesh.Nt):
            B = mesh.matrice_B(p)
            bTb = np.dot(np.transpose(B),B)
            
            'Browses each vertex of triangle p  and calculates integral coeff_d * int2d_p grad(phi_I)*grad(phi_J)'
            for i in range(0, 3):
                I = mesh.Triangles[p].sommets[i]
                for j in range(0, 3):
                    J = mesh.Triangles[p].sommets[j]
                    this.D[I-1,J-1] += mesh.aire_element(p+1) * np.dot( np.transpose(mesh.grad_phi_ref[j]) ,np.dot(bTb, mesh.grad_phi_ref[i]))
        
        
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
        mesh=this.mesh
        this.b=np.zeros((NB,mesh.Ns))
        psi=0.5
        
        for m in range(NB): 
            'With coag '
            this.b[m,:]=np.dot(this.M.toarray(),this.dt*Q[m,:] + this.Uold[m,:])
            
#            'Only diffusion'
#            this.b[m,:]+=np.dot(this.M.toarray(),this.Uold[m,:])
            
        'Condition neumann bord int fonction constante /uniquement pour U[0,:]=0.5' 
#        for p in range(0,np.size(mesh.Bord_4)):
#            taille=mesh.aire_seg(p+1,4)
#            
#            p1=mesh.Bord_4[p].sommets[0]
#            p2=mesh.Bord_4[p].sommets[1]
        for p in range(0,np.size(mesh.Bord_2)):
            taille=mesh.aire_seg(p+1,2)
            
            p1=mesh.Bord_2[p].sommets[0]
            p2=mesh.Bord_2[p].sommets[1]
            
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
        print("dt", this.dt)
        this.Am=this.matrices_Am()
        
        b=this.vector_b(Q)        
        for m in range(this.NB):

            this.U[m,:] = np.linalg.solve(this.Am[m].toarray(),b[m,:])
            #print ('U({})={},  Q({})={}'.format(m,sum(this.Uold[m,:]),m,sum(Q[m,:])))
            
        this.Uold=np.array(this.U)

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
    