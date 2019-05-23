# -*- coding: utf-8 -*-
"""
Created on Tue May  7 12:56:42 2019

@author: Home
"""

from scipy.sparse import lil_matrix
#from scipy.sparse.linalg import spsolve
import numpy as np
from math import floor, sqrt


class FE_method :
    "------------------------------------------------------------------------------------------------------------------------------------------------------"
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
                
        this.U = np.array(this.Uold) 

        return this.Uold, this.U
    
    "------------------------------------------------------------------------------------------------------------------------------------------------------"
    """ Find equilibrum , steady state
    """
    def equilibrium(this,U,Uold,prec=1e-010):
        norm=np.linalg.norm(U-Uold,np.inf)
        if norm<=prec:
            return True
        return False
    
    "------------------------------------------------------------------------------------------------------------------------------------------------------"
    """ Calculates time step : 
        condition dt < U/Qloss
    """

    def dtcalc(this,Qloss,dtmax = 0.2):
        dt = dtemp = dtmax
        U = this.Uold
        NB = this.NB
        
        for m in range( NB-1 ):       
            if np.all(Qloss[m,:] != 0):
                dtemp = min( U[m,:] / Qloss[m,:] )   
                
                if (dtemp < dt):
                    dt = round_down(dtemp)
                    
        return dt
    
    "------------------------------------------------------------------------------------------------------------------------------------------------------"
    """ dt calc for stage = 3
        dt < min_{m,id} (U1 / delta*Qloss1[m,id] + (1-delta)*Qloss2[m,id])
        
    """
    def dtcalc3(this,Qloss1,Qloss2,dtmax = 0.2):
#    """
#    Variante avec Q au lieu de Qloss 
#    """
#    def dtcalc3(this,Q1,Q2,dtmax = 0.2):

        dt = dtemp = dtmax
        U1 = this.Uold
        NB = this.NB
        Ns = this.mesh.Ns
        bool = False

        gamma = 1 - sqrt(2)/2.0
        delta = 1 - 1/(2*sqrt(gamma))
        print("In calc dt3")
        for m in range(NB-1) :
            for id in range(Ns) :
                if (U1[m,id]>0) : 
                    #"IN U1 if"
                    denom = delta*Qloss1[m,id] + (1-delta)*Qloss2[m,id]
#                    denom = delta*Q1[m,id] + (1-delta)*Q2[m,id]

                    if (denom > 0) :
                        dtemp = U1[m,id]/denom
                        if (dtemp < dt) :
                            dt = round_down(dtemp) 
                            bool = True
                            print("in")
        return dt,bool
    "------------------------------------------------------------------------------------------------------------------------------------------------------"
    def dtcalc2(this,U,stage=0,U2=None,dtmax = 0.2 ):
        dt = dtmax
        NB = this.NB
        Ns = this.mesh.Ns
        ' Coagulation coefficient (matrix) : a_ij=alpha/(i*j) '
        alpha = 10.0
        a = np.array([[alpha/i*j for i in range(1,NB+1)] for j in range(1,NB+1)])
        
        if (stage == 0) :
            for m in range(NB-1):       
                for id in range(Ns):
                    if ( U[m,id] > 0 ):
                        som = np.dot(U[:,id],a[m,:]) #som = sum(U[:,i]*a[m,:])
                        
                        if ( dt*som >= 1 ):
                            dt = round_down(1/som)
                            
        if (stage == 2) :
            gamma = 1 - sqrt(2)/2
            dt = dtmax   # dt2 = gamma*dt < dt3 or dtmax

            for m in range(NB-1):       
                for id in range(Ns):
                    if (U[m,id]>0):
                        som = np.dot(U[:,id],a[m,:]) # som = sum(U[:,i]*a[m,:])
                        
                        if ( gamma*dt*som >= 1 ):
                            print("IN 2 ")
                            dt = round_down( 1 / (gamma*som))
        return dt
    "------------------------------------------------------------------------------------------------------------------------------------------------------"
    """ Coagulation term :
        - Matrix NB*Ns with all the U vectors 
        
    """
    def Qcalc(this,U):
        mesh=this.mesh
        NB=this.NB #Total number of clusters / max size of a cluster
        
        ' Coagulation coefficient (matrix) : a_ij=alpha/(i*j) '
        alpha=10.0
        a=np.array([[alpha/i*j for i in range(1,NB+1)] for j in range(1,NB+1)])

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
    "------------------------------------------------------------------------------------------------------------------------------------------------------"
    """
    Matrix A's calculation of cluster m : 
        stage=0 -> no runge kutta M+ dt*coeff*D
        stage=2 pu 3 ->  M - gamma*dt*coeff*D
    """
    def matrice_A(this,m=0,stage=0):
        mesh=this.mesh
        gamma= 1 - sqrt(2)/2
        A = lil_matrix((mesh.Ns, mesh.Ns))
        
        
        if (stage==0):
            A = this.M + this.dt*this.coeff_d[m]*this.D
        else :
            A = this.M + gamma*this.dt*this.coeff_d[m]*this.D

        A=A.tocsr()
        return A
    """
    List of the NB matrices A  (coeff_d changes) : 
    """
    def matrices_Am(this,stage=0):
        this.Am=[]
        for m in range(this.NB):
            this.Am.append(this.matrice_A(m,stage=stage))
        return this.Am
    
    "------------------------------------------------------------------------------------------------------------------------------------------------------"
    """
    Mass matrix :   
        - form : int2d_Omega phi_I*phi_J
    'lumped mass matrix'
    """
    def matrice_mass(this): 
        mesh=this.mesh
        this.M = lil_matrix((mesh.Ns, mesh.Ns))
        
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
                        
        this.M = this.M.tocsr() 
        return this.M
    "------------------------------------------------------------------------------------------------------------------------------------------------------"
    """
    Rigidity matrix :
            - coeff_d(m) : diffusion coeff of cluster m
            - form :coeff_d(m) * int2d_Omega grad(phi_I)*grad(phi_J)
            
    """
    def matrice_rigidite(this,m=0):
        mesh = this.mesh
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
    "------------------------------------------------------------------------------------------------------------------------------------------------------"
    """
    func psi different for each neuron
    """
    def func_psi(this,idb):
        return idb+1/this.mesh.NbBords
    "------------------------------------------------------------------------------------------------------------------------------------------------------"   
    """
    Vector b  :
        b = int2d_Omega (Uold*v) + dt*int2d_Omega (Q(Uold)*v) + Neumann cond
        b = M*(Uold + dt*Q) + Neumann cond
        Neumann cond calculated with quadrature rule 
        stage=0 -> No runge kutta : M*(Uold + dt*Q) + Neumann cond
        
        stage=2 -> b = M*(Uold + dt*gamma*Q(Uold)) + Neumann cond
        stage=3 -> b = M*(U2 +  dt*delta*Q(Uold) +  dt*(1-delta)*Q(U2)) - (1-gamma)*dt*coeff_d*D*U2 + Neumann cond
        
        Returns : Matrix NB*Ns
    """
    
    def vector_b(this,Q,stage=0,U2=None):
        NB=this.NB
        mesh=this.mesh
        this.b=np.zeros((NB,mesh.Ns))
#        psi=0.5
        if (stage == 0):
            for m in range(NB): 
                this.b[m,:]=np.dot(this.M.toarray(),this.dt*Q[m,:] + this.Uold[m,:])
                
        elif (stage == 2):
            gamma = 1 - sqrt(2)/2
            for m in range(NB): 
                
                this.b[m,:]=np.dot(this.M.toarray(),gamma*this.dt*Q[m,:] + this.Uold[m,:])
                
        elif (stage == 3):
            gamma = 1 - sqrt(2)/2.0
            delta = 1 - 1/(2*sqrt(gamma))
            
            # U2 passed in arg
            Q2 = this.Qcalc(U2)[0]
            for m in range(NB): 
                this.b[m,:]=np.dot(this.M.toarray(),this.Uold[m,:] + delta*this.dt*Q[m,:] \
                + (1-delta)*this.dt*Q2[m,:]) \
                - np.dot((1 - gamma)*this.dt*this.D.toarray(),U2[m,:])
                
        'Condition neumann bord int fonction constante /uniquement pour U[0,:]=0.5' 
        for idb in range(mesh.NbBords):
            for p in range(0,np.size(mesh.Bords[idb])):
                taille=mesh.aire_seg(p+1,idb)
                
                p1=mesh.Bords[idb][p].sommets[0]
                p2=mesh.Bords[idb][p].sommets[1]
                
#                this.b[0,p1-1]+=(taille/2)*this.dt*this.func_psi(idb)
#                this.b[0,p2-1]+=(taille/2)*this.dt*this.func_psi(idb)
                if (stage == 2):
                    this.b[0,p1-1]+=(taille/2)*this.dt*gamma*0.5#this.func_psi(idb)
                    this.b[0,p2-1]+=(taille/2)*this.dt*gamma*0.5#this.func_psi(idb)
                else:
                    this.b[0,p1-1]+=(taille/2)*this.dt*0.5#this.func_psi(idb)
                    this.b[0,p2-1]+=(taille/2)*this.dt*0.5#this.func_psi(idb)

        
        return this.b
    "------------------------------------------------------------------------------------------------------------------------------------------------------"
    """ Vector U :
        Solves system of the form AU=b to find U
        stage=0 -> no runge kutta 
        stage=2 -> calc U2
        stage=3 -> calc U3=U
    """
    
    def vector_U(this,method=0):

        mesh=this.mesh
        Q,Qloss=this.Qcalc(this.Uold)

        if ((method==0) or (this.t ==0)):
            this.dt = this.dtcalc(Qloss)
            this.Am=this.matrices_Am()  
            b=this.vector_b(Q)
        
        else :
            stage = 2
            dtmax = 0.2
#            i = 1
            while True : #do while
                # dt2
                this.dt = this.dtcalc2(this.Uold,stage=2,dtmax = dtmax)

                'Calculating U2 with curent dt'
                U2=np.zeros([this.NB,mesh.Ns])
                this.Am=this.matrices_Am(stage=stage)
                b=this.vector_b(Q,stage)
                
                for m in range(this.NB):
                    U2[m,:] = np.linalg.solve(this.Am[m].toarray(),b[m,:])
                
                'Calculate dt for stage = 3'
                Qloss2 = this.Qcalc(U2)[1]
#                Q2 = this.Qcalc(U2)[0]

                dt3, bool = this.dtcalc3(Qloss,Qloss2,dtmax = this.dt)
#                dt3, bool = this.dtcalc3(Q,Q2,dtmax = this.dt)
                
                'if < dt2 -> calculate U2 with smaller dt2'# HOW SMALLER ???????
                if bool :  #value changed
                    dtmax = this.dt*0.8 #  0.8 ????
                else : 
                    print('-----Stop-----')
                    print("dt stage 2 :", this.dt)
                    break;

            
            'Calculating b for stage = 3'
            stage = 3
            print("dt stage 3 :",this.dt)
            this.Am=this.matrices_Am(stage=stage)  
            b=this.vector_b(Q,stage=stage,U2=U2)
        
        'Calculating U^n+1'
        for m in range(this.NB):
            this.U[m,:] = np.linalg.solve(this.Am[m].toarray(),b[m,:])
            
        this.Uold=np.array(this.U)
            
        return this.U
    "------------------------------------------------------------------------------------------------------------------------------------------------------"    
    """ Creates the matrices M=mass, D=rigidity, A= M + dt*D
    """
    
    def maj_matrices(this):
        this.matrice_mass()
        this.matrice_rigidite()
        this.matrices_Am()
        return
"------------------------------------------------------------------------------------------------------------------------------------------------------"    
"""
    Rounds down the number n with a certain nb of decimals 
    Used to round down the time step dt
"""
def round_down(n, decimals=4):
    multiplier = 10 ** decimals
    return floor(n * multiplier) / multiplier


################### Archives ##################################
#        if (stage == 3) :
#            gamma = 1 - sqrt(2)/2.0
#            delta = 1 - 1/(2*sqrt(gamma))
#            dt = dtmax / (1 - delta)    # dt3 = (1-delta)*dt < dtmax
#            for m in range(NB-1):
#                for id in range(Ns):
#                    if (U2[m,id]>0):
#                        som = np.dot(U2[:,id],a[m,:])
#                        
#                        if ( (1-delta)*dt*som >= 1 ):
#                            dt = round_down( 1/som ) 
    