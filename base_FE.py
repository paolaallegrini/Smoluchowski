# -*- coding: utf-8 -*-
"""
Mesh project script containing mesh class
Mahshid Khezri Nejad Paola Allegrini
Prof: Bertrand Thierry
"""

import numpy as np
from math import cos, sqrt, sin
from cmath import exp


"""
Class Mesh: Contains the essential information of a mesh: Triangles, Nodes, Exterieur Bord and interieur bord
It also aontains the necessery functions that calculate "Masse", "Rigidite" ->A matrices and the vector b 
The function vector_U calculates the solution to the given problem and puts it in the vector U
"""
class Mesh:
    def __init__(this, Format_, Ns,Nodes , Nt, Triangles, b_ext_size, Bord_exts, b_int_size, Bord_int, coeff_d, dt):
        this.Format = Format_
        this.Nodes = Nodes
        this.Triangles = Triangles
        this.Ns = Ns
        this.Nt = Nt
        this.grad_phi_ref = np.array([[-1, -1], [1, 0], [0, 1]])

        this.coeff_d = coeff_d
        this.dt = dt
        
        this.b_ext_size = b_ext_size
        this.Bord_exts = Bord_exts

        this.b_int_size = b_int_size
        this.Bord_int = Bord_int

        this.Nodes_inter = []
        this.find_int_nodes()
        this.matrice_mass()
        this.matrice_rigidite()

        this.vector_b()
        this.matrice_A()
        #this.vector_U()

    """
    finds internal bound's nodes' id_s
    """
    def find_int_nodes(this):
        this.Nodes_inter = []
        
        for i in range(0, this.b_int_size):
            for s in this.Bord_int[i].sommets:
                if not(s in this.Nodes_inter):
                    this.Nodes_inter.append(s)


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
            p1 = this.Nodes[this.Bord_exts[id-1].sommets[0]- 1]
            p2 = this.Nodes[this.Bord_exts[id-1].sommets[1]- 1]
        else:
            p1 = this.Nodes[this.Bord_int[id-1].sommets[0]- 1]
            p2 = this.Nodes[this.Bord_int[id-1].sommets[1]- 1]

        return (sqrt( (p1.x - p2.x)**2 + (p1.y - p2.y)**2 ))

    """
    u_inc function
    - x : first coordinate of a point
    - y : second coordinate of a point
    """
    def u_inc(this, x, y):
        return exp(np.complex(0, 1)*this.k) * (x*cos(this.alpha) + y*sin(this.alpha))

    def vector_b(this):
        this.b = np.zeros(this.Ns, dtype = complex)

        for id_s in this.Nodes_inter:
            p = this.Nodes[id_s-1]
            this.b[id_s-1] = - this.u_inc(p.x, p.y)
        return this.b

    """
    Matrix A's calculation: 
    """
    def matrice_A(this):
        # this.A = this.M + this.D
        this.A = np.zeros((this.Ns, this.Ns), dtype = np.complex)

        for i in range(1, this.Ns):
            for j in range(1, this.Ns):
                this.A[i][j] = this.M[i][j] + this.D[i][j]

#        for id_s in this.Nodes_inter: #no internal border
#            this.A[int(id_s) -1][:] = 0
#            this.A[int(id_s) -1][int(id_s) -1] = 1
        return this.A
    
    """
    Matrix B's calculation, this matrix is to be used in calculating matrix D
    """
    def matrice_B(this, p):

        this.B = np.zeros((2, 2), dtype = complex)

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
        this.M = np.zeros((this.Ns, this.Ns))#, dtype = np.complex)

        for p in range(0, this.Nt):
            for i in range(0, 3):
                I = this.Triangles[p].sommets[i]
                for j in range(0, 3):
                    J = this.Triangles[p].sommets[j]
                    if i == j : # 2
                        this.M[I-1][J-1] += this.aire_element(p+1)/6.0
                    else: # 1
                        this.M[I-1][J-1] += this.aire_element(p+1)/12.0
                        
# conditon ext border
#        for p in range(0, this.b_ext_size):
#            for i in range(0, 2):
#                I = this.Bord_exts[p].sommets[i] 
#                for j in range(0, 2):
#                    J = this.Bord_exts[p].sommets[j]
#                
#                    if i == j : # 2 
##                       this.M[I-1][J-1] += np.complex(0, -1) * (this.k) * this.aire_seg( p+1, 1 ) /3.0 # * (this.k) * this.aire_element(p)/6.0
#                        this.M[I-1][J-1] += 20 * this.aire_seg( p+1, 1 ) /3.0
#                    else : # 1
##                       this.M[I-1][J-1] += np.complex(0, -1) * (this.k) * this.aire_seg( p+1, 1 )/6.0
#                        this.M[I-1][J-1] += 20 * this.aire_seg( p+1, 1 )/6.0

        return this.M

    def matrice_rigidite(this ):
        this.D = np.zeros((this.Ns, this.Ns))#, dtype = np.complex)

        for p in range(0, this.Nt):
            B = this.matrice_B(p)
            bTb = np.matmul(B, np.transpose(B))

            for i in range(0, 3):
                I = this.Triangles[p].sommets[i]
                for j in range(0, 3):
                    J = this.Triangles[p].sommets[j]
                    this.D[I-1][J-1] += -(this.coeff_d*this.dt)*(this.aire_element(p+1) ) * np.matmul( np.transpose(this.grad_phi_ref[j]) ,np.matmul(bTb, this.grad_phi_ref[i]))
        return this.D

    def vector_U(this):
        this.U = np.linalg.solve(this.A, this.b)
        print(this.U)
        for n in this.Nodes:
            this.U[n.id -1 ] = np.abs(this.U[n.id-1] + this.u_inc(n.x, n.y))
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

