# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 10:01:43 2019

@author: Home
"""
"""
Reading a gmsh file and filling a mesh object
"""
from base_FE import Mesh, Node, Element, Triangle, Segment
import numpy as np

def read_file(filename):
    Nodes = np.empty((100000, 4), dtype = float)
    MeshFormat = np.empty(3, dtype = int) #[None] * 3
    Number0fNodes = 0
    NumberOfTr = 0
    NumberOfSeg = 0
    Number0fElems = 0
    coeff_d=0.5
    dt=0.01
    
    cnt_ext = cnt_inter = 0
    with open(filename) as f:
        content = f.readlines()

    for i in range(0, len(content)):
        line = content[i]

        if line[0] == '$': # reading a property
            
            property = line[1:-1]
            if property == "MeshFormat":
                i+= 1
                line = content[i][0:-1]

                formats = np.asarray( line.split(" ") );
                MeshFormat = formats # check if it works 


            elif property == "Nodes":
                i+=1
                line = content[i]
            
                Number0fNodes = (int)(content[i][0:-1].split(" ")[0])

                i += 1;
                
                for j in range(i, i + Number0fNodes):
                    Nodes[j - i] = np.asarray([content[j][0:-1].split(" ")[:4]] )

            elif property == "Elements":
                i+=1
                line = content[i]

                Number0fElems = (int)(content[i][0:-1].split(" ")[0])

                i += 1;
                #lent = content[i][0:-1].split(" ")
                Elems = list() #np.empty((Number0fElems ,0 ), dtype = list)

                cnt = 0;
                for j in range(i, i + Number0fElems): # todo i doesnt change

                    type = (int)(content[j][0:-1].split(" ")[1])
                    bntg = (int)(content[j][0:-1].split(" ")[2])

                    vertice_n = 4;

                    if type == 2: # triangle
                        vertice_n = 3;
                        NumberOfTr += 1
                        
                    elif type == 1:
                        vertice_n = 2;
                        NumberOfSeg +=1

                    Elems.append(np.asarray([content[j][0:-1].split(" ")[1: ( vertice_n +bntg + 3)]]))  #[cnt] = np.asarray([content[j][0:-1].split(" ")[1:(vertice_n * 2 + 2)]])

                    if Elems[-1][0][0] == '1' and Elems[-1][0][2] == '1':
                        cnt_ext += 1
                    if Elems[-1][0][0] == '1' and Elems[-1][0][2] == '2':
                        cnt_inter += 1

                    cnt = cnt +1
            else:
                a = 2

    Nodes_ = np.empty(Number0fNodes, dtype = Node) #[Node]* (Number0fNodes+1);
    Elems_ = np.empty(Number0fElems, dtype = Element)
    Trs_ = np.empty(NumberOfTr, dtype = Triangle)
    

    segs_ext = np.empty(cnt_ext, dtype = Segment)
    segs_int = np.empty(cnt_inter, dtype = Segment)

    cnt = 0
    for i in range(0, Number0fNodes):
        Nodes_[ i ] = Node(i,Nodes[i][1], Nodes[i][2], Nodes[i][3])

    ide_ext = ide_int = 0


    cntT = 0
    for i in range(0, Number0fElems):

        Elems_i = Elems[i][0].astype(int)

        nbTag = int(Elems[i][0][1])
        Elems[i] = Elems[i].astype( int)
        
        Elems_[ i ] = Element(i, Elems_i[0], Elems_i[2:2+nbTag], Elems_i[2+nbTag:] ) #id type tags sommets


        if Elems_i[0] == 2:
            Trs_[cntT] = Triangle(i, Elems_i[2:nbTag],  Elems_i[2+nbTag:])
            cntT +=1


        elif Elems_i[0] == 1:

            if int(Elems_i[2]) == 1: #ext
                segs_ext[ide_ext] = Segment(i, Elems_i[2+nbTag:]) #seg_s;
                ide_ext += 1
                
            else:
                if Elems_i[2] == 2: #int
                    segs_int[ide_int] = Segment(i, Elems_i[2+nbTag:])#seg_s
                    ide_int +=  1

    princeMesh = Mesh(MeshFormat, Number0fNodes, Nodes_, NumberOfTr , Trs_, cnt_ext, segs_ext , cnt_inter, segs_int,coeff_d, dt) #def __init__(this, Format_, Ns,Nodes , Nt, Triangles):
    return princeMesh