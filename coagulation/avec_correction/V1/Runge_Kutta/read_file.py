# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 10:01:43 2019

@author: Home
"""
"""
Reading a gmsh file and filling a mesh object
"""
from base_FE_Q import Mesh, Node, Element, Triangle, Segment
import numpy as np

def read_file(filename):
    Nodes = np.empty((100000, 4), dtype = float)
    Number0fNodes = 0
    NumberOfTr = 0
    NumberOfSeg = 0
    Number0fElems = 0
    NumberOfBorders =0
    Cnt_bord=0
    
    with open(filename) as f:
        content = f.readlines()

    for i in range(0, len(content)):
        line = content[i]

        if line[0] == '$': # reading a property
            
            property = line[1:-1]
            if property =="PhysicalNames":
                i+=1
                line = content[i]
            
                NumberOfBorders = (int)(content[i][0:-1].split(" ")[0]) -1
                #print("Nb bords : {}".format(NumberOfBorders))
                NumberOfSeg = NumberOfBorders
                Cnt_bord=[0]*NumberOfSeg
                i += 1;

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

                    Elems.append(np.asarray([content[j][0:-1].split(" ")[1: ( vertice_n +bntg + 3)]]))  #[cnt] = np.asarray([content[j][0:-1].split(" ")[1:(vertice_n * 2 + 2)]])
                    
                    for c in range(NumberOfSeg) :
                        if Elems[-1][0][0] == '1' and Elems[-1][0][2] == str(c+1):
#                        cnt_1 += 1
                            Cnt_bord[c] +=1
                            break
                    cnt = cnt +1
            else:
                a = 2

    Nodes_ = np.empty(Number0fNodes, dtype = Node) #[Node]* (Number0fNodes+1);
    Elems_ = np.empty(Number0fElems, dtype = Element)
    Trs_ = np.empty(NumberOfTr, dtype = Triangle)
    
    Segs=[np.empty(Cnt_bord[c], dtype = Segment) for c in range(NumberOfSeg)]

    cnt = 0
    for i in range(0, Number0fNodes):
        Nodes_[ i ] = Node(i,Nodes[i][1], Nodes[i][2], Nodes[i][3])


    ide=[0]*(NumberOfSeg)
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
            for c in range(NumberOfSeg) :
                if int(Elems_i[2]) == (c+1): #bord_ext
                    Segs[c][ide[c]] = Segment(i, Elems_i[2+nbTag:]) #seg_s;
                    ide[c] += 1                            
                    break
    princeMesh = Mesh(Number0fNodes, Nodes_, NumberOfTr , Trs_,Segs,Cnt_bord) #def __init__(this, Format_, Ns,Nodes , Nt, Triangles):
    return princeMesh