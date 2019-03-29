# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 13:23:05 2019

@author: Home
"""

"""
Reading a gmsh file and filling a mesh object
"""
from base_FE_source import Mesh, Node, Element, Triangle, Segment
import numpy as np

def read_file(filename):
    Nodes = np.empty((100000, 4), dtype = float)
    Number0fNodes = 0
    NumberOfTr = 0
    NumberOfSeg = 0
    Number0fElems = 0
    #NumberOfBorders =0
    cnt_1 = cnt_2 = cnt_3 = tr_4= 0 # compteurs seg dans chaque Bord
    
    with open(filename) as f:
        content = f.readlines()

    for i in range(0, len(content)):
        line = content[i]

        if line[0] == '$': # reading a property
            
            property = line[1:-1]
            if property =="PhysicalNames":
                i+=1
                line = content[i]
            
                #NumberOfBorders = (int)(content[i][0:-1].split(" ")[0]) -1
                #print("nb bords : {}".format(NumberOfBorders))

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
                        NumberOfSeg +=1

                    Elems.append(np.asarray([content[j][0:-1].split(" ")[1: ( vertice_n +bntg + 3)]]))  #[cnt] = np.asarray([content[j][0:-1].split(" ")[1:(vertice_n * 2 + 2)]])

                    if Elems[-1][0][0] == '1' and Elems[-1][0][2] == '1':
                        cnt_1 += 1
                        
                    if Elems[-1][0][0] == '1' and Elems[-1][0][2] == '2':
                        cnt_2 += 1
                        
                    if Elems[-1][0][0] == '1' and Elems[-1][0][2] == '3':
                        cnt_3 += 1
                    
                    if Elems[-1][0][0] == '2' and Elems[-1][0][2] == '4': #triangles source
                        tr_4 += 1
                        
                    cnt = cnt +1
            else:
                a = 2

    Nodes_ = np.empty(Number0fNodes, dtype = Node) #[Node]* (Number0fNodes+1);
    Elems_ = np.empty(Number0fElems, dtype = Element)
    Trs_ = np.empty(NumberOfTr, dtype = Triangle)
    Trs_source =  np.empty(tr_4, dtype = Triangle)

    segs_1 = np.empty(cnt_1, dtype = Segment)
    segs_2 = np.empty(cnt_2, dtype = Segment)
    segs_3 = np.empty(cnt_3, dtype = Segment)

    cnt = 0
    for i in range(0, Number0fNodes):
        Nodes_[ i ] = Node(i,Nodes[i][1], Nodes[i][2], Nodes[i][3])


    ide_1 = ide_2 = ide_3 = idt_4 = 0
    cntT = 0
    for i in range(0, Number0fElems):

        Elems_i = Elems[i][0].astype(int)

        nbTag = int(Elems[i][0][1])
        Elems[i] = Elems[i].astype( int)
        
        Elems_[ i ] = Element(i, Elems_i[0], Elems_i[2:2+nbTag], Elems_i[2+nbTag:] ) #id type tags sommets


        if Elems_i[0] == 2:
            Trs_[cntT] = Triangle(i, Elems_i[2:nbTag],  Elems_i[2+nbTag:])
            cntT +=1
            
            if int(Elems_i[2]) == 4: #Source 
                Trs_source[idt_4] = Triangle(i, Elems_i[2:nbTag],  Elems_i[2+nbTag:]) #Triangle;
                idt_4 += 1


        elif Elems_i[0] == 1:

            if int(Elems_i[2]) == 1: #Mur ou Ext 
                segs_1[ide_1] = Segment(i, Elems_i[2+nbTag:]) #seg_s;
                ide_1 += 1
                
            elif Elems_i[2] == 2: # Interieur ou Gauche
                segs_2[ide_2] = Segment(i, Elems_i[2+nbTag:])#seg_s
                ide_2 +=  1
            elif Elems_i[2] == 3: # Droite
                segs_3[ide_3] = Segment(i, Elems_i[2+nbTag:])#seg_s
                ide_3 +=  1
                
    Segs=[segs_1,segs_2,segs_3]
    Cnt_a=[tr_4,cnt_1,cnt_2,cnt_3]
    princeMesh = Mesh(Number0fNodes, Nodes_, NumberOfTr , Trs_,Trs_source,Segs,Cnt_a) #def __init__(this, Format_, Ns,Nodes , Nt, Triangles):
    return princeMesh