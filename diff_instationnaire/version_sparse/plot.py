# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 14:11:19 2019

@author: Home
"""

#import numpy as np 
#from base_FE2 import Mesh, Node, Element, Triangle, Segment
#from read_file import read_file

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.cm as cm
import matplotlib.animation as animation

def plot_mesh(mesh,U):
    Ns=mesh.Ns
    Nt=mesh.Nt
    Nodes=mesh.Nodes
    Triangles=mesh.Triangles
    
    X =[]
    for i in range(Ns):
        X.append([Nodes[i].x,Nodes[i].y])
    X=np.asarray(X)
    
    triang=[]
    for i in range(Nt):
        triang.append(Triangles[i].sommets -[1,1,1])

    "Create triangulation"
    trig= mtri.Triangulation(X[:,0],X[:,1], triangles=triang)
    
    "Plot mesh"
    plt.triplot(trig, lw=0.5, color='white')
    
    "plot solution"
    cmap = cm.get_cmap(name='coolwarm', lut=None)
    plt.tricontourf(trig,U,cmap=cmap)

    return

def plot_animation(mesh,Uit) :
    Ns=mesh.Ns
    Nt=mesh.Nt
    Nodes=mesh.Nodes
    Triangles=mesh.Triangles
    
    X =[]
    for i in range(Ns):
        X.append([Nodes[i].x,Nodes[i].y])
    X=np.asarray(X)
    
    triang=[]
    for i in range(Nt):
        triang.append(Triangles[i].sommets -[1,1,1])

    "Create triangulation"
    trig= mtri.Triangulation(X[:,0],X[:,1], triangles=triang)
    
    im_strs = []
    fig_s, ax_s = plt.subplots()
    ax_s.set_aspect('equal')
    ax_s.set_axis_off()
    
    Itf=np.shape(Uit)[0]
    cmap = cm.get_cmap(name='coolwarm', lut=None)
    
    for it in range(Itf):
        # plot new U
        z=Uit[it,:]
        C = ax_s.tricontourf(trig, z, 10, cmap=cmap)
        im_strs.append(C.collections)
    
    
    ani_strs = animation.ArtistAnimation(fig_s, im_strs, interval=150)
    ani_strs.save('Solution_Uit.gif', writer='imagemagick', bitrate=300)
    print('Solution_Uit.gif')
    return
