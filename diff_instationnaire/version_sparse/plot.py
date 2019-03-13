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

def plot_mesh(mesh,Uit):
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

    trig= mtri.Triangulation(X[:,0],X[:,1], triangles=triang)
    #plt.triplot(trig, marker="o")
    #plt.tricontourf(trig, mesh.U)
    
    plt.ion()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    line1, = ax.tricontourf(trig, Uit[0,:])
    
    for it in range(np.size(Uit)):
        line1.set_zdata(Uit[it,:])
        fig.canvas.draw()
    
#    # Set up the figure
#    fig, axs = plt.subplots(nrows=2, ncols=2)
#    axs = axs.flatten()
#    
#    # Plot the triangulation.
#    axs[0].tricontourf(trig, mesh.U)
#    axs[0].triplot(trig, 'ko-')
#    axs[0].set_title('Triangular grid')

    return
