# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 12:19:57 2019

@author: Home
"""
import numpy as np
import math



""" Diffusion dirichlet
Mesh: rectangle ->  L=100, l=L/10
coeff_d=1, Tinit=0, Tleft(t)=10
Time : dt=dt=1, Itf=1000
T exact=  (1 -erf(x/(2*sqrt(coeff_d*t))))*Tleft

"""
'err=max(abs(u_h - u))'
h=[8, 16, 32, 64]
error=[0.5842094801978515, 0.5816900328066532, 0.5673926298383765, 0.5635901759943946]

'error L2'
error1=[3.035734694803761, 2.6055220498150957, 3.571498058581578, 5.540660100148751]

'error L2 avec dt/2 aussi '
error2=[3.035734694803761, 2.5717560247881877, 3.4993463276403673, 5.399063691929157]

P=[]
for e in range(np.size(error)-1):
    P.append(math.log1p(error[e]/error[e+1])/math.log(2))
print('P=',P)

""" Diffusion dirichlet
Mesh: rectangle ->  L=1000, l=L/10
coeff_d=1, Tinit=Texact(t=200), Tleft(t)=10
Time : dt=dt=1, Tf=3000
T exact=  (1 -erf(x/(2*sqrt(coeff_d*t))))*Tleft

"""

'with h/2'
err200=[78.80696081328382, 6.265330296055058, 0.5522186187757253, 0.0600672454362108]

'with h/2 and dt/2 aussi'
errt200=[78.80696081328382, 6.274065660558867, 0.5574413721030831, 0.06157132321084041]

P200=[]
Pdt200=[]
for e in range(np.size(err200)-1):
    P200.append(math.log1p(err200[e]/err200[e+1])/math.log(2))
    Pdt200.append(math.log1p(errt200[e]/errt200[e+1])/math.log(2))
print('P200=',P200)
print('Pdt200=',Pdt200)




'error L2 in time, pour h=L/32, Tf=1000'
edt=[3.571498058581578, 3.523259809930193, 3.4993463276403673, 3.487441018399278, 3.481501224439862]

Pdt=[]
for e in range(np.size(edt)-1):
    Pdt.append(math.log1p(edt[e]/edt[e+1])/math.log(2))
print('Pdt=',Pdt)