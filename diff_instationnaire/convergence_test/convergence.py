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
#print('P=',P)

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
for e in range(np.size(err200)):
    err200[e]=math.sqrt(err200[e])
    errt200[e]=math.sqrt(errt200[e])
    
for e in range(np.size(err200)-1):
    P200.append(math.log1p(err200[e]/err200[e+1])/math.log(2))
    Pdt200.append(math.log1p(errt200[e]/errt200[e+1])/math.log(2))
print('Order of convergence for h=[8,16,32,64]')
print('\n The unsteady diffusion equation (dirichlet condition on one border)')
print('\t P200 : ',P200)
#print('Pdt200=',Pdt200)

"Unsteady diffusion no lump "
ee=[10.874150214047107, 2.167811430499547, 0.6428646694153766, 0.3697269066770792]
et=[10.874150214047107, 2.1658730065527183, 0.6404040801319402, 0.3629706879269922]
P300=[]
Pdt300=[]
for e in range(np.size(ee)-1):
    P300.append(math.log1p(ee[e]/ee[e+1])/math.log(2))
    Pdt300.append(math.log1p(et[e]/et[e+1])/math.log(2))
#print('P300',P300)
print('\t P300',Pdt300)

###########################################################################
'error L2 in time (Implicit Euler), pour h=L/32, Tf=1000'
edt=[3.571498058581578, 3.523259809930193, 3.4993463276403673, 3.487441018399278, 3.481501224439862]

Pdt=[]
for e in range(np.size(edt)):
    edt[e]=math.sqrt(edt[e])
    
for e in range(np.size(edt)-1):
    Pdt.append(math.log1p(edt[e]/edt[e+1])/math.log(2))
print('\nIn time for the unsteady diffusion equation')
print('\t Pdt=',Pdt)

###########################################################################
"Laplace equation :L=10,  U(x,y)=x^3 - 3*x*y^2"

h=[8,16,32,64,128]
errlap=[181.25446246252406, 28.10973969757918, 3.509560464112835, 0.4433517756352217, 0.05779743518189447]
Plap=[]

for e in range(np.size(errlap)):
    errlap[e]=math.sqrt(errlap[e])
    
for e in range(np.size(errlap)-1):  
    Plap.append(math.log1p(errlap[e]/errlap[e+1])/math.log(2))
print('\nLaplace equation')
print('\t P=',Plap)
###########################################################################

""" Diffusion dirichlet
Mesh: rectangle ->  L=1000, l=L/10
coeff_d=1, Tinit=Texact(t=300), Tleft(t)=10
Time : dt=dt=1, Tf=2000
T exact=  (1 -erf(x/(2*sqrt(coeff_d*t))))*Tleft

Cranck-Nicholson scheme
"""

'with h/2'
ec=[16.101430022920887, 3.2568608296069517, 0.9573560043688844, 0.40872309983496746]

'with h/2 and dt/2 aussi'
ect=[16.101430022920887, 3.249946225369745, 0.9461857211184315, 0.3904070397235896]
Pc=[]
Pdtc=[]

for e in range(np.size(ec)-1):
    Pc.append(math.log1p(ec[e]/ec[e+1])/math.log(2))
    Pdtc.append(math.log1p(ect[e]/ect[e+1])/math.log(2))
print('Order of convergence for h=[8,16,32,64]')
print('\n The unsteady diffusion equation (dirichlet condition on one border)')
print('\t Pc : ',Pc)
print('\t Pdtc : ',Pdtc)
#print('Pdt200=',Pdt200)



ec=[(30, 2.09818121768437), (60, 0.7062195766207844), (90, 0.5021214661804264), (120, 0.3626066756110393), (140, 0.3301806218534679)]
ect=[(30, 2.09818121768437), (60, 0.6713310929363228), (90, 0.4727368243282289), (120, 0.335702876035873), (140, 0.31294585587211887)]
dt=[1, 1/2,30/90,90/120,120/140]
Pc=[]
Pdtc=[]
for e in range(4):
    Pc.append(math.log1p(ec[e][1]/ec[e+1][1])/math.log(ec[e+1][0]/ec[e][0]))
    Pdtc.append(math.log1p(ect[e][1]/ect[e+1][1])/math.log(ect[e+1][0]/ec[e][0]))
print('Order of convergence for L/h,  h=[30,60,90,120,140]')
print('\n The unsteady diffusion equation (dirichlet condition on one border)')
print('\t Pc : ',Pc)
print('\t Pdtc : ',Pdtc)




################################################
L=1000
h=[L/30,L/60,L/120,L/240]
ecc=[2.09818121768437, 0.7062195766207844, 0.3626066756110393,0.2619384329466815]
ecct=[2.09818121768437,0.6713310929363228,0.335702876035873]
Pcc=[]
for e in range(np.size(ecct)-1):
    Pcc.append(math.log1p(ecct[e]/ecct[e+1])/math.log(2))

print('Pcc : ',Pcc)
h=[L/7.5,L/15,L/30,L/60,L/120,L/240]

ef=[11.778210526028754,3.939936729398552, 1.461702120700374, 0.5161867362656795, 0.25359384492253234, 0.1093028573150669]
Pf=[]
for e in range(np.size(ef)-1):
    Pf.append(math.log1p(ef[e]/ef[e+1])/math.log(2))

print('Pf : ',Pf)


#################################################
"""Convergence in time Cranck Nicholson
L=1000; h=L/120, Tf=1000, tinit=300
"""
dt=[160,80,40,20,10]
e120=[9.584708769184232, 4.961219678948797, 2.5306367994173864, 1.3572413743625682, 0.7704366605002304]
Pt=[math.log1p(e120[e]/e120[e+1])/math.log(2) for e in range(4)]
print('Pt :', Pt)
#################################################




#e1=[0.5631527342862813, 0.1392650950012512, 0.07922321976162881]
#P=[]
#for e in range(np.size(e1)-1):  
#    P.append(math.log1p(e1[e]/e1[e+1])/math.log(2))
#    
#print("P300=",P)