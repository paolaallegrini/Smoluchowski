# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 12:23:43 2019

@author: Home
"""
import numpy as np
import math

"space"
# ref sol n=160 // dt=0.1 // Tf=0.8
# apres n=5 : n*2 : n=80 
error = [0.00815798, 0.00201872, 0.000506228, 0.000121633, 2.45617e-5]
PCNs=[math.log(error[e]/error[e+1])/math.log(2) for e in range(np.size(error)-1)]
print(" P in space ", np.around(PCNs,decimals=3))


"time"
# dt=h 
#error=[0.00315156, 0.00230992, 0.00119982, 0.000446297]

#ref sol n=80 // dt=0.0125 // Tf=1.0
#h=L/80 dt=0.2 : dt*0.2 : dt=0.025
error = [0.00922934, 0.00395955, 0.00164619, 0.000540758]
PCN=[math.log(error[e]/error[e+1])/math.log(2) for e in range(np.size(error)-1)]
print(" P exp : ", np.around(PCN,decimals=3))


#ref sol n=160 // dt=0.0125 // Tf=1.0
#h=L/160 dt=0.2 : dt*0.2 : dt=0.025
error = [0.00922993, 0.00396014, 0.00164668, 0.000541066]
PCN=[math.log(error[e]/error[e+1])/math.log(2) for e in range(np.size(error)-1)]
#print(" P exp : ", np.around(PCN,decimals=3))


#ref sol n=80 // dt=0.0125 // Tf=0.8 icard_eps =0.1
#h=L/160 dt=0.2 : dt*0.2 : dt=0.025
error = [0.00484407, 0.00191119, 0.000760457, 0.000244824]
PCN=[math.log(error[e]/error[e+1])/math.log(2) for e in range(np.size(error)-1)]
print(" P picard : ", np.around(PCN,decimals=3))





#Diffusion script diffusion CN 
error = [0.000348569, 0.00158072, 0.00484633, 0.0137782]
PCN=[math.log(error[e+1]/error[e])/math.log(2) for e in range(np.size(error)-1)]
print(" P diff CN : ", np.around(PCN,decimals=3))


error = [0.00688912, 0.00242316, 0.0007903358, 0.000174205]
PCN=[math.log(error[e]/error[e+1])/math.log(2) for e in range(np.size(error)-1)]
print(" P diff CN dt=0.0125 : ", np.around(PCN,decimals=3))



# my heat_eq psi=0.5 Tf=10 n=100  coeffd=0.1   ref sol with dt=0.00625
#dt=0.2:/2:0.0125
error=[0.00349127, 0.00120488, 0.000369723, 5.8835e-5,4.27232e-7, 2.02684e-010]#, 2.00649e-10]
PCN=[math.log(error[e]/error[e+1])/math.log(2) for e in range(np.size(error)-1)]
print(" P diff CN dt=0.00625: ", np.around(PCN,decimals=3))



##ref sol n=80 // dt=0.0125 // Tf=0.8 semi-implicit
##h=L/160 dt=0.2 : dt*0.2 : dt=0.025
#error = [0.011594, 0.00499357, 0.00205587, 0.000671907]
#PCN=[math.log(error[e]/error[e+1])/math.log(2) for e in range(np.size(error)-1)]
#print(" P semi_imp  : ", np.around(PCN,decimals=3))

##avec u_n+1 * u_n
#error = [0.0126218, 0.00558442, 0.00232012, 0.000760052]
#PCN=[math.log(error[e]/error[e+1])/math.log(2) for e in range(np.size(error)-1)]
#print(" P lin : ", np.around(PCN,decimals=3))




#heat eq example freefem
#error=[0.132463,0.0351687, 0.00893047]
## at TF 
##error=[0.0355335,0.00885413,0.0022132]	
## max with lumped
##error=[0.101606,0.0208001 ,0.00500231]
##at Tf with lumped 
##error=[0.0496353,0.0114092,0.00279499]
#PCN=[math.log(error[e]/error[e+1])/math.log(2) for e in range(np.size(error)-1)]
#print(" P diff CN : ", np.around(PCN,decimals=3))
