%% Laplace equation de freefem m=1*factor uex= x^3 -2*x*y^2 %%
h=[0.2 0.1 0.05 0.025 0.0125]
err=[0.01 0.0025 0.000625 0.00015625 3.90625e-5]

display('Laplace')
n=4;
for i=1:n-1
    PL2inf(i)=log2(err(i)/err(i+1));
end

PL2inf




%% CN heat equation de freefem m=1*factor, T=3, dt=h %%
h=[0.2 0.1 0.05 0.025 0.0125]
eL2=[0.2933731, 0.0496353, 0.0114092, 0.00279499, 0.000695239]
einfL2=[0.390066, 0.101606, 0.0208001, 0.00500231, 0.00125064]

n=5
for i=1:n-1
    PL2(i)=log2(eL2(i)/eL2(i+1));
    PL2inf(i)=log2(einfL2(i)/einfL2(i+1));
end
display('CN heat equation')
PL2
PL2inf


%% CN heat equation Neuman on border 2 m=80, h=1/m (hmax=0.018), T=3, dt=0.2, 0.1, 0.05, 0.025 %%
eL2=[0.00277706,0.000605764,0.000210822]
for i=1:2
    P(i)=log2(eL2(i)/eL2(i+1))
end
%% dt=h%
e1=0.000982814
e2=0.00021876

e1=0.000246869
e2=6.03889e-005

e1=6.61052e-005
e2=1.73307e-005

p=log2(e1/e2)


%% avec hole
e1=0.146285
e2=0.039237

e1=0.0387469
e2=0.00988092

e1=0.00983485
e2=0.00247218

p=log2(e1/e2)

