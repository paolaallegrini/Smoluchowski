

Moo=10
Mo=@(t) 2*Moo./(2+Moo*t)

f=@(t,x) Mo(t).*Mo(t).*exp(-Mo(t)*x)

t=1:0.1:3
x=1
f(t,x)
plot(t,f(t,x))