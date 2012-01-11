clear;
jacob;
clear zzz;
B=Jac;
B(B!=0)=1;
nU = 159;
nP = 53;
nT = nU + nP;
K=Jac(1:nU,1:nU);
G=Jac(1:nU,nU+1:nT);
D=Jac(nU+1:nT,1:nU);
E=Jac(nU+1:nT,nU+1:nT);

rhs;
f=res(1:nU);
g=res(nU+1:nT);
