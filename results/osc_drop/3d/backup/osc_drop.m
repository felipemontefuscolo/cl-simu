clear;

mu = 0.001;
sig = 72*0.001;
rho = 1000;
R = 2.5*1e-6;
D=2*R;

% = Ampĺitude/Raio
AR = 0.02;

Ohnesorge = mu/sqrt(rho*sig*R)

U = sqrt( 8*sig*0.02^2/(R*rho) )

Weber = rho*D*U/sig

Re = rho*U*D/mu

tr = 2*pi*sqrt(rho*R^3/(8*sig))
