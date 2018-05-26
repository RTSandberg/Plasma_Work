function atotvec = odef_tracer(x,ode_params) %xtracer,f0vec,c1,c2,L)
%%% right-hand side of Vlasov 1dES Lagrangian particle method ODEs 
%%% account for tracer with a parameter Ntr
%%% input:
%%%     x:  (2*Ns+2*Ntr,1) vector, contains position,velocity of source and
%%%         tracer in form [xsource; vsource; xtracer; vtracer]
%%%     ode_params:     struct, contains parameter list
%%%         fields:     
%%%                     function: string, indicates name of file with right-hand side of ODE
%%%                     Ntr : scalar, number of tracers
%%%                     parameters specific to Lagrangian:
%%%                     f0vec:  vector, contains f0 weights
%%%                     c1:     scalar, needed for 1dES Lagrangian ode function
%%%                     c2:     scalar, "" 
%%%                     L:      scalar, indicates size of computational domain
f0vec = ode_params.f0vec; 
c1 = ode_params.c1; 
c2 = ode_params.c2; 
L = ode_params.L;

Ns = length(f0vec);
Ntr = ode_params.Ntr;

 

xsvec = x(1:Ns); vsvec = x(Ns+1:2*Ns);
xtrvec = x(2*Ns+1:2*Ns+Ntr); vtrvec = x(2*Ns+Ntr+1:end);
alpha = c1/L * xsvec'*f0vec;
xtvec = [xsvec;xtrvec]; vtvec = [vsvec; vtrvec];

[xj,xi] = meshgrid(xsvec,xtvec);
Kmat = .5*sign(xi-xj);
avec = c1*Kmat*f0vec + alpha + c2 * xtvec;

atotvec = [vtvec;avec];