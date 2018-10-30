function [E,overlap] = odef_uniformf0(x,ode_params)


%%% right-hand side of Vlasov 1dES Lagrangian particle method ODEs 
%%% assumes particles initialized so that weighting f0 is uniform
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

%%% 

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

%temp = sort(xsvec);
[~, indx,ind_sorted] = unique(xsvec);
%[~,ind_sorted] = sort(xsvec);
overlap = 0;
E_interact = zeros(size(indx));
if length(indx)==Ns % if no points are on top of eachother
    f0_sorted = f0vec(indx);
    E_interact(1) = sum(f0_sorted(2:end));
    
else
   f0_sorted = zeros(size(indx));
   for ii = 1:Ns
       f0_sorted(ind_sorted(ii)) = f0_sorted(ind_sorted(ii)) + f0vec(ii);
   end
    overlap = 1;
    'You have sheet overlap!'
end

E_interact(1) = -sum(f0_sorted(2:end));
for ii = 2:length(indx)
    E_interact(ii) = E_interact(ii-1) + f0_sorted(ii)+f0_sorted(ii-1);
end

E = .5*c1*E_interact(ind_sorted) + alpha + c2*xtvec;
