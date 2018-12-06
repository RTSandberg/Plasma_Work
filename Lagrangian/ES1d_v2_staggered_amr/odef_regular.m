function [E,phi] = odef_regular(x,ode_params)


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
rhobar = ode_params.rhobar; 
L = ode_params.L;

Ns = length(f0vec);
Ntr = ode_params.Ntr;

deltasq = (ode_params.deltasq); 

xsvec = x(1:Ns); vsvec = x(Ns+1:2*Ns);
xtrvec = x(2*Ns+1:2*Ns+Ntr); vtrvec = x(2*Ns+Ntr+1:end);
alpha = c1/L * xsvec'*f0vec;
xtvec = [xsvec;xtrvec]; vtvec = [vsvec; vtrvec];


E_interact = zeros(size(xtvec));
for ii = 1:Ns+Ntr

    E_interact(ii) = f0vec'* ...
        ((xtvec(ii)-xsvec)./sqrt((xtvec(ii)-xsvec).^2 + deltasq));
        
end
    

E = .5*c1*E_interact + alpha + rhobar*xtvec;

if ode_params.get_potential
        phi = zeros(size(xtvec));
        
        for ii = 1:Ns+Ntr
            tempvec = (xtvec(ii)-xsvec').^2;
            phi(ii) = -c1/2*(sqrt(tempvec+deltasq)-tempvec/L)*f0vec;
        end
        
%         phi = phi - 3;
%         phi = -c1/2*(phi_interact - phi_interact.^2/L)*f0vec;
        %xtvec.*xtvec*rhobar/2 - alpha*xtvec;
        
else
    phi = zeros(size(xtvec));
end