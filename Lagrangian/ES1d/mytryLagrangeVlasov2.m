%%%mytryLagranageVlasov2
%%%(needs to be renamed)

%%% implement a semi-lagrangian approach to solving 
%%% 1x1p collisionless electrostatic Vlasov equation.

%%% Ryan Sandberg
%%% started March 2018
%%% last updated May 10, 2018

%%% Vlasov is f' + v. gradx f + a. gradv f = 0
%%% f' + v dfdx + q/m (E) dfdv = 0

%%% Krasny's steps
%%% 1) initialize grid <-> get initial positions of particles
%%%     particles follow phase space trajectories xi(t), vi(t)
%%% 2) solve ode's: xi' = vi, vi' = F(xi) =
%%%     q^2 sum (k(xi,xj) + xj) f0(aj,bj) delta a delta b + q rhobar xi
%%% rhobar = background ion charge density (i.e. this is an electron
%%%     vlasov solver


% desired outline:
% generate input deck
% generate f0vec
% load input deck
% generate derived data

function mytryLagrangeVlasov2(input_deck)
load(input_deck)
input_data = load(input_deck); %this command is stupid but I don't know how to make command like ceil work otherwise


figure_name = ['../output_files/' run_day '/' run_name '/' run_name '_'];
movie_name = [figure_name 'phase_space.avi'];
key_params = {};


%% set up mesh
 delx = L/Nx; key_params = [key_params,'xmin','L','Nx','delx'];
xmesh = xmin + .5*delx : delx : xmin + L;
key_params = [key_params,'delv'];



%% Physical Constants

key_params = [key_params,'m','q'];

x0 = xmin+.5*delx : delx : xmin + L;
v0 = -vmax + .5*delv : delv : vmax; Nv = length(v0);

N = Nx*Nv;
[xvec0,vvec0] = meshgrid(x0, v0);
xvec0 = reshape(xvec0,[N,1]);
vvec0 = reshape(vvec0,[N,1]);
f0vec = reshape(f0vec,[N,1]);


rhobar = -q/L*delx*delv*sum(f0vec);
key_params = [key_params,'vth','k','rho0','deln','Nv','N'];


save('output_data')
prerun()


Nt = ceil(input_data.tf/input_data.delt); 
key_params = [key_params,'tf','delt', 'Nt'];
soln = zeros(2*N,Nt+1);
soln(1:N,1) = xvec0; soln(N+1:end,1) = vvec0; 
%% f ( xvec, vvec)

%initialize diagnostics
densitytot = zeros([Nx,Nt+1]);
Etot = zeros([Nx,Nt+1]);
phitot  = zeros([Nx,Nt+1]);
edensity = xweight(xvec0, f0vec, xmesh, delx,xmin,delv,q);
density = edensity + rhobar;
densitytot(:,1) = density;

if save_movie
    Lagrangev = VideoWriter(movie_name);
    open(Lagrangev)
end


xvec = xmesh';
X = xvec * ones([1,Nx]);
Y = X';
D = X-Y;
K = delx*(.5*sign(D) + Y/L);
M = -delx*(.5*abs(D) + 1/L*(xvec*xvec'));

E = K*density;
phi = M*density;
Etot(:,1) = E;
phitot(:,1) = phi;

ode_params.function = 'odef_tracer';
ode_params.f0vec = f0vec;
ode_params.c1 = q^2*delx*delv / m;
ode_params.c2 = rhobar*q/m;
ode_params.L = L;
ode_params.Ntr = 0;

plot_data = struct('pointsize',pointsize,'f0vec',f0vec,'xmin',xmin,...
    'L',L,'delt',delt, 'figure_font',figure_font,'xmesh',xmesh,...
    'N',N,'delv',delv,'delx',delx,'xvec0',xvec0,'vvec0',vvec0,'Nx',Nx,'Nv',Nv);



for ii = 1:Nt
    x = soln(:,ii);
    
    x = ode_int(x,ode_params,method_params);

    soln(:,ii+1) = x;
    
    % calculate diagnostic information to save for later

    edensity = xweight(x(1:N), f0vec, xmesh, delx,xmin,delv,q);
    density = edensity + rhobar;
    densitytot(:,ii+1) = density;
    E = K*density;
    phi = M*density;
    Etot(:,ii+1) = E;
    phitot(:,ii+1) = phi;

    %plot diagnostics
    plot_data.x = x; plot_data.E = E; plot_data.density = density;
    plot_data.time = ii*delt;
    inrun(plot_in_run,save_movie,subplot_array,plot_data)


    pause(.01)
        
end

if save_movie
    % close panel movie writer
    close(Lagrangev);
    
end

save('output_data')

postrun
