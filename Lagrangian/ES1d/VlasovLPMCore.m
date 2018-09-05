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

function VlasovLPMCore(input_deck)
tic
load(input_deck)
input_data = load(input_deck); %this command is stupid but I don't know how to make command like ceil work otherwise



topic = '2_stream';
run_day = 'July_2_2018';
run_name = 'equilibrium';
figure_name = ['../../../output_s/' topic '/' run_day '/' run_name '_'];
movie_name = [figure_name 'phase_space.avi'];
Lagrangev = 1;
key_params = {};


%% set up mesh
 delx = L/Nx; key_params = [key_params,'xmin','L','Nx','delx'];
xmesh = xmin + .5*delx : delx : xmin + L;
key_params = [key_params,'delv'];



%% Physical Constants

key_params = [key_params,'m','q'];

x0 = xmin+.5*delx : delx : xmin + L;
v0 = vmin + .5*delv : delv : vmax; Nv = length(v0);

N = Nx*Nv;
[xvec0,vvec0] = meshgrid(x0, v0);
if shear.shearQ % shear initial distribution to offset spikes in E field
    xvec0 = xvec0 + shear.slope*vvec0;
    xvec0 = mod(xvec0 - xmin,L)+xmin;
end

xvec0 = reshape(xvec0,[N,1]);
vvec0 = reshape(vvec0,[N,1]);
f0vec = reshape(f0vec,[N,1]);
N = length(f0vec);

rhobar = -q/L*delx*delv*sum(f0vec);
% key_params = [key_params,'vth','k','rho0','deln','Nv','N'];



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




ode_params.function = 'odef_tracer';
ode_params.f0vec = f0vec;
ode_params.c1 = q^2*delx*delv / m;
ode_params.c2 = rhobar*q/m;
ode_params.L = L;
ode_params.Ntr = 0;
ode_params.rhobar = rhobar;

potential_params = ode_params;
potential_params.function = 'potential_tracer';
potential_params.Ntr = Nx;
potential_params.c1 = q*delx*delv;
potential_params.c2 = rhobar;

xvec = xmesh';
X = xvec * ones([1,Nx]);
Y = X';
D = X-Y;
K = delx*(.5*sign(D) + Y/L);
% M = -delx*(.5*abs(D) + 1/L*(xvec*xvec'));

E = K*density;
% phi = M*density;
% phi = potential_tracer([xvec0;vvec0;xvec;zeros([Nx,1])],potential_params);
Etot(:,1) = E;
% phitot(:,1) = phi(N+1:end);

plot_data = struct('pointsize',pointsize,'f0vec',f0vec,'xmin',xmin,...
    'L',L,'delt',delt, 'figure_font',figure_font,'xmesh',xmesh,...
    'N',N,'delv',delv,'delx',delx,'xvec0',xvec0,'vvec0',vvec0,'Nx',Nx,...
    'Nv',Nv,'ode_params',ode_params,'potential_params',potential_params);



save('../../../big_simulation_data/output_data')
prerun(figure_name, save_movie, prerun_subplot_array,plot_data)

%%% deprecated for now - thoughts on keeping movies to about 300 frames
% diagnostic_increment = 1;
% if Nt > 500
%     diagnostic_increment = floor(Nt/100);
% end
diagnostic_increment = min(Nt, diagnostic_increment);
secondinitialtime = toc

tic

for ii = 1:Nt
    x = soln(:,ii);
    
    [x,v] = ode_int(x,ode_params,method_params);

    soln(:,ii+1) = x;
    
    % calculate diagnostic information to save for later

    edensity = xweight(x(1:N), f0vec, xmesh, delx,xmin,delv,q);
    density = edensity + rhobar;
    densitytot(:,ii+1) = density;
%     E = K*density;
    E = interp1(x(1:N),m/q*v(N+1:end),xmesh);
%     phi = M*density;
%     phi = potential_tracer([x; xvec;zeros([Nx,1])],potential_params);
    Etot(:,ii+1) = E;
%     phitot(:,ii+1) = phi(N+1:end);

    %plot diagnostics
    if ( mod(ii, diagnostic_increment) == 0) && plot_in_run 
        plot_data.x = x; plot_data.E = E; plot_data.density = density;
        plot_data.time = ii*delt;
        inrun(Lagrangev,plot_in_run,save_movie,inrun_subplot_array,plot_data)
    end

    pause(.01)
        
end
actualruntime = toc
tic
if save_movie
    % close panel movie writer
    close(Lagrangev);
    print([figure_name 'phase_final'],'-dpng')
%     savefig([figure_name 'phase_final'])
end

save('../../../big_simulation_data/output_data')
postrunmytryL = toc