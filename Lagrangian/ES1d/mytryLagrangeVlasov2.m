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
data = load(input_deck); %this command is stupid but I don't know how to make command like ceil work otherwise


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
%% declare initial distribution
% 
% %finite particles
% N = 2;
% rhobar = -N*q/L;
% x10 = 2.1; x20 = 3; 
% x10p = 0; x20p = 1.2;
% xvec0 = [x10;x20]; 
% vvec0 = [x10p;x20p];
% f0vec = 1/delx/delv*[1;1];
% % end finite particles
% 
% %continuous distribution
% N = Nx;
% 
% lambda = L;
% k = 2*2*pi/lambda;
% deln = .001;
% rho0 = 1;
% vp = 1/k;
% % 
% beam = beams(1);
% rho0 = beam.n0; deln = beam.amplitude; k = 1; vp = 1/k;
% 
% xvec0 = xmesh+.12*delx; % so the 1st order density calculation has reduced spikes
% xvec0 = xvec0'; % get column vector
% vvec0 = 1*vp*deln/delx/(rho0/q)*sin(k*xvec0); 

% vvec0 = 1*vp*deln/rho0*sin(k*xvec0); % tweaking again to get traveling wave
% vvec0 = [.25*ones([N,1])]; 
% f0vec = rho0/q/delv+1/delx/delv * deln*sin(k * xvec0);
% f0vec = rho0/delv+1/delx/delv * deln*sin(k * xvec0); % tweaking again
% % f0vec = rho0/q/delx/N/delv*ones([N,1]);
% % end continuous distribution
% 
% 
% % thermal
% vth = .05;
% lambda = L;
% k = 2*pi/lambda;
% rho0 = 1;
% deln = .001;
% % 
% xlist = xmesh;%+.12*delx; % so the 1st order density calculation has reduced spikes
% vlist= -.1:delv:.1;
% Nv = length(vlist);
% [xvec0,vvec0] = meshgrid(xlist,vlist);
% N = numel(xvec0);
% % trying a shear or offset initialization
% xvec0 = xvec0 + 0*vvec0;
% xvec0 = mod(xvec0-xmin,L)+xmin;
% %
% f0vec = 1/vth/sqrt(pi)*exp(-vvec0.^2/vth^2).*(rho0/q/L);%+1/delx * deln*sin(k * xvec0));
% % f0vec = rho0/q/delv/delx/N*ones([N,1]);
% 
% xvec0 = reshape(xvec0,[N,1]);
% vvec0 = reshape(vvec0,[N,1]);
% f0vec = reshape(f0vec,[N,1]);
% % end thermal initialization

rhobar = -q/L*delx*delv*sum(f0vec);
key_params = [key_params,'vth','k','rho0','deln','Nv','N'];




prerun(xvec0,vvec0,f0vec,pointsize,xmin,L,figure_font,figure_name,save_movie)

%% Solve
% need solver to respect periodic boundary conditions
% use RK 2nd order for now
%un+1 = un + h(b f(un) + c f(un + a h f(un) )
% c = 1/2a, b = 1-c
a = 1; c = 1/2/a; b = 1-c;
%constants for solver
c1 = q^2*delx*delv / m;
c2 = rhobar*q/m;

% Ntfloat = tf/delt; Nt = ceil(Ntfloat);
Nt = ceil(data.tf/data.delt); 
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

for ii = 1:Nt
    x = soln(:,ii);
    fn = odef(x,f0vec,c1,c2,L);
    x2 = x + a*delt*fn;
    x2(1:N) = mod(x2(1:N)-xmin,L)+xmin;
    f2 = odef(x2,f0vec,c1,c2,L);
    v = (b*fn + c*f2);
    
    x = x + delt*v;
    x(1:N) = mod(x(1:N)-xmin,L)+xmin;
    soln(:,ii+1) = x;
    
    %diagnostics
    if plot_in_run
        edensity = xweight(x(1:N), f0vec, xmesh, delx,xmin,delv,q);
        density = edensity + rhobar;
        densitytot(:,ii+1) = density;
        E = K*density;
        phi = M*density;
        Etot(:,ii+1) = E;
        phitot(:,ii+1) = phi;
        
%         edistribution = xvweight(x(1:N), x(N+1:end), f0vec, xmesh, delx,xmin,delv, .3);
        figure(1)
        subplot(2,1,1)
%         plot(x(1:floor(N/10):N),x(N+1:floor(N/10):end),'o' )
%         
%         title('Single particles in phase space')
%         xlim([0,L])
%         xlabel('x\omega_p/v_0')
%         ylabel('v/v_0')
%         ylim([-.1,.1])
%         set(gca,'fontsize', figure_font)
        
        
    scatter(x(1:N),x(N+1:end),pointsize, f0vec);
    xlim([xmin,xmin+L])
    ylim([-.5,.5])
    title(sprintf('Phase space at time = %.02f',ii*delt));
    colorbar()
    xlabel('x')
    ylabel('v')
    set(gca,'fontsize',figure_font)
    
    if save_movie
        frame = getframe(gcf);
        writeVideo(Lagrangev,frame)
    end
%         
%         subplot(3,1,2)
%             scatter(x(1:N),x(N+1:end),pointsize, f0vec);
%     xlim([2,2.5])
%     title('Phase space zoom in')
%     colorbar()
%     xlabel('x')
%     ylabel('v')
%     set(gca,'fontsize',figure_font)
        
%         imagesc([xmin,xmin+L],[-.3,.3],edistribution)
%         title('phase space')
%         xlabel('x\omega_p/v_0')
%         ylabel('v/v_0')
%         colorbar()
%         set(gca,'fontsize', figure_font, 'YDir','normal')
% 
        subplot(2,1,2)
        plot(xmesh,density,xmesh, E)%,xmesh,phi)
        title(sprintf('Time = %f',ii*delt));
       %ylim([-deln/delx,deln/delx])
%        ylim([-.1,.1])
        legend('density','field')%,'potential')
        xlabel('x \omega_p/v_0')
        
        set(gca,'fontsize', figure_font)
        
        pause(.01)
        
    end
end
tlist = delt*(0:Nt);
if save_movie
    % close panel movie writer
    close(Lagrangev);
    
end

if plot_dephi
    % plot density
    figure(2)
    subplot(2,1,1)
    imagesc([0,data.tf], [data.xmin,data.xmin+data.L],densitytot)
    colorbar()
    title('Charge Density \cdot e\omega_p/v_0')
    xlabel('t\omega_p')
    ylabel('x v_0/\omega_p')
    set(gca,'fontsize', figure_font,'YDir','normal')

    subplot(2,1,2)
    imagesc([0,data.tf], [data.xmin,data.xmin+data.L],Etot)
    colorbar()
    title('E Field \cdot ')
    xlabel('t\omega_p')
    ylabel('x v_0/\omega_p')
    set(gca,'fontsize', figure_font,'YDir','normal')
% 
%     subplot(3,1,3)
%     imagesc([0,tf], [xmin,xmin+L],phitot)
%     colorbar()
%     title('Potential \cdot ')
%     ylabel('x v_0/\omega_p')
% 
%     xlabel('t\omega_p')
    set(gca,'fontsize', figure_font,'YDir','normal')
    
    if save_movie   
        print([figure_name 'DE'],'-dpng')
    end
end
    
    
%plot particles
if plot_part
figure(3)
plot(tlist,soln(1:N,:),'.')%,delt*1:Nt,soln(2,:),'o')
ylim([0,L])
xlim([0,data.tf])
title('Phase space particle positions')
xlabel('t\omega_p')
ylabel('x v_0/\omega_p')
set(gca,'fontsize', figure_font)
% % two particle test

if plot_two && ~isempty(spots)

omega = q*sqrt(1/m/L);
c3 = (spots(2).v + spots(1).v) / 2;
c1 = (spots(1).v - c3)/omega;
c2 = 1/2*(spots(1).x - spots(2).x + L/2);
c4 = spots(1).x - c2;

p1 = c1*sin(omega*tlist) +c2*cos(omega*tlist) + c3*tlist + c4;
p2 = L/2 + c3*tlist + c4 - c1*sin(omega*tlist) - c2 * cos(omega*tlist);
p1 = mod(p1-xmin,L)+xmin; p2 = mod(p2-xmin,L)+xmin;
    
    
    hold on
plot(tlist,p1,'o',tlist,p2,'o')
hold off

legend('p1 computed', 'p2 computed', 'p1 analytic', 'p2 analytic')
end
set(gca,'fontsize', figure_font)
% % end two particle test
end


if plot_phase
    figure(4)
    plot(soln(1,:),soln(N+1,:),'.')
    for n = 1:N
        hold on
        plot(soln(n,:),soln(N+n,:),'.')
        hold off
    end
    ylim([-.005,.005])
    xlim([0,L])
    title('phase space orbits')
    ylabel('v/v_0')
    xlabel('x v_0/\omega_p')
    set(gca,'fontsize', figure_font)
end

if save_movie
   % anytime we save a figure or movie, we also want a file with key parameters
   filename = [figure_name 'key_parameters.txt'];
   fileID = fopen(filename,'w');
   for ii = 1:length(key_params)
       param_name = key_params{ii};
       fprintf(fileID,[param_name ' = %.02f\n'], eval(param_name) );
   end
   fclose(fileID);
end