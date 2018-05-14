%%% implement a semi-lagrangian approach to solving 
%%% 1x1p collisionless electrostatic Vlasov equation.
%%% Then put in a tree code

%%% Vlasov is f' + v. gradx f + a. gradv f = 0
%%% f' + v dfdx + q/m (E) dfdv = 0

%%% Krasny's steps
%%% 1) initialize grid <-> get initial positions of particles
%%%     particles follow phase space trajectories xi(t), vi(t)
%%% 2) solve ode's: xi' = vi, vi' = F(xi) =
%%%     q^2 sum (k(xi,xj) + xj) f0(aj,bj) delta a delta b + q rhobar xi
%%% rhobar = background ion charge density (i.e. this is an electron
%%%     vlasov solver

function [part1,part1an] = LagrangeVlasovFun(L,Nx,delv,tf,delt, x10,x10p,x20,x20p)

figure_font = 22;

%% set up mesh
xmin = 0; %L = 2*pi; Nx = 50; 
delx = L/Nx;
xmesh = xmin + .5*delx : delx : xmin + L;
%delv = .001; 
vmax = .5;

plot_in_run = 0;

%% Physical Constants
m = 1;
q = 1;



%% declare initial distribution

%finite particles
N = 2;
rhobar = -N*q/L;
% x10 = .3; x20 = .7; 
% x10p = -.6; x20p = 0;
xvec0 = [x10;x20]; 
vvec0 = [x10p;x20p];
f0vec = 1/delx/delv*[1;1];

%continuous distribution
% N = Nx;
% 
% lambda = L;
% k = 2*2*pi/lambda;
% deln = .001;
% rho0 = 1;
% vp = 1/k;
% 
% xvec0 = xmesh+.12*delx; % so the 1st order density calculation has reduced spikes
% xvec0 = xvec0'; % get column vector
% vvec0 = 1*vp*deln/delx/(rho0/q)*sin(k*xvec0); %zeros([N,1]); 
% f0vec = rho0/q/delv+1/delx/delv * deln*sin(k * xvec0);
% % f0vec = 1/delv + 1/delx/delv*ones([N,1]);

% thermal
% vth = .01;
% lambda = L;
% k = 2*pi/lambda;
% rho0 = 1;
% deln = .001;
% 
% xlist = xmesh+.12*delx; % so the 1st order density calculation has reduced spikes
% vlist= -.03:delv:.03;
% Nv = length(vlist);
% [xvec0,vvec0] = meshgrid(xlist,vlist);
% f0vec = exp(-vvec0.^2/vth^2).*(rho0/q/delv+1/delx/delv * deln*sin(k * xvec0));
% N = numel(f0vec);
% xvec0 = reshape(xvec0,[N,1]);
% vvec0 = reshape(vvec0,[N,1]);
% f0vec = reshape(f0vec,[N,1]);

rhobar = -q/L*delx*delv*sum(f0vec);


%% Solve
% need solver to respect periodic boundary conditions
% use RK 2nd order for now
%un+1 = un + h(b f(un) + c f(un + a h f(un) )
% c = 1/2a, b = 1-c
a = 1; c = 1/2/a; b = 1-c;
%constants for solver
c1 = q^2*delx*delv / m;
c2 = rhobar*q/m;

% tf = 5*2*pi;
%delt = .1;
Nt = ceil(tf/delt);
soln = zeros(2*N,Nt+1);
soln(1:N,1) = xvec0; soln(N+1:end,1) = vvec0; 
%% f ( xvec, vvec)

%initialize diagnostics
densitytot = zeros([Nx,Nt+1]);
Etot = zeros([Nx,Nt+1]);
phitot  = zeros([Nx,Nt+1]);
edensity = xweight(xvec0, f0vec, xmesh, delx,xmin,delv,q);
densitytot(:,1) = edensity+rhobar;

xvec = xmesh';
X = xvec * ones([1,Nx]);
Y = X';
D = X-Y;
K = .5*sign(D) + Y;
M = -.5*abs(D) - xvec*xvec';

for ii = 1:Nt
    x = soln(:,ii);
    fn = odef(x,f0vec,c1,c2,N,L);
    x2 = x + a*delt*fn;
    x2(1:N) = mod(x2(1:N)-xmin,L)+xmin;
    f2 = odef(x2,f0vec,c1,c2,N,L);
    v = (b*fn + c*f2);
    
    x = x + delt*v;
    x(1:N) = mod(x(1:N)-xmin,L)+xmin;
    soln(:,ii+1) = x;
    
    %diagnostics
    if plot_in_run
        edensity = xweight(x(1:N), f0vec, xmesh, delx,xmin,delv,q);
        density = edensity + rhobar;
        densitytot(:,ii) = density;
        E = K*density*delx;
        phi = M*density*delx;
        Etot(:,ii) = E;
        phitot(:,ii) = phi;
        
        edistribution = xvweight(x(1:N), x(N+1:end), f0vec, xmesh, delx,xmin,delv, vmax);
        figure(1)
        title(sprintf('Time = %f',ii*delt));
        subplot(3,1,1)
        title('Single particles in phase space')
        plot(x(1:10:N),x(N+1:10:end),'o' )
        xlim([0,.5])
        ylim([-vmax,vmax])
        set(gca,'fontsize', figure_font)
        
        
        subplot(3,1,2)
        title('phase space')
        imagesc([xmin,xmin+L],[-vmax,vmax],edistribution)
        colorbar()
        set(gca,'fontsize', figure_font, 'YDir','normal')
        pause(.01)

        subplot(3,1,3)
        plot(xmesh,density,xmesh, E)
%        ylim([-deln/delx,deln/delx])

        legend('density','field')
        set(gca,'fontsize', figure_font)
        
    end
end
tlist = delt*(0:Nt);


% plot density
% figure(2)
% subplot(3,1,1)
% imagesc([0,tf], [xmin,xmin+L],densitytot)
% colorbar()
% set(gca,'fontsize', figure_font,'YDir','normal')
% title('Charge Density \cdot e\omega_p/v_0')
% ylabel('x v_0/\omega_p')
% 
% subplot(3,1,2)
% imagesc([0,tf], [xmin,xmin+L],Etot)
% colorbar()
% set(gca,'fontsize', figure_font,'YDir','normal')
% title('E Field \cdot ')
% ylabel('x v_0/\omega_p')
% 
% subplot(3,1,3)
% imagesc([0,tf], [xmin,xmin+L],phitot)
% colorbar()
% set(gca,'fontsize', figure_font,'YDir','normal')
% title('Potential \cdot ')
% ylabel('x v_0/\omega_p')
% 
% xlabel('t\omega_p')
% 
% %plot particles
% figure(3)
% plot(tlist,soln(1:N,:),'.')%,delt*1:Nt,soln(2,:),'o')
% ylim([0,L])
% xlim([0,tf])
% title('Phase space particle positions')
% xlabel('t\omega_p')
% ylabel('x v_0/\omega_p')
% set(gca,'fontsize', figure_font)
% % % two particle test
% 
omega = q*sqrt(2/m/L);
c3 = (x20p + x10p) / 2;
c1 = (x10p - c3)/omega;
c2 = 1/2*(x10 - x20 + L/2);
c4 = x10 - c2;

p1 = c1*sin(omega*tlist) +c2*cos(omega*tlist) + c3*tlist + c4;
p2 = L/2 + c3*tlist + c4 - c1*sin(omega*tlist) - c2 * cos(omega*tlist);
p1 = mod(p1-xmin,L)+xmin; p2 = mod(p2-xmin,L)+xmin;

part1 = x(1);
part1an = p1(end);

% hold on
% plot(tlist,p1,'o',tlist,p2,'o')
% hold off
% legend('p1 computed', 'p2 computed', 'p1 analytic', 'p2 analytic')
% set(gca,'fontsize', figure_font)
% % % end two particle test