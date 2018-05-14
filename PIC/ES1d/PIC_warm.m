%%% PIC for warm plasma
%%% 1d PIC code
% Ryan
% based on Alec's notes
%%%
% Initialize particles
% mesh particles
% solve - DFT
% push field to particles
% advect (1 time step each particles and velocities)


% begin cycle:
% 1) mesh particles
% 2) solve for potential, then field
% 3) push field to particles
% 4) move particles
% 5) repeat

%%% constants

%%% initialize
% N_part
% xlim
% delx_part (meant to vary)
% initial density
% weights
% particle list
% velocity list, x and y directions
% 
% N_mesh
% delx_mesh
% E1, E2, B fields, density
clear('all')
close('all')

plot_init = 1; plot_running = 1; plot_two = 0; plot_energy = 0; plot_rep = 1;
figure_font = 22;

% set mesh
L = 2*pi; xmin = 0; xmax = xmin + L;
N_mesh = 64; %p  %number of subintervals of mesh
delx_mesh = L/(N_mesh);
% set mesh points, these are the centers of the subintervals
x_mesh = delx_mesh*(1:N_mesh) - .5*delx_mesh;



xlist = (xmin+delx_mesh*.5:delx_mesh:xmin+L)';

% for a distribution
ppcell = 64;
Np = N_mesh*ppcell;
delxp = L/Np;
deln = .01;
lam = 2*pi;
k = 2*pi/lam/1;
vp = 1/k;
rho0 = 1.;
vthermal = .01;

%
omega = sqrt(1+3*k^2*vthermal^2);
tf = 2*pi; delt = .1; Nt = ceil(tf/delt);
tlist = (0:delt:tf)';

density0 = rho0 + deln*sin(k*xlist);
charges = (rho0*delxp*ones([Np,1]) );
masses = charges;
positions =  pi*.1*delx_mesh+delxp* (0:Np-1)';
positions = positions + deln/k/rho0*cos(k*positions);
positions = mod(positions - xmin,L) + xmin;
%velocities = vp*deln/rho0*sin(k*positions);
velocities = zeros([Np,1]);
velocities = velocities + gaussian(Np,vthermal);


%half step backward to stagger velocities
density = weight(positions, velocities, charges, N_mesh, xmin, delx_mesh);
[E,phibar,dbar]= field_solve(density, delx_mesh);
accel = field2accel(E, positions, charges, masses, xmin, delx_mesh,N_mesh);
velocities = velocities - .5*delt*accel;

%checking density
if plot_init
    figure(2)
        fxv = xvdistribution(positions, velocities, charges, N_mesh, xmin, delx_mesh,.04);
        imagesc([xmin+1*delx_mesh,xmax-0*delx_mesh],[min(velocities),max(velocities)],fxv');
        title('Initial Distribution')
        xlabel('x')
        ylabel('v')
        colorbar()
        set(gca,'YDir','normal')
end


%initialize diagnostics
positionstot = zeros([Np,Nt]);
velocitiestot = zeros([Np,Nt]);
densitytot = zeros([N_mesh,Nt]);
Etot = zeros([N_mesh,Nt]);
phitot = zeros([N_mesh,Nt]);
kinetic = zeros([Nt,1]);
potential = zeros([Nt,1]);
chargetot = zeros([Nt,1]);
maxdens = zeros([Nt,1]);


% run
for count = 1:Nt

    positionstot(:,count) = positions;
    velocitiestot(:,count) = velocities;
    
    density = weight(positions, velocities, charges, N_mesh, xmin, delx_mesh);
    [E,phibar,dbar]= field_solve(density, delx_mesh);
    phi = real(ifft(phibar));
    accel = field2accel(E, positions, charges, masses, xmin, delx_mesh,N_mesh);
%     if mod(count,125) == 0
%         figure(2)
%         hold on
%         plot(xlist,density)
%         hold off
%     end

    
    %plot during run
    
    if( mod(count,10) == 0 && plot_running)
        figure(1)
        subplot(2,1,1)
        plot(xlist,density-rho0,xlist,E,xlist,phi)
        title(sprintf('time = %f',count*delt))
        xlim([xmin,xmin+L])
        ylim([-2,2])
        xlabel('xv_0/\omega_p')
        legend('\rho v_0/(\omega_p e)','E e/(m_e v_0 \omega_p','\phi e/(m_ev_0^2)')
        set(gca,'fontsize', figure_font)

        subplot(2,1,2)
        fxv = xvdistribution(positions, velocities, charges, N_mesh, xmin, delx_mesh,.04);
        imagesc([xmin+1*delx_mesh,xmax-0*delx_mesh],[min(velocities),max(velocities)],fxv');
        hold on
        plot(positions(1:floor(Np/5):end),velocities(1:floor(Np/5):end),'ro')
        plot(positions(end),velocities(end),'ro')
        hold off
        xlim([xmin,xmin+L])
        xlabel('x\omega_p/v_0')
        ylabel('v/ v_0')
        set(gca,'fontsize', figure_font,'YDir','normal')
        pause(.01)
    end
    
    densitytot(:,count) = density;
    Etot(:,count) = E;
    phitot(:,count) = phi;
    chargetot(count) = delx_mesh*sum(density-rho0);
    maxdens(count) = max(abs(density));
    
    velocities = velocities + delt*accel;
    positions = positions + delt*velocities;
    positions = mod(positions - xmin,L) + xmin;
    
    potential(count) = abs(.5/L*dbar'*conj(phibar));
    kinetic(count) = .5*sum(masses.*velocities.*velocitiestot(:,count));
    
end



if plot_rep
    figure
    subplot(3,1,1)
    imagesc(tlist, xlist, densitytot)
    title('Charge Density \rho v_0/(\omega_p e)')
    xlabel('t/ \omega_p')
    ylabel('x v_0/\omega_p')
    colorbar
    set(gca,'fontsize', figure_font,'YDir','normal')
    
    subplot(3,1,2)
    imagesc(tlist, xlist, Etot)
    title('Electric field E e/(m_ev_0\omega_p)')
    xlabel('t /\omega_p')
    ylabel('x v_0/\omega_p')
    colorbar
    set(gca,'fontsize', figure_font,'YDir','normal')

    subplot(3,1,3)
    imagesc(tlist, xlist, phitot)
    title('Electric potential \phi e/(m_ev_0^2)')
    xlabel('t /\omega_p')
    ylabel('x v_0/\omega_p')
    colorbar
    set(gca,'fontsize', figure_font,'YDir','normal')
end

if plot_two
% x1 = positions(1)-c3 + c3*cos(omega*tlist);
% x2 = L/2.+ positions(1)-c3 + c3*(-masses(1)/masses(2))*cos(omega*tlist);
% x1 = mod(x1-xmin,L)+xmin; x2 = mod(x2-xmin,L)+xmin;

%figure
% plot(tlist,positionstot')
% title('Particle positions over time')
% xlabel('t \omega_p')
% ylabel('x v_0/\omega_p')
% hold on
% plot(tlist,x1,tlist,x2)
% hold off
% legend('P 1 computed','P2 computed','P1 analytic','P2 analytic')
end

if plot_energy
    figure
%     plot(tlist,kinetic)
    plot(tlist,potential,tlist,kinetic,tlist,potential+kinetic)
    title('Energy \epsilon_0/(m_e v_0^2)')
    legend('potential', 'kinetic','total')
    xlabel('t /\omega_p')
    set(gca,'fontsize', figure_font)
    
%     figure
%     plot(tlist,maxdens)
%     title('max of density')
%     xlabel ('t /\omega_p')
%     set(gca,'fontsize', figure_font)
end
