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

plot_init = 1; plot_running = 1; plot_two = 0; plot_energy = 1; plot_rep = 1;
figure_font = 22;

% set mesh
L = 2.*pi; xmin = 0; xmax = xmin + L;
N_mesh = 32; %p  %number of subintervals of mesh
delx_mesh = L/(N_mesh);
% set mesh points, these are the centers of the subintervals
x_mesh = delx_mesh*(1:N_mesh) - .5*delx_mesh;

%for visualization, set vmax
vmax = .2;

delt = .1; tf = 5*pi;  Nt = floor(tf/delt);
tlist = (0:delt:tf-delt)';
if(length(tlist)<Nt)
    tlist = [tlist; tlist(end)+delt];
end
xlist = (xmin+delx_mesh*.5:delx_mesh:xmin+L)';

few_parts = 1; % otherwise do continuous distribution
if few_parts
% % %for a few particles
    Np = 2;
    charges = [1,-1]'; masses = [1,1]'; positions = [2.1,4]'; velocities = [0,0]';
%     charges = [1,1]'; masses = [1,1]'; positions = [2.1,3]'; velocities = [0,1.2]';
    rho0 = sum(charges)/L;

    %two particle test
    if Np > 1
    r0 = positions(1) - positions(2);
    omega = sqrt(charges(1)*charges(2)/L*(masses(1)+masses(2))/masses(1)/masses(2) );
%     omega = 1;
    c2 = (positions(1)-positions(2)+L/2)/(1+masses(1)/masses(2)) ;
    c4 = positions(1) - c2;
    c3 = (masses(2)*velocities(2) + masses(1)*velocities(1))/(masses(1)+masses(2));
    c1 = (velocities(1) - c3)/omega;
    tf = 1.5*2*pi; %/omega; 
    delt = .05; Nt = floor(tf/delt);
    tlist = (0:delt:tf-delt)';
    if(length(tlist)<Nt)
        tlist = [tlist; tlist(end)+delt];
    end
    end
end

% for a distribution
if ~ few_parts
    ppcell = 2;
    Np = N_mesh*ppcell;
    delxp = L/Np;
    deln = .01;
    lam = 2*pi;
    k = 2*2*pi/lam;
    vp = 1/k;
    rho0 = 1.;
    
    density0 = rho0 + deln*sin(k*xlist);
    charges = (rho0*delxp*ones([Np,1]) );
    masses = charges;
    positions =  .337*delx_mesh+delxp* (0:Np-1)';
    positions = positions + deln/k/rho0*cos(k*positions);
    positions = positions ;
    positions = mod(positions - xmin,L) + xmin;
    velocities = zeros([Np,1]);
    velocities = vp*deln/rho0*sin(k*positions);
end


%half step backward to stagger velocities
density = weight(positions, velocities, charges, N_mesh, xmin, delx_mesh);
[E,phibar,dbar]= field_solve(density, delx_mesh);
accel = field2accel(E, positions, charges, masses, xmin, delx_mesh,N_mesh);
velocities = velocities - .5*delt*accel;

%checking density
if plot_init
    figure(1)
    plot(xlist,density,'o')
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
        figure(2)
        subplot(2,1,1)
        plot(xlist,density-rho0,xlist,E,xlist,phi)
        title(sprintf('time = %f',count*delt))
        xlim([xmin,xmin+L])
        ylim([-1,3])
        xlabel('xv_0/\omega_p')
        legend('\rho v_0/(\omega_p e)','E e/(m_e v_0 \omega_p','\phi e/(m_ev_0^2)')
        set(gca,'fontsize', figure_font)

        subplot(2,1,2)
%         fxv = xvdistribution(positions, velocities, charges, N_mesh, xmin, delx_mesh,vmax);
%         imagesc([xmin+1*delx_mesh,xmax-0*delx_mesh],[min(velocities),max(velocities)],fxv');
        hold on
        plot(positions(1:40:end),velocities(1:40:end),'ro')
        plot(positions(end),velocities(end),'ro')
        hold off
        title('phase space')
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
    
    potential(count) = -.5*real(delx_mesh/N_mesh*dbar'*conj(phibar));
    kinetic(count) = .5*sum(masses.*velocities.*velocitiestot(:,count));
    
end



if plot_rep
    figure(5)
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
%     x1 = positions(1)-c2 + c2*cos(omega*tlist) ;
%     x2 = L/2.+ positions(1)-c2 + c2*(-masses(1)/masses(2))*cos(omega*tlist);
%     x1 = mod(x1-xmin,L)+xmin; x2 = mod(x2-xmin,L)+xmin;
    
    %did I get these wrong?
    x1 = c1*sin(omega*tlist) + c2*cos(omega*tlist) + c3*tlist + c4;
    x2 = L/2. + c3*tlist + c4 - c1*masses(1)/masses(2)*sin(omega*tlist) -c2*masses(1)/masses(2)*cos(omega*tlist);    
    x1 = mod(x1-xmin,L)+xmin; x2 = mod(x2-xmin,L)+xmin;    

    figure(3)
    plot(tlist,positionstot','.')
    title('Particle positions over time')
    xlabel('t \omega_p')
    ylabel('x v_0/\omega_p')
    hold on
    plot(tlist,x1,'o',tlist,x2,'o')
    hold off
    ylim([xmin, xmin + L])
    legend('P 1 computed','P2 computed','P1 analytic','P2 analytic')
    set(gca,'fontsize', figure_font)
end

if plot_energy
    figure(4)
    plot(tlist,potential,tlist,kinetic,tlist,potential+kinetic)
    title('Energy \epsilon_0/(m_e v_0^2)')
    legend('Potential','kinetic','total')
    xlabel('t /\omega_p')
    set(gca,'fontsize', figure_font)
    
%     figure
%     plot(tlist,maxdens)
%     title('max of density')
%     xlabel ('t /\omega_p')
%     set(gca,'fontsize', figure_font)
end
