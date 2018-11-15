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
clear
close('all')

tic
plot_init = 1; plot_running = 0; plot_two = 0; plot_energy = 0; plot_rep = 1;
plot_first = 0; plot_periodic = 0; compute_error = 0; oscillator_error = 0;
plot_parts = 0;
figure_font = 22;


topic = 'developing_LPM/1d_test_cases/cold_langmuir_wavebreaking/pic_early_break';
run_day = 'Nov_2_2018';
run_name = 'pic_standing_tf_2tp_Nx_118_ppc_8_Nt_1600_n1_p85';
figure_name = ['../../../PlasmaResearch/output_s/' topic '/'  run_day '/' run_name];
movie_name = [figure_name '_phase_space.avi'];

Nt=1600;

save_fig = 0;
save_final_data = 0;
save_movie = 0;

if save_movie
picmovie = VideoWriter(movie_name);
    open(picmovie)
end

q = -1;
m = 1;
% set mesh
L = 2*pi; xmin = 0; xmax = xmin + L;
N_mesh =284; %p  %number of subintervals of mesh
delx_mesh = L/(N_mesh);
% set mesh points, these are the centers of the subintervals
x_mesh = delx_mesh*(1:N_mesh) - .5*delx_mesh;


convergence_name = [figure_name sprintf('delt_p1_Nx_%i_',N_mesh)];

%for visualization, set vmax
vmax = 1;

tf = 2*2*pi;  delt = tf/Nt;
tlist = (0:Nt)*delt;
% if(length(tlist)<Nt)
%     tlist = [tlist; tlist(end)+delt];
% end
xlist = (xmin+delx_mesh*.5:delx_mesh:xmin+L)';

few_parts = 0; % otherwise do continuous distribution
if few_parts
% % %for a few particles
    Np = 20;
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
    nstreams = 1;
    ppcell = 20;
    Np = nstreams * N_mesh*ppcell;
    delxp = L/N_mesh/ppcell;
    deln = .8;
    lam = L;
    k = 1*2*pi/lam;
    vp = 1/k;
    rho0 = 1.;
    v0 = 1;
    
    density0 = rho0 + deln*sin(k*xlist);
    
    positions0 =  0*delx_mesh+delxp* (0:Np/nstreams-1)';
    positions1 = positions0;%+deln/rho0/k*(cos(k*positions0));%...
%         +sqrt(3)*cos(k*positions0));
    positions2 = positions0;%+ deln/rho0/k*(-sin(k*positions0)...
%         -sqrt(3)*cos(k*positions0));
    positions = [positions1];%;positions2];
%     positions = positions ;
    positions = mod(positions - xmin,L) + xmin;
    velocities = zeros([Np,1]);
%     velocities = v0 * [ones([Np/nstreams,1]);-ones([Np/nstreams,1])];
    velocities = -deln/rho0*cos(k*positions0);
    velocities1 = v0+deln/rho0*(-cos(k*positions0)...
        -1/sqrt(3)*sin(k*positions0));
    velocities2 = -v0+deln/rho0*(-cos(k*positions0)...
        +1/sqrt(3)*sin(k*positions0));
%     velocities = [velocities1;velocities2];

charges = q*rho0*delxp*ones(size(positions));% + deln*delxp*sin(k*positions0 );
% charges = q*rho0*delxp*(1+deln*sin(k*positions)).*(1-deln*sin(k*positions0));
% charges = [charges; charges];
    
    masses = m/q*charges;

% a second stream
    
end


%start computing diagnostics
density = weight(positions, velocities, charges, N_mesh, xmin, delx_mesh);
[E,phibar,dbar]= field_solve(density, delx_mesh);
phi = real(ifft(phibar));
current = weight(positions, velocities, velocities.*charges, N_mesh, xmin, delx_mesh);

%checking density
if plot_init
    figure(1)
    subplot(2,1,1)
    plot(positions,velocities,'.')
%         plot(positions(1:Np/nstreams),velocities(1:Np/nstreams),'r.')
%         hold on
%         plot(positions(Np/nstreams+1:end),velocities(Np/nstreams+1:end),'b.')
%         hold off
    title('Phase space')
    xlabel('x/(v_0\omega_p)')
    ylabel('v/v_0')
    set(gca,'fontsize',figure_font)
    
    subplot(2,1,2)
    plot(xlist, density-nstreams*q/L)
    title('Initial density')
    xlabel('x/(v_0\omega_p)')
    ylabel('v/v_0')
    set(gca,'fontsize',figure_font)
end


%initialize diagnostics
positionstot = zeros([Np,Nt+1]);
velocitiestot = zeros([Np,Nt+1]);
densitytot = zeros([N_mesh,Nt+1]);
currenttot = zeros([N_mesh, Nt+1]);
Etot = zeros([N_mesh,Nt+1]);
phitot = zeros([N_mesh,Nt+1]);
kinetic = zeros([Nt+1,1]);
potential = zeros([Nt+1,1]);
chargetot = zeros([Nt+1,1]);
maxdens = zeros([Nt+1,1]);

    densitytot(:,1) = density-q/L;
    currenttot(:,1) = current;
    Etot(:,1) = E;
    phitot(:,1) = phi;
    chargetot(1) = delx_mesh*sum(density-rho0);
    maxdens(1) = max(abs(density));
    

    
    positionstot(:,1) = positions;
    velocitiestot(:,1) = velocities;
    
%     potential(1) = -.5*delx_mesh/N_mesh*real(dbar'*phibar);
    potential(1) = .5*delx_mesh*(E'*E);
%     [velocities, velocitiestot(:,count)];
    kinetic(1) = .5*sum(masses.*velocities.*velocitiestot(:,1));
    %%% done initializing diagnostics

%take half step backward to stagger positions, velocities    
accel = field2accel(E, positions, charges, masses, xmin, delx_mesh,N_mesh);
velocities = velocities - .5*delt*accel;
velocitiestot(:,1) = velocities;

    
    thetas = 0:.1:2*pi;
    crossed = 0;
    
    
picpreproc = toc
tic
% run
for count = 2:Nt+1

    
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
    
    if( count >= 1 && mod(count,1) == 0 && plot_running)
        figure(2)
        set(gcf,'Position',[7 91 1059 894])
%         set(gcf,'Position',[192 86 840 617]);
        subplot(3,1,1)
        plot(xlist,density-q*nstreams*rho0,'*',xlist,E,'o');%,xlist,phi)
%         plot(xlist, E)
        title(sprintf('time = %f',count*delt))
        xlim([xmin,xmin+L])
        ylim([-.5,.5])
        xlabel('xv_0/\omega_p')
        legend('\rho v_0/(\omega_p e)','E e/(m_e v_0 \omega_p)');%,'\phi e/(m_ev_0^2)')
%         legend('E e/(m_e v_0\omega_p)')
        set(gca,'fontsize', figure_font)

        subplot(3,1,2)
%         fxv = xvdistribution(positions, velocities, charges, N_mesh, xmin, delx_mesh,vmax);
%         imagesc([xmin+1*delx_mesh,xmax-0*delx_mesh],[min(velocities),max(velocities)],fxv');
        
        scatter(positions,velocities,10,charges)
        % see how many particles in a cell
%         hold on
%         for ii=1:N_mesh
%            plot((ii-.5)*delx_mesh*[1 1],[-1,1],'k') 
%         end
%         hold off
        %
%         hold on
%         plot(positions(1:1000:end),velocities(1:1000:end),'r*')
%         hold off
%         plot(positions(1:Np/nstreams),velocities(1:Np/nstreams),'r.')
%         hold on
%         plot(positions(Np/nstreams+1:end),velocities(Np/nstreams+1:end),'b.')
%         hold off
        %         hold on
%         plot(positions(end),velocities(end),'ro')
%         hold off
        title('phase space')
        xlim([xmin,xmin+L])
        ylim([-1.2,1.2])
        xlabel('x\omega_p/v_0')
        ylabel('v/ v_0')
        set(gca,'fontsize', figure_font,'YDir','normal')
        
        subplot(3,1,3)
        % spatial fourier spectrum of E
        freqs = 2*pi*(-N_mesh/2:N_mesh/2-1)/delx_mesh/N_mesh;
        Ebar = fft(E);
        plot(freqs,fftshift(abs(Ebar)))
        title('Fourier modes of E')
        xlim([-5,5])
        xlabel('k')
        set(gca,'fontsize',figure_font)
        
    % using first particle as proxy for wave
%         plot(positions(1)-positions0(1),velocities(1),'*')
%         hold on
%         plot(deln*cos(thetas),deln*sin(thetas))
%         hold off
%         xlim([-2*deln,2*deln])
%         ylim([-2*deln,2*deln])
        
if save_movie
        frame = getframe(gcf);
        writeVideo(picmovie,frame)
end
%         pause(.01)
    end
    
    densitytot(:,count-1) = density-q/L;
    Etot(:,count-1) = E;
    phitot(:,count-1) = phi;
    chargetot(count-1) = delx_mesh*sum(density-rho0);
    maxdens(count-1) = max(abs(density));
    
    velocities = velocities + delt*accel;
    positions = positions + delt*velocities;
    positions = mod(positions - xmin,L) + xmin;
    
    positionstot(:,count) = positions;
    velocitiestot(:,count) = velocities;
    velocity_count = (velocities + velocitiestot(:,count-1))/2;
    currenttot(:,count-1) = weight(positions, (velocity_count), ...
        velocity_count.*charges, N_mesh, xmin, delx_mesh);
    
%     potential(count) = -.5*delx_mesh/N_mesh*real(dbar'*phibar);
    potential(count-1) = .5*delx_mesh*(E'*E);
%     [velocities, velocitiestot(:,count)];
    kinetic(count-1) = .5*sum(masses.*velocities.*velocitiestot(:,count-1));
    
    % check for crossing\
    if length(find(positions(2:end)-positions(1:end-1)<0))>1
        sprintf('you had sheet crossing at time %.02f',count*delt)
        break
    end
    
end

% if crossed
%     'you had sheet crossing'
% end

if(save_final_data)
    solnfinal = [positions;velocities];
    save([figure_name '_solnfinal'], 'solnfinal');
end
    
    % advance velocity on more step to get the average
    
    density = weight(positions, velocities, charges, N_mesh, xmin, delx_mesh);
    [E,phibar,dbar]= field_solve(density, delx_mesh);
    accel = field2accel(E, positions, charges, masses, xmin, delx_mesh,N_mesh);
    vnew = velocities + delt*accel;
    vcount = (vnew + velocities)/2;
    currenttot(:,count) = weight(positions, vcount, vcount.*charges, N_mesh, xmin, delx_mesh);
    
    
    densitytot(:,count) = density-q/L;
    Etot(:,count) = E;
    phitot(:,count) = phi;
    chargetot(count) = delx_mesh*sum(density-rho0);
    maxdens(count) = max(abs(density));
    
    
runtime = toc

tic

if(save_fig)
    solnfinal = [positions;velocities];
    save([convergence_name 'solnfinal'], 'solnfinal');
end

figure
subplot(2,1,1)
        plot(xlist,density-nstreams*rho0,xlist,E)
%         plot(xlist, E)
        title(sprintf('time = %f',count*delt))
        xlim([xmin,xmin+L])
%         ylim([-.005,.005])
        xlabel('xv_0/\omega_p')
        legend('\delta\rho v_0/(\omega_p e)','E e/(m_e v_0 \omega_p)');%,'\phi e/(m_ev_0^2)')
%         legend('E e/(m_e v_0\omega_p)')
        set(gca,'fontsize', figure_font)

subplot(2,1,2)
plot(positions,velocities,'.')
%     plot(positions(1:Np/nstreams),velocities(1:Np/nstreams),'r.')
%     hold on
%     plot(positions(Np/nstreams+1:end),velocities(Np/nstreams+1:end),'b.')
%     hold off
    %         hold on
%         plot(positions(end),velocities(end),'ro')
%         hold off
    title('final phase space')
    xlim([xmin,xmin+L])
    ylim([-3,3])
    xlabel('x\omega_p/v_0')
    ylabel('v/ v_0')
    set(gca,'fontsize', figure_font,'YDir','normal')

%     figname = ['../../../PlasmaResearch/output_s/2_stream/July_25_2018/'...
%         figname '_ppc_' int2str(ppcell)...
%         '_cells_' int2str(N_mesh)];
%     if savefig
%     print(figname,'-dpng')
% end

if plot_rep
    figure(5)
%     subplot(2,1,1)
    imagesc(tlist, xlist, densitytot)
    title('Charge Density \rho v_0/(\omega_p e)')
    xlabel('t/ \omega_p')
    ylabel('x v_0/\omega_p')
    colorbar
    set(gca,'fontsize', figure_font,'YDir','normal')
    
%     subplot(2,1,2)
%     imagesc(tlist, xlist, Etot)
%     title('Electric field E e/(m_ev_0\omega_p)')
%     xlabel('t /\omega_p')
%     ylabel('x v_0/\omega_p')
%     colorbar
%     set(gca,'fontsize', figure_font,'YDir','normal')

%     subplot(3,1,3)
%     imagesc(tlist, xlist, phitot)
%     title('Electric potential \phi e/(m_ev_0^2)')
%     xlabel('t /\omega_p')
%     ylabel('x v_0/\omega_p')
%     colorbar
%     set(gca,'fontsize', figure_font,'YDir','normal')
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
    legend('Potential','kinetic','total','Location','best')
    xlabel('t \omega_p')
    set(gca,'fontsize', figure_font)
    
    figure
    subplot(2,1,1)
    vdrift = mean(velocitiestot(1:Np/nstreams,:));
    plot(tlist,vdrift)
    title('drift velocity')
    set(gca,'fontsize', figure_font)
    subplot(2,1,2)
    vsqav = mean(velocitiestot(1:Np/nstreams,:).^2) - vdrift.^2;
    plot(tlist, sqrt(vsqav))
    title('thermal velocity')
    xlabel('t \omega_p')
    set(gca,'fontsize', figure_font)
    
    figure
    
    vdrift2 = mean(velocitiestot(Np/nstreams+1:end,:));
    vsqav2 = mean(velocitiestot(Np/nstreams+1:end,:).^2) - vdrift2.^2;
    
    
    plot(tlist,potential,tlist,2*L*(vdrift.^2+vdrift2.^2),tlist,2*L*(vsqav+vsqav2))
    legend('potential','drift','thermal','Location','best')
    title('energy over time')
    xlabel('t\omega_p')
    set(gca,'fontsize', figure_font)
    
    figure
    maxdens = abs(max(densitytot - nstreams*rho0));
    subplot(2,1,1)
    plot(tlist,maxdens)
    title('max of density perturbation')
    xlabel ('t /\omega_p')
    xlim([0,tlist(end)])
    set(gca,'fontsize', figure_font)
    
    subplot(2,1,2)
    plot(tlist,log(maxdens))
    xlim([0,tlist(end)])
    title('ln(|max \delta density|)')
    xlabel ('t /\omega_p')
    
    % regression to get slope
    Xmat = [ones(size(tlist)), tlist];
    [growth_fit] = Xmat\log(maxdens)';
    hold on
    plot(tlist,growth_fit(1) + growth_fit(2)*tlist);
    hold off
    legend('ln(|max \delta \rho|)',sprintf('%.02f + %.02f t',...
        growth_fit(1),growth_fit(2)),'Location','best')
    set(gca,'fontsize', figure_font)
    
    figure
    Ebar = fft(Etot);
    Ebar = max(abs(Ebar));
        subplot(2,1,1)
    plot(tlist,Ebar)
    title('max of |fft of E|')
    xlabel ('t /\omega_p')
    xlim([0,tlist(end)])
    set(gca,'fontsize', figure_font)
    
    subplot(2,1,2)
    plot(tlist,log(Ebar))
    xlim([0,tlist(end)])
    title('ln(max|fft of E|)')
    xlabel ('t /\omega_p')
    
    % regression to get slope
    Xmat = [ones(size(tlist)), tlist];
    [growth_fit] = Xmat\log(Ebar)';
    hold on
    plot(tlist,growth_fit(1) + growth_fit(2)*tlist);
    hold off
    legend('ln(|max |fft E|)',sprintf('%.02f + %.02f t',...
        growth_fit(1),growth_fit(2)),'Location','best')
    set(gca,'fontsize', figure_font)
    
    
end

if plot_first
    figure
    plot(tlist, positionstot(1,:))
    title('Position of left-most particle over time')
    xlabel('t\omega_p')
end

if plot_parts
    figure
    plot(tlist,positionstot(1:16:end,:),'.','MarkerSize',20)
    title('particle positions')
    xlabel('t/(1/\omega_p)')
    ylabel('x/(v_0/\omega_p)')
    set(gca,'fontsize',24)
end


if plot_periodic
%     tlist = delt*(0:Nt);

    delom = 2*pi/Nt/delt;
    
    if mod(Nt,2)==0
        n_left = -Nt/2;
        n_right = Nt/2;
    else
        n_left = -ceil(Nt/2);
        n_right = floor(Nt/2);
    end
    delk = 2*pi/L;
    if mod(N_mesh,2) == 0
        nk_left = -N_mesh/2;
        nk_right = N_mesh/2-1;
    else
        nk_left = -floor(N_mesh/2);
        nk_right = floor(N_mesh/2);
    end
    figure
    subplot(2,1,1)
    plot(tlist,positionstot(1,:))
    xlabel('t\omega_p')
    ylabel('x /v_0/\omega_p')
    title(sprintf('Position of leftmost particle over time, k=%.02f',k))
    set(gca,'fontsize', figure_font)
    xlim([0,delt*Nt])
    
    subplot(2,1,2)
    omegas = delom * (n_left:n_right);
    y = fft(positionstot(1,:)-mean(positionstot(1,:)));
    m = abs(fftshift(y));
    plot(omegas,m)
    xlim([0,5])
    xlabel('\omega/\omega_p')
    ylabel('x /v_0/\omega_p')
    title('|FFT(Position of leftmost particle over time)|')
    set(gca,'fontsize', figure_font)
    
%     
%     figure
%     subplot(2,1,1)
%     plot(tlist,soln(N+1,:))
%     xlabel('t\omega_p')
%     ylabel('v /v_0')
%     title('Velocity of leftmost particle over time')
%     xlim([0,delt*Nt])
%     set(gca,'fontsize', figure_font)
%     
%     
%     
% %     figure
%     subplot(2,1,2)
%     y = fft(soln(N+1,:)-mean(soln(N+1,:)));
%     m = abs(fftshift(y));
%     plot(omegas,m)
%     xlim([0,5])
%     xlabel('\omega')
%     ylabel('v /v_0')
%     title('|FFT(Velocity of leftmost particle over time)|')
%     set(gca,'fontsize', figure_font)
%     
    
    y = fft2(Etot - mean(mean(Etot)));
    m = fftshift(abs(y));
    
    figure
    subplot(1,2,1)
    imagesc(delom*[n_left,n_right],delk*[nk_left,nk_right],m)
    ylabel('k v_0 omega_p')
    set(gca,'FontSize',figure_font, 'YDir','normal')
    subplot(1,2,2)
    imagesc(2*pi/Nt/delt*[n_left,n_right],2*pi/L*[nk_left,nk_right],m)
%     title('zoom in')
    title(sprintf('2d Fourier transform of E in x, t, with k=%.02f',k))
    xlim([-2,2])
    ylim([-3,3])
    ylabel('k v_0 \omega_p')
    xlabel('\omega / \omega_p')
    set(gca,'FontSize',22,'YDir','normal')
end

if oscillator_error
    % single particle
    
    single_vec = [positionstot(1,:) - positions0(1);...
        (velocitiestot(1,:)+[velocitiestot(1,2:end),vnew(1)])/2];
%    single_vec = [xnew(1) - alpha(1),(vnew(1)+vnewnew(1))/2];
    initial_single = [positions1(1)-positions0(1);0]*ones([1,Nt+1]);
    initial_amp = sqrt(initial_single(1,:).^2+initial_single(2,:).^2);
    single_amp = sqrt(single_vec(1,:).^2+single_vec(2,:).^2);
    
    % amplitude error
    single_rel_amp_error = abs(single_amp-initial_amp)...
        ./initial_amp;
    
    
    %phase error
    initial_phase = -atan2(initial_single(2,:),initial_single(1,:));
    final_phase = -atan2(single_vec(2,:),single_vec(1,:));
    final_phase = unwrap(final_phase);
    expected_phase = mod(tlist-initial_phase,2*pi);
    expected_phase = unwrap(expected_phase);
    single_rel_phase_error = (final_phase - expected_phase)...
        ;
    
    figure
    subplot(3,1,1)
    plot(tlist,single_amp,tlist,initial_amp)
    title('amplitude of left most particle')
    legend('computed amplitude','expected','location','best')
    set(gca,'fontsize',figure_font)
    subplot(3,1,2)
    plot(tlist,final_phase,tlist,expected_phase);
    title('phase of left-most particle')
    legend('computed phase','expected','location','best')
    set(gca,'fontsize',figure_font)
    subplot(3,1,3)
    plot(single_vec(1,:),single_vec(2,:))
    title('phase plane orbit of left-most particle')
    set(gca,'fontsize',figure_font)
  
    % wave oscillator
    wave_amp = sqrt(max(abs(phitot)).^2 + max(abs(currenttot)).^2);
    wave_phase = -atan2(phitot(1,:).*max(abs(phitot)),...
        currenttot(1,:).*max(abs(currenttot)));
    wave_phase = unwrap(wave_phase);
    
    
    figure
    subplot(3,1,1)
    plot(tlist,wave_amp,tlist,wave_amp(1)*ones(size(tlist)))
    title('amplitude of current,density wave')
    legend('computed amplitude','expected','location','best')
    set(gca,'fontsize',figure_font)
    subplot(3,1,2)
    plot(tlist,wave_phase,tlist,tlist - wave_phase(1));
    title('phase of wave oscillator')
    legend('computed phase','expected','location','best')
    set(gca,'fontsize',figure_font)
    subplot(3,1,3)
    title('phase plane orbit of wave_oscillator')
    plot(currenttot(1,:).*max(abs(currenttot)),phitot(1,:).*max(abs(phitot)))
    set(gca,'fontsize',figure_font)
    
    figure
    
%     plot(tlist,single_rel_amp_error,
    plot(tlist,single_rel_phase_error)
    title('Single particle oscillator  error')
    legend('phase','location','best')
    xlabel('t/(1/\omega_p)')
    set(gca,'fontsize',figure_font)
    tfreal = Nt*delt;
    oscillator_errors=...%[single_rel_amp_error(end),single_rel_phase_error(end),...
        ...%[(final_phase(end)-initial_phase(1))/tf,...
        [(final_phase(end)-final_phase(1))/tfreal-1,max(single_rel_amp_error)]
end

if compute_error
   % by fourier transform
   y0 = fftshift(fft(Etot(:,1)));
   m0 = abs(y0)/N_mesh;
   [max0,maxind] = max(m0);
   p0 = angle(y0);
   initial_phase = p0(maxind);
   
   
   yf = fftshift(fft(Etot(:,end)));
   mf = abs(yf)/N_mesh;
   [maxf,maxind] = max(mf);
   pf = angle(yf);
   final_phase = pf(maxind);
end


if save_movie
    close(picmovie);
end
postruntime = toc
