%%%VlasovLPMCore

%%% implement a semi-lagrangian approach to solving 
%%% 1x1p collisionless electrostatic Vlasov equation.

%%% Ryan Sandberg
%%% started March 2018
%%% last updated August 9, 2018

%%% Vlasov is f' + v. gradx f + a. gradv f = 0
%%% f' + v dfdx + q/m (E) dfdv = 0

%%% Krasny's steps
%%% 1) initialize grid <-> get initial positions of particles
%%%     particles follow phase space trajectories xi(t), vi(t)
%%% 2) solve ode's: xi' = vi, vi' = F(xi) =
%%%     q^2 sum (k(xi,xj) + xj) f0(aj,bj) delta a delta b + q rhobar xi
%%% rhobar = background ion charge density (i.e. this is an electron
%%%     vlasov solver



% function VlasovLPMCore
clear
close all

tic

topic = 'LPM_eigenfunctions_cold_plasma';
run_day = 'Oct_18_2018';
run_name = 'rk4/tf_tp_Nx_32_delt_p001_n1_p6';
figure_name = ['../../../PlasmaResearch/output_s/' topic '/' run_day '/' run_name '_'];
movie_name = [figure_name 'phase_space.avi'];

Lagrangev = 1;
key_params = {};


save_movie = 0; 
save_figs = 0;
save_final = 1;

tf = 2*pi;
delt = .001; 
Nt = ceil(tf/delt); 
key_params = [key_params,'tf','delt', 'Nt'];
% 
xmin = 0; 
%L = 7.2552;
L = 2*pi;
Nx = 32;
delx = L/Nx;
% 
delv = .2;
vmin = -.5;
vmax = .5; 
key_params = [key_params,'xmin','L','Nx','delx','delv'];
convergence_name = [figure_name sprintf('Nx_%d_delt_p001_',Nx)];

m = 1;
q = -1;
key_params = [key_params,'m','q'];

% set up diagnostic mesh
Nxd = Nx/2;
 delxd = L/Nxd; 
xmesh = xmin + .5*delxd : delxd : xmin + L; xmesh = xmesh';




% pic-like initialization
k0 = 2*pi/L;
k = k0;
n0 = 1;
n1 = .6;
v0 = 0;

alpha = xmin+0*delx : delx : xmin + L-.1*delx;
xvec0 = alpha + n1/n0/k*cos(k*alpha);
xvec0 = mod([xvec0],L)';
% xvec0 = [L/8;2*L/8];
% [xvec0,sortind,indc] = unique(xvec0);
vvec0 = ones(size(alpha))';
vvec0 = 0*vvec0;
% vvec0 = v0*[vvec0;-vvec0];
% vvec0 = vvec0(sortind);
% vvec0 = [0;0];

Nv = 1;
f0vec = 1/delv*ones(size(xvec0));
f0vec = 1/delv*ones(size(xvec0));
% not using the first order approximation:
% f0vec = 1/delv*(1+n1*sin(k*xvec0)).*(1-n1*sin(k*alpha'));

rhobar = -q/L*delx*delv*sum(f0vec);
N = length(f0vec);
% key_params = [key_params,'vth','k','rho0','deln','Nv','N'];

% 
figure_font = 22; 
pointsize =16;
% 
method_params = struct('method','rk4','delt', delt, 'periodic',1,'xmin',0,'period',L,'a',1);
ode_params = struct('smooth',0);
%

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

%
% in run - phase space particles, weighted phase space, macro E, density, potential
%           micro E, density, potential
%   choose 1 - 4 options out of 
%   1) plot_phase_space_part, plots phase space coordinates
%   2) plot_interpolated_phase_space
%   3) alt_plot_interpolated_phase_space : use MATLAB interpolant
%   4) plot_Edp, plots E, density, and potential
%   5) aperiodicity,  
%   6) plot_micro_E,
%   7) periodic_plot_micro_E
diagnostic_increment = 2;
start_plot_in_run =1; 
plot_in_run = 0;


num_inrun=0; inrun_subplot_array  = struct([]);
num_inrun=num_inrun+1;
inrun_subplot_array = [inrun_subplot_array, struct('p', num_inrun, ...,
    'plot_feature', 'plot_phase_space_part',...
    'setx',1,'xlim',[xmin,xmin+L],'delxvis',delx,'sety',1,...
    'ylim',[2*vmin-.01,2*vmax+.01],'delvvis',10*delv,'micro',0,'macro',1)];
% num_inrun=num_inrun+1;
% inrun_subplot_array = [inrun_subplot_array, struct('p', num_inrun, ...,
%     'plot_feature', 'aperiodicity',...
%     'setx',1,'xlim',[xmin,xmin+delx],'delxvis',delx,'sety',1,...
%     'ylim',[5*vmin,5*vmax],'delvvis',1*delv,'micro',1,'macro',0)];
% num_inrun=num_inrun+1;
% inrun_subplot_array = [inrun_subplot_array, struct('p', num_inrun, ...,
%     'plot_feature', 'periodic_plot_micro_E',...
%     'setx',1,'xlim',[xmin,xmin+delx],'delxvis',delx,'sety',0,...
%     'ylim',[5*vmin,5*vmax],'delvvis',1*delv,'micro',1,'macro',0)];


% num_inrun=num_inrun+1;
% inrun_subplot_array = [inrun_subplot_array, struct('p', num_inrun, ...
%     'plot_feature', 'alt_plot_interpolated_phase_space',...
%     'setx',1,'xlim',[0,L],'delxvis',delx,'sety',1,...
%     'ylim',[2*vmin,2*vmax],'delvvis',.01*delv,'micro',1,'macro',0)];
% num_inrun=num_inrun+1; inrun_subplot_array = [inrun_subplot_array, struct('p',...
%     num_inrun, 'plot_feature', 'plot_micro_E',...
%     'setx',1,'xlim',[xmin,xmin+L],'delxvis',delx,'sety',1,...
%     'ylim',[-.03,.03],'delvvis',11*delv,'micro',1,'macro',1)];
% num_inrun=num_inrun+1; inrun_subplot_array = [inrun_subplot_array, struct('p',...
%     num_inrun, 'plot_feature', 'plot_micro_phi',...
%     'setx',1,'xlim',[xmin,xmin+L],'delxvis',delx,'sety',1,...
%     'ylim',[-.031,.031],'delvvis',11*delv,'micro',0,'macro',1)];
num_inrun=num_inrun+1;
inrun_subplot_array = [inrun_subplot_array,  struct('p', num_inrun, ...
    'plot_feature', 'plot_Edp','setx',1,'xlim',[xmin,xmin+L],...
    'delxvis',delx,'sety',1,'ylim',[-2,2],'delvvis',.01,'micro',1,'macro',0)];
% % num_inrun=num_inrun+1;
% inrun_subplot_array = [inrun_subplot_array,  struct('p', num_inrun, ...
%     'plot_feature', 'plot_spectrum','setx',1,'xlim',[-5,5],...
%     'delxvis',delx,'sety',0,'ylim',[0,.03],'delvvis',.01,'micro',1,'macro',0)];
% 
% num_inrun=num_inrun+1;
% inrun_subplot_array = [inrun_subplot_array,  struct('p', num_inrun, ...
%     'plot_feature', 'plot_flow_map','setx',1,'xlim',[xmin,xmin+L],...
%     'delxvis',delx,'sety',1,'ylim',[-1,1],'delvvis',.01,'micro',1,'macro',0)];
% 
% num_inrun=num_inrun+1;
% inrun_subplot_array = [inrun_subplot_array,  struct('p', num_inrun, ...
%     'plot_feature', 'plot_oscillator_phase','setx',1,'xlim',[xmin,xmin+L],...
%     'delxvis',delx,'sety',0,'ylim',[-1,1],'delvvis',.01,'micro',1,'macro',0)];

plot_rows = num_inrun; plot_cols = 1;
for ii = 1:num_inrun
    inrun_subplot_array(ii).m = plot_rows;
    inrun_subplot_array(ii).n = plot_cols;
end
%
% post run diagnostics - E, density, potential; v vs x; x vs t; 
% 2 particle case; analytic; cold case: analyticplot_initial =1;
plot_dephi =1; 
plot_part= 1; 
plot_two = 0; 
plot_phase= 0; 
inter_particle_separation = 0;
normE = 0;
plot_periodic = 0;
plot_energy = 0;
plot_pos1 = 0;
thermalization =0;
oscillator_error = 0;

need_all_data = 1;
% 
% end input parameters


soln = zeros(2*N,Nt+1);
soln(1:N,1) = xvec0; soln(N+1:end,1) = vvec0; 

%initialize diagnostics
densitytot = zeros([Nxd,Nt+1]);
currenttot = zeros([Nxd,Nt+1]);
Etot = zeros([Nxd,Nt+1]);
phitot  = zeros([Nxd,Nt+1]);

edensity = xweight(xvec0, q*f0vec, xmesh, delxd,xmin,delv);
density = edensity + rhobar;
densitytot(:,1) = density;
currentdens = xweight(xvec0,q*f0vec.*vvec0,xmesh,delxd,xmin,delv);
currenttot(:,1) = currentdens;

% one way of calculating E, phi
% X = xmesh * ones([1,Nx]);
% Y = X';
% D = X-Y;
% K = delx*(.5*sign(D) + Y/L);
% M = -delx*(.5*abs(D) + 1/L*(xvec*xvec'));

% E = K*density;
% phi = M*density;
% phi = potential_tracer([xvec0;vvec0;xvec;zeros([Nx,1])],potential_params);

[~,v] = ode_int([xvec0;vvec0],ode_params,method_params);
% E = interp1(xvec0(1:N),m/q*v(N+1:end),xmesh);
 Emesh = xweight(xvec0(1:N), m/q*v(N+1:end)/delv, xmesh, delxd,xmin,delv);
Etot(:,1) = Emesh;
% phitot(:,1) = phi(N+1:end);





plot_data = struct('pointsize',pointsize,'f0vec',f0vec,'xmin',xmin,...
    'L',L,'delt',delt, 'figure_font',figure_font,'xmesh',xmesh,...
    'N',N,'delv',delv,'delx',delx,'xvec0',xvec0,'vvec0',vvec0,'Nx',Nx,...
    'Nv',Nv,'ode_params',ode_params,'potential_params',potential_params);



    plot_data.x = [xvec0;vvec0]; plot_data.E = Emesh; plot_data.density = density;
    plot_data.current = currentdens; plot_data.time = 0; plot_data.deln = n1;
    inrun(Lagrangev,start_plot_in_run,save_movie,inrun_subplot_array,plot_data)


if save_movie
    Lagrangev = VideoWriter(movie_name);
    open(Lagrangev)
end

initialtime = toc

tic

for ii = 1:Nt
    x = soln(:,ii);
    
    [x,v] = ode_int(x,ode_params,method_params);

    soln(:,ii+1) = x;
    
    % calculate diagnostic information to save for later

    edensity = xweight(x(1:N), q*f0vec, xmesh, delxd,xmin,delv);
    density = edensity + rhobar;
    densitytot(:,ii+1) = density;
    currentdens = xweight(x(1:N),q*(x(N+1:end)).*f0vec,xmesh,delxd,xmin,delv);
    currenttot(:,ii+1) = currentdens;

%      E = interp1(x(1:N),m/q*v(N+1:end),xmesh,'pchip');
 Emesh = xweight(xvec0(1:N), m/q*v(N+1:end)/delv, xmesh, delxd,xmin,delv);
     Etot(:,ii+1) = Emesh;

    %plot diagnostics
    if ( mod(ii, diagnostic_increment) == 0) && ii >= start_plot_in_run 
        plot_data.x = x; plot_data.E = Emesh; plot_data.density = density;
        plot_data.current = currentdens;    plot_data.time = ii*delt;
        inrun(Lagrangev,plot_in_run,save_movie,inrun_subplot_array,plot_data)
    end
        
end
actualruntime = toc
tic
if save_movie
    % close panel movie writer
    close(Lagrangev);
%     savefig([figure_name 'phase_final'])
end

% if need_all_data
% save('../../../big_simulation_data/output_data')
% postrun;
% end

% having all analysis in a script at the end of the run is a 
% a terrible hack, but I don't have a better way to get everything to post
% run that doesn't involve slow file transfers.

if(save_final)
    solnfinal = soln(:,end);
    save([figure_name 'solnfinal'], 'solnfinal');
end

analytics = zeros([9,1]);

analytics(1:4) = [norm(Etot(:,end),Inf ), norm(Etot(:,end),1),...
    mean(soln(N+1:end,end)), std(soln(N+1:end,end))];


tlist = delt*(0:Nt);

if plot_dephi
    % plot density
    figure
    subplot(2,1,1)
    imagesc([0,tf], [xmin,xmin+L],densitytot)
    colorbar()
    title('Charge Density \cdot e\omega_p/v_0')
    xlabel('t\omega_p')
    ylabel('x v_0/\omega_p')
    set(gca,'fontsize', figure_font,'YDir','normal')

    subplot(2,1,2)
    imagesc([0,tf], [xmin,xmin+L],Etot)
    colorbar()
    title('E Field \cdot ')
    xlabel('t\omega_p')
    ylabel('E\epsilon_0 v_0^2/(|e|\omega_p^2)')
    set(gca,'fontsize', figure_font,'YDir','normal')

%     subplot(3,1,3)
%     imagesc([0,input_data.tf], [input_data.xmin,input_data.xmin+input_data.L],phitot)
%     colorbar()
%     title('Potential \cdot ')

    xlabel('t\omega_p')
    set(gca,'fontsize', figure_font,'YDir','normal')
    
    if save_movie   
        print([figure_name 'DE'],'-dpng')
    end
end
    
if plot_energy
    figure
    KE = delx*delv*.5*f0vec'*soln(N+1:end,:).^2;
    PE = .5*sum(Etot.^2)*delx;
    subplot(2,1,1)
    plot(tlist,KE,tlist,PE,tlist,KE+PE)
    legend('kinetic','potential','total')
    xlabel('t\omega_p')
    title('Energy')
    set(gca, 'fontsize',figure_font)
    subplot(2,1,2)
    plot(tlist,PE)
end
    
%plot particles
if plot_part
figure
plot(tlist,soln(1:N,:),'.')%,delt*1:Nt,soln(2,:),'o')
ylim([0,L])
xlim([0,tf])
title('Phase space particle positions')
xlabel('t\omega_p')
ylabel('x v_0/\omega_p')
set(gca,'fontsize', figure_font)
% % two particle test

if plot_two 

omega = q*sqrt(1/m/L);
c3 = (vvec0(2) + vvec0(1)) / 2;
c1 = (vvec0(1) - c3)/omega;
c2 = 1/2*(xvec0(1) - xvec0(2) + L/2);
c4 = xvec0(1) - c2;

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
    figure
    set(gcf, 'Position', [1 57 939 648])
%     plot(soln(1,:),soln(N+1,:),'.')
%     for n = 1:N
%         hold on
%         plot(soln(n,:),soln(N+n,:),'.')
%         hold off
%     end
%     ylim([-.005,.005])
%     xlim([0,L])
%     title('phase space orbits')
%     ylabel('v/v_0')
%     xlabel('x v_0/\omega_p')
%     set(gca,'fontsize', figure_font)
%     
    subplot(3,1,1)
%     scatter(soln(1:N,end),soln(N+1:end,end),pointsize,f0vec);
    plot(soln(1:2:N,end),soln(N+1:2:end,end),'b.')
    hold on
    plot(soln(2:2:N,end),soln(N+2:2:end,end),'r.')
    hold off
    xlim([0,L])
    ylim([-2,2])
    title(sprintf('Phase space particles at time = %.02f',tlist(end)));
    
    xlabel('x\omega_p/v_0')
    ylabel('v/v_0')
    set(gca,'fontsize',figure_font)
    
    subplot(3,1,2)
    plot(xmesh, E)%,xmesh,phi)
    title(sprintf('E at time %.02f',tlist(end)));
    xlim([0,L])
%         legend('density',
        legend('E');%,'potential')
        xlabel('x \omega_p/v_0')
        
        set(gca,'fontsize', figure_font)
        
        subplot(3,1,3)
        
    klist = L/delxd/Nxd*(-Nxd/2:Nxd/2-1);
    % to deal with E interpolations that may have fringe nans:
    idnan = isnan(E);
    E(idnan) = zeros(size(E(idnan)));
    plot(klist,abs(fftshift(fft(E))))
    title(sprintf('Fourier spectrum of E at time = %f',tlist(end)));
    xlabel('k*L/2\pi= mode number')
    xlim([-3,3])
        set(gca,'fontsize', figure_font)

    if save_figs
        print([convergence_name 'final_phase'],'-dpng')
    end

end

if inter_particle_separation
    figure
    xseparation = soln(1:N,:);
    xseparation = reshape(xseparation,[Nv,Nx,Nt+1]);
    
    xseparation = xseparation(:,2:Nx,:) - xseparation(:,1:Nx-1,:);
    xseparation = reshape(xseparation,[Nv*(Nx-1),Nt+1]);
    xseparation = L/2/pi*unwrap((xseparation)'*2*pi/L)';
    xseparation = xseparation - xseparation(:,1)*ones([1,Nt+1]);

    plot(tlist,max(abs(xseparation) ))
    title('Separation Error')
    xlabel('t\omega_p')
    ylabel('\Delta x \omega_p/v_0')
    
    set(gca,'fontsize', figure_font)
    
    analytics(5:6) = [norm(xseparation,Inf), norm(xseparation,1)];
end

if normE
    figure
    maxE = max(abs(Etot));
    plot(tlist, maxE)
    hold on
    plot(tlist,sum(abs(Etot)));
    hold off
    legend('inf norm (E)', 'l1 norm E')
    xlabel('t\omega_p')
    ylabel('E\epsilon_0 v_0^2/(|e|\omega_p^2)')
    title('Magnitudes of E over time')
    set(gca,'fontsize', figure_font)
end

if plot_periodic
    tlist = delt*(0:Nt);

    delom = 2*pi/Nt/delt;
    
    if mod(Nt,2)==0
        n_left = -Nt/2;
        n_right = Nt/2;
    else
        n_left = -ceil(Nt/2);
        n_right = floor(Nt/2);
    end
    delk = 2*pi/L;
    if mod(Nx,2) == 0
        nk_left = -Nx/2;
        nk_right = Nx/2-1;
    else
        nk_left = -floor(Nx/2);
        nk_right = floor(Nx/2);
    end
    figure
    subplot(2,1,1)
    plot(tlist,soln(1,:))
    xlabel('t\omega_p')
    ylabel('x /v_0/\omega_p')
    title(sprintf('Position of leftmost particle over time, k=%.02f',2*pi/beams.wavelength))
    set(gca,'fontsize', figure_font)
    xlim([0,delt*Nt])
    
    subplot(2,1,2)
    omegas = delom * (n_left:n_right);
    y = fft(soln(1,:)-mean(soln(1,:)));
    m = abs(fftshift(y));
    plot(omegas,m)
    xlim([0,5])
    xlabel('t\omega_p')
    ylabel('x /v_0/\omega_p')
    title('|FFT(Position of leftmost particle over time)|')
    set(gca,'fontsize', figure_font)
    
    
    figure
    subplot(2,1,1)
    plot(tlist,soln(N+1,:))
    xlabel('t\omega_p')
    ylabel('v /v_0')
    title('Velocity of leftmost particle over time')
    xlim([0,delt*Nt])
    set(gca,'fontsize', figure_font)
    
    
    
%     figure
    subplot(2,1,2)
    y = fft(soln(N+1,:)-mean(soln(N+1,:)));
    m = abs(fftshift(y));
    plot(omegas,m)
    xlim([0,5])
    xlabel('t\omega_p')
    ylabel('v /v_0')
    title('|FFT(Velocity of leftmost particle over time)|')
    set(gca,'fontsize', figure_font)
    
    
    y = fft2(Etot - mean(mean(Etot)));
    m = fftshift(abs(y));
    
    figure
    subplot(1,2,1)
    imagesc(delom*[n_left,n_right],delk*[nk_left,nk_right],m)
    ylabel('k v_0 omega_p')
    set(gca,'FontSize',22, 'YDir','normal')
    subplot(1,2,2)
    imagesc(2*pi/Nt/delt*[n_left,n_right],2*pi/L*[nk_left,nk_right],m)
%     title('zoom in')
    title(sprintf('2d Fourier transform of E in x, t, with k=%.02f',2*pi/beams.wavelength))
    xlim([-2,2])
    ylim([-3,3])
    ylabel('k v_0 \omega_p')
    xlabel('\omega / \omega_p')
    set(gca,'FontSize',22,'YDir','normal')
end

if plot_pos1
    plot(tlist,soln(1,:))
    xlabel('t\omega_p')
    ylabel('x /v_0/\omega_p')
    title(sprintf('Position of leftmost particle over time, Nx = %i, del t = %.3f', Nx, delt))
    set(gca,'FontSize',22)
end

if thermalization
    vdrift1 = delx*delv*f0vec(2:2:end)'*soln(N+2:2:end,:) / L;
    % <v>/<1> = int v f / int f
    vdrift2 = delx*delv*f0vec(1:2:end)'*soln(N+1:2:end,:) / L;
    vsq1 = delx*delv*f0vec(2:2:end)'*soln(N+2:2:end,:).^2 / L;
    vsq2 = delx*delv*f0vec(1:2:end)'*soln(N+1:2:end,:).^2 / L;
    
    driftene = .5 * 2*L * (vdrift1.^2 + vdrift2.^2);
    thermalene = .5 * 2*L * (vsq1-vdrift1.^2 + vsq2-vdrift2.^2);
    PE = .5*sum(Etot.^2)*delx;
    figure
    plot(tlist,driftene,tlist,thermalene,tlist,PE)
    title('thermal v drift energy')
    legend('drift','thermal','potential')
    xlabel('t\omega_p')
    set(gca,'fontsize',figure_font)
    
    
    maxE = max(abs(Etot));
    
    figure
    subplot(2,1,1)
    plot(tlist,maxE)
    xlabel('t\omega_p')
    title('max(|E|)')
    set(gca,'fontsize',figure_font)

    subplot(2,1,2)
    plot(tlist,log(maxE))
    % regression to get slope
    Xmat = [ones(size(tlist')), tlist'];
    [growth_fit] = Xmat\log(maxE)';
    hold on
    plot(tlist,growth_fit(1) + growth_fit(2)*tlist);
    hold off
    xlabel('t\omega_p')
    title('log(max(|E|))')
    legend('log(max(|E|))',sprintf('%.02f + %.02f t',...
        growth_fit(1),growth_fit(2)),'Location','best')
    set(gca,'fontsize',figure_font)
    
    figure
    Ebar = max(abs(fft(Etot)));
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
    Xmat = [ones(size(tlist')), tlist'];
    [growth_fit] = Xmat\log(Ebar)';
    hold on
    plot(tlist,growth_fit(1) + growth_fit(2)*tlist);
    hold off
    legend('ln(|max |fft E|)',sprintf('%.02f + %.02f t',...
        growth_fit(1),growth_fit(2)),'Location','best')
    set(gca,'fontsize', figure_font)

end

if oscillator_error
    % single particle
    single_vec = [(soln(N/2,:) - alpha(N/2))/max(soln(N/2,:)-alpha(N/2));...
        soln(N+N/2,:)/max(soln(N+N/2,:))];
%    single_vec = [xnew(1) - alpha(1),(vnew(1)+vnewnew(1))/2];
%     initial_single = [xvec0(1)-alpha(1);0]*ones([1,Nt+1]);
    initial_single = single_vec(:,1)*ones([1,Nt+1]);
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
    subplot(3,1,2)
    plot(tlist,final_phase,tlist,expected_phase);
    legend('computed phase','expected','location','best')
    subplot(3,1,3)
    plot(single_vec(1,:),single_vec(2,:))
  
    
    
    %current -v density
    %amplitude error
    wave_amp = sqrt(max(currenttot).^2 + max(densitytot).^2);
    wave_amp_init = wave_amp(1);
    rel_wave_amp_error = abs(wave_amp - wave_amp_init )/wave_amp_init ;
    % phase error
    wave_phase = -atan2(sign(currenttot(1,:)).*max(currenttot),...
        sign(densitytot(1,:)).*max(densitytot));
    wave_phase = unwrap(wave_phase);
    expected_phase = tlist +wave_phase(1);
    
    
    figure
    
    subplot(3,1,1)
    plot(tlist,wave_amp,tlist,wave_amp_init*ones(size(tlist)))
    title('amplitude of wave oscillator')
    legend('computed amplitude','expected','location','best')
    subplot(3,1,2)
    plot(tlist,wave_phase,tlist,expected_phase);
    legend('computed phase','expected','location','best')
    
    subplot(3,1,3)
    plot(sign(densitytot(1,:)).*max(densitytot) ,...
        sign(currenttot(1,:)).*max(currenttot))
    
    figure
    subplot(2,1,1)
    %plot(tlist,rel_wave_amp_error,
    plot(tlist,abs(wave_phase-expected_phase))
    legend('phase','location','best')
    title('Wave oscillator  error')
    xlabel('t/(1/\omega_p)')
    set(gca,'fontsize',figure_font)
    
    
    subplot(2,1,2)
%     plot(tlist,single_rel_amp_error,
    plot(tlist,single_rel_phase_error)
    title('Single particle oscillator  error')
    legend('phase','location','best')
    xlabel('t/(1/\omega_p)')
    set(gca,'fontsize',figure_font)
    tfreal = Nt*delt;
    oscillator_errors=[single_rel_amp_error(end),...%single_rel_phase_error(end),...
        ...%[(final_phase(end)-initial_phase(1))/tf,...
        (final_phase(end)-final_phase(1))/tfreal-1,...%(wave_phase(end)-wave_phase(1))/tf,...
        (wave_phase(end)-wave_phase(1))/tfreal-1,]
end


if 0
   % anytime we save a figure or movie, we also want a file with key parameters
   filename = [figure_name 'key_parameters.txt'];
   fileID = fopen(filename,'w');
   for ii = 1:length(key_params)
       param_name = key_params{ii};
       fprintf(fileID,[param_name ' = %.02f\n'], eval(param_name) );
   end
   fclose(fileID);
end



postrun = toc