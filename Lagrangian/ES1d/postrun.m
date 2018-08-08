function analytics = postrun
%%%
%%% post run options:
%%%     plot_dephi
%%%     plot_part
%%%     plot_two
%%%     plot_phase
%%%     save_movie
%%%     inter_particle_separation

%%% output: analytics: max force (max(E), l1(E)
load ../../../big_simulation_data/output_data
analytics = zeros([9,1]);

analytics(1:4) = [norm(Etot(:,end),Inf ), norm(Etot(:,end),1),...
    mean(soln(N+1:end,end)), std(soln(N+1:end,end))];


tlist = delt*(0:Nt);

if plot_dephi
    % plot density
    figure
    subplot(2,1,1)
    imagesc([0,input_data.tf], [input_data.xmin,input_data.xmin+input_data.L],densitytot)
    colorbar()
    title('Charge Density \cdot e\omega_p/v_0')
    xlabel('t\omega_p')
    ylabel('x v_0/\omega_p')
    set(gca,'fontsize', figure_font,'YDir','normal')

    subplot(2,1,2)
    imagesc([0,input_data.tf], [input_data.xmin,input_data.xmin+input_data.L],Etot)
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
xlim([0,input_data.tf])
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
    scatter(soln(1:N,end),soln(N+1:end,end),pointsize,f0vec);
%     plot(soln(1:2:N,end),soln(N+1:2:end,end),'b.')
%     hold on
%     plot(soln(2:2:N,end),soln(N+2:2:end,end),'r.')
%     hold off
    xlim([0,L])
    ylim([-2,2])
    title(sprintf('Phase space particles at time = %.02f',tlist(end)));
    colorbar()
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
        
    klist = 2*pi/delx/Nx*(-Nx/2:Nx/2-1);
    plot(klist,abs(fftshift(fft(E))))
    title(sprintf('Fourier spectrum of E at time = %f',tlist(end)));
    xlabel('k')
    xlim([-2,2])
        set(gca,'fontsize', figure_font)
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
