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
load output_data
analytics = zeros([6,1]);
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
% 
%     subplot(3,1,3)
%     imagesc([0,tf], [xmin,xmin+L],phitot)
%     colorbar()
%     title('Potential \cdot ')
% 
%     xlabel('t\omega_p')
    set(gca,'fontsize', figure_font,'YDir','normal')
    
    if save_movie   
        print([figure_name 'DE'],'-dpng')
    end
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