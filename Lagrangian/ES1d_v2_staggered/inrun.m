function inrun(Lagrangev,plot_in_run,save_movie, subplot_array,plot_data,inrun_window)
%%% compute and display in-run diagnostics
%%% input:
%%%     plot_in_run
%%%     plot_info: struct array of subplot parameters
%%%         m: number of rows
%%%         n: number of columns
%%%         p: number in subplot array
%%%         plot_feature: string, indicates which diagnostic function to 
%%%                         call in order to plotplot
%%%     plot_data:  struct of things to use in plotting
%%%             fields:
%%%                 x, pointsize,f0vec,xmin,L,delt, figure_font,xmesh,density,E,N,time,delv
%%%
%%%         options for plot_feature:
%%%             -plot_phase_space_part,  plots points in phase space,
%%%             colored by f value
%%%             -aperiodicity,  
%%%             -plot_interpolated_phase_space
%%%             -alt_plot_interpolated_phase_space,     uses built-in
%%%             interpolation function
%%%             - plot_Edp
%%%             - plot_micro_E,
%%%             - periodic_plot_micro_E



    if plot_in_run
        
        figure(inrun_window)
        set(gcf,'Position',[7 91 1059 894])
        
        for plot_info = subplot_array
            subplot(plot_info.m,plot_info.n,plot_info.p)
            eval(strcat(plot_info.plot_feature, '(plot_data,plot_info)' ) )
            set(gca,'fontsize',plot_data.figure_font);
        end       
    
    if save_movie
        frame = getframe(gcf);
        writeVideo(Lagrangev,frame)
    end
% 
    end

end

% in run - phase space particles, weighted phase space, macro E, density, potential
%           micro E, density, potential
%   choose 1 - 4 options

function plot_phase_space_part(plot_data,plot_info)
    scatter(plot_data.x(1:plot_data.N),plot_data.x(plot_data.N+1:end),...
        plot_data.pointsize,plot_data.f0vec);
%     plot(plot_data.x(1:2:plot_data.N),plot_data.x(plot_data.N+1:2:end),'b.','MarkerSize',plot_data.pointsize)
%     hold on
%     plot(plot_data.x(2:2:plot_data.N),plot_data.x(plot_data.N+2:2:end),'r.','MarkerSize',plot_data.pointsize)
%     hold off
    if plot_info.setx
    xlim([plot_info.xlim(1),plot_info.xlim(2)])
    end
    if plot_info.sety
    ylim([plot_info.ylim(1),plot_info.ylim(2)])
    end
    title(sprintf('Phase space particles at time = %.02f',plot_data.time));
%     c = colorbar();
%     c.Label.String = 'particle weight';
    xlabel('x\omega_p/v_0')
    ylabel('v/v_0')
    set(gca,'fontsize',plot_data.figure_font)
end
function plot_interpolated_phase_space(plot_data,plot_info)
        edistribution = xvweight(plot_data.x(1:plot_data.N), ...
            plot_data.x(plot_data.N+1:2*plot_data.N), plot_data.f0vec, ...
            plot_data.xmesh, plot_data.delx,plot_data.xmin,plot_info.delvvis, plot_info.ylim(2));

        imagesc([plot_info.xlim(1),plot_info.xlim(2)],[plot_info.ylim(1),plot_info.ylim(2)],edistribution)
        title('Phase space')
        xlabel('x\omega_p/v_0')
        ylabel('v/v_0')
        colorbar()
        set(gca,'fontsize', plot_data.figure_font, 'YDir','normal')
end

function alt_plot_interpolated_phase_space(plot_data,plot_info)

        xvis = plot_info.xlim(1):plot_info.delxvis:plot_info.xlim(2);
        vvis = plot_info.ylim(1):plot_info.delvvis:plot_info.ylim(2);
        [X,V] = meshgrid(xvis,vvis);
        xperiodic = plot_data.x(1:plot_data.N);
        vperiodic = plot_data.x(plot_data.N+1:end);
        xperiodic = [xperiodic - plot_data.L; xperiodic; xperiodic+plot_data.L];
        vperiodic = [vperiodic; vperiodic; vperiodic];
        f0vecperiodic = plot_data.f0vec;
        f0vecperiodic = [f0vecperiodic;f0vecperiodic;f0vecperiodic];
        
        finterp = griddata(xperiodic,...
            vperiodic,f0vecperiodic,...
            X,V);
        imagesc(plot_info.xlim,plot_info.ylim,finterp)
        title('Phase space')
        xlabel('x\omega_p/v_0')
        ylabel('v/v_0')
        colorbar()
        caxis([min(plot_data.f0vec-.001),max(plot_data.f0vec+.001)]);
        set(gca,'fontsize', plot_data.figure_font, 'YDir','normal')
end


function plot_Edp(plot_data,plot_info)
    plot(plot_data.xmesh,plot_data.density,plot_data.xmesh,...
        plot_data.current,plot_data.x(1:plot_data.N),plot_data.E,'o')
%       plot( plot_data.xmesh,plot_data.density,plot_data.xmesh, plot_data.E,'o')%,xmesh,phi)
        title(sprintf('fields at time %.02f',plot_data.time));
    if plot_info.setx
    xlim([plot_info.xlim(1),plot_info.xlim(2)])
    end
    if plot_info.sety
    ylim([plot_info.ylim(1),plot_info.ylim(2)])
    end
%         legend('density',
        legend('density','current','E');%,'potential')
        xlabel('x \omega_p/v_0')
        
end

function plot_micro_E(plot_data,plot_info)
    x = plot_data.x;
    N = plot_data.N;
    
    fine_mesh = plot_data.xmin:plot_data.L/4000:plot_data.xmin+plot_data.L;
    Nmesh = length(fine_mesh);
    x_incl_mesh = [x; fine_mesh'; zeros(size(fine_mesh))'];
    plot_data.ode_params.Ntr = Nmesh;
    acceltot = eval(strcat(plot_data.ode_params.function, '(x_incl_mesh,plot_data.ode_params)'));
    Emacro = acceltot(N+Nmesh+1:2*N+Nmesh);
    Emicro = acceltot(2*N + Nmesh + 1:end);
    
    micro_title = '';
    macro_title = '';
    plot_list = [];
    plot(plot_data.xmin-20,0);
    if plot_info.micro
        micro_plot = plot(fine_mesh,Emicro,'.','MarkerSize',6,'DisplayName','Micro field');
        micro_title = 'micro';
        plot_list = micro_plot;
    end
    if plot_info.macro
        hold on
        macro_plot = plot(x(1:N),Emacro,'o','DisplayName','Macro field');
        hold off
        macro_title = 'macro';
        plot_list = [plot_list, macro_plot];
    end
    if plot_info.macro && plot_info.micro
        micro_title = 'micro and ';
    end
        
%     end

    

    if plot_info.setx
    xlim([plot_info.xlim(1),plot_info.xlim(2)])
    end
    if plot_info.sety
    ylim([plot_info.ylim(1),plot_info.ylim(2)])
    end
    title(sprintf('Time t=%.02f, %s%s field',plot_data.time,micro_title, macro_title));
    legend(plot_list)
    xlabel('x\omega_p/v_0')
    ylabel('E\epsilon_0 v_0^2/(|e|\omega_p^2)')
end

function plot_micro_phi(plot_data,plot_info)
    xvec0 = plot_data.xvec0;
    vvec0 = plot_data.vvec0;
    N = plot_data.N;
    
    fine_mesh = plot_data.xmin:plot_data.L/4000:plot_data.xmin+plot_data.L;
    Nmesh = length(fine_mesh);
    x_incl_mesh = [xvec0; vvec0; fine_mesh'; zeros(size(fine_mesh))'];
    plot_data.potential_params.Ntr = Nmesh;
    
    phi = eval(strcat(plot_data.potential_params.function, '(x_incl_mesh,plot_data.potential_params)'));
    phimacro = phi(1:N);
    phimicro = phi(N+1:end);
    
    micro_title = '';
    macro_title = '';
    plot(-20,1)
    plot_list = [];
    if plot_info.micro
        micro_plot = plot(fine_mesh,phimicro,'.','MarkerSize',6,'DisplayName','Micro \phi');
        micro_title = 'micro';
        plot_list = micro_plot;
    end
    if plot_info.macro
        hold on
        macro_plot = plot(xvec0(1:N),phimacro,'o','DisplayName','Macro \phi');
        hold off
        macro_title = 'macro';
        plot_list = [plot_list, macro_plot];
    end
    if plot_info.macro && plot_info.micro
        micro_title = 'micro and ';
    end
        
%     end

    if plot_info.setx
    xlim([plot_info.xlim(1),plot_info.xlim(2)])
    end
    if plot_info.sety
    ylim([plot_info.ylim(1),plot_info.ylim(2)])
    end
    title(sprintf('Initial %s%s potential',micro_title, macro_title));
    legend(plot_list)
    xlabel('x\omega_p/v_0')
%     ylabel('E\epsilon_0 v_0^2/(|e|\omega_p^2)')
    
end
function plot_spectrum(plot_data,plot_info)
    
    Nx = plot_data.Nx;
    L = plot_data.L;
    klist = 2*pi/plot_data.delx/Nx*(-Nx/2:Nx/2-1);
    % to deal with E interpolations that may have fringe nans:
    E = plot_data.Emesh;
    idnan = isnan(E);
    E(idnan) = zeros(size(E(idnan)));
    plot(klist,abs(fftshift(fft(plot_data.density)))/Nx)
    if plot_info.setx
    xlim([plot_info.xlim(1),plot_info.xlim(2)])
    end
    title(sprintf('Fourier spectrum of density at time = %.02f',plot_data.time));
    xlabel('k')
end

function plot_flow_map(plot_data,plot_info)
    % specific to 2 stream

    x = plot_data.x;
    N = plot_data.N;
    delx = plot_data.delx;
    L = plot_data.L;
    x1 = x(1:2:N);
    x2 = x(2:2:N);
    alphas = .5*delx:delx:L;
%     plot(alphas,x1,'b.','MarkerSize',plot_data.pointsize)
%     hold on
%     plot(alphas,x2,'r.','MarkerSize',plot_data.pointsize)
%     hold off
    
    v1 = x(N+1:2:end);
    v2 = x(N+2:2:end);
    x1 = [x1(end);x1;x1(1)];
    x2 = [x2(end);x2;x2(1)];
    v1 = [v1(end);v1;v1(1)];
    v2 = [v2(end);v2;v2(1)];
    flowmapcheck1 = norm([x1(3:end) - x1(1:end-2),v1(3:end)-v1(1:end-2)]);
    flowmapcheck2 = norm([x2(3:end) - x2(1:end-2),v2(3:end)-v2(1:end-2)]);
    plot(alphas,flowmapcheck1,'b.',alphas,flowmapcheck2,'r.')
title('flow map')
xlim([0,L])
% ylim([0,plot_data.L])
xlabel('\alpha')
ylabel('x \omega_p/v_0')

        set(gca,'fontsize', plot_data.figure_font)
end

function plot_oscillator_phase(plot_data,plot_info)
    plot(sign(plot_data.density(1))*max(plot_data.density),...
        sign(plot_data.current(1))*max(plot_data.current),'*')
    hold on
    plot(-plot_data.deln*cos(plot_data.time),plot_data.deln*sin(plot_data.time),'r*')
    plot(plot_data.deln*cos(0:.1:2*pi),plot_data.deln*sin(0:.1:2*pi) )
    hold off
    xlim([-plot_data.deln*(1.5),plot_data.deln*(1.5)])
    ylim([-plot_data.deln*(1.5),plot_data.deln*(1.5)])
    xlabel('sign(\rho(1))*(max(\rho)')
    ylabel('sign(j(1))*max(j)')
    title('current vs density oscillation')
    hold on
    single_amp = abs(plot_data.xvec0-.37*plot_data.delx);
    plot(plot_data.x(1)-.37*plot_data.delx ,plot_data.x(plot_data.N+1),'o')
    phase_init = atan2(0,plot_data.xvec0(1)-.37*plot_data.delx);
    plot(single_amp*cos(plot_data.time-phase_init),-single_amp*sin(plot_data.time-phase_init),'ro')
    hold off
    legend('wave phase','expected wave','single phase','expected single')
end

%other
% single particles

%         plot(x(1:floor(N/10):N),x(N+1:floor(N/10):end),'o' )
%         
%         title('Single particles in phase space')
%         xlim([0,L])
%         xlabel('x\omega_p/v_0')
%         ylabel('v/v_0')
%         ylim([-.1,.1])
%         set(gca,'fontsize', figure_font)

%zoom in:
%             scatter(x(1:N),x(N+1:end),pointsize, f0vec);
%     xlim([2,2.5])
%     title('Phase space zoom in')
%     colorbar()
%     xlabel('x\omega_p/v_0')
%     ylabel('v/v_0')
%     set(gca,'fontsize',figure_font)
