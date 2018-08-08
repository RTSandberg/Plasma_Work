function inrun(Lagrangev,plot_in_run,save_movie, subplot_array,plot_data)
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
        
        figure(2)
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
        plot_data.pointsize, plot_data.f0vec);
% %     plot(plot_data.x(1:2:plot_data.N),plot_data.x(plot_data.N+1:2:end),'b.')
% %     hold on
% %     plot(plot_data.x(2:2:plot_data.N),plot_data.x(plot_data.N+2:2:end),'r.')
%     hold off
    if plot_info.setx
    xlim([plot_info.xlim(1),plot_info.xlim(2)])
    end
    if plot_info.sety
    ylim([plot_info.ylim(1),plot_info.ylim(2)])
    end
    title(sprintf('Phase space particles at time = %.02f',plot_data.time));
    colorbar()
    xlabel('x\omega_p/v_0')
    ylabel('v/v_0')
    set(gca,'fontsize',plot_data.figure_font)
end

function aperiodicity(plot_data,plot_info)

    index_list_ii =  1 : plot_data.Nv;
    xtemp = mod(plot_data.x(index_list_ii)-plot_data.xmin,plot_data.delx);
    plot(xtemp,plot_data.x(plot_data.N + index_list_ii),'.','MarkerSize',12)
    
    for ii = 2:plot_data.Nx
        index_list_ii = (ii-1)*plot_data.Nv + 1 : ii*plot_data.Nv;
        xtemp = mod(plot_data.x(index_list_ii)-plot_data.xmin,plot_data.delx);
        hold on
        plot(xtemp,plot_data.x(plot_data.N + index_list_ii),'.','MarkerSize',12)
        hold off
    end
    if plot_info.setx
    xlim([plot_info.xlim(1),plot_info.xlim(2)])
    end
    if plot_info.sety
    ylim([plot_info.ylim(1),plot_info.ylim(2)])
    end
    title(sprintf('Quotient Phase space at time = %.02f',plot_data.time));
    legend('1st period','2nd period','3rd period','4th period')
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
%       plot(plot_data.xmesh,plot_data.density, 
      plot(plot_data.xmesh, plot_data.E)%,xmesh,phi)
        title(sprintf('fields at time %.02f',plot_data.time));
    if plot_info.setx
    xlim([plot_info.xlim(1),plot_info.xlim(2)])
    end
    if plot_info.sety
    ylim([plot_info.ylim(1),plot_info.ylim(2)])
    end
%         legend('density',
        legend('field');%,'potential')
        xlabel('x \omega_p/v_0')
        
        set(gca,'fontsize', plot_data.figure_font)
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
    set(gca,'fontsize',plot_data.figure_font)
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
    set(gca,'fontsize',plot_data.figure_font)
end

function periodic_plot_micro_E(plot_data,plot_info)
    x = plot_data.x;
    N = plot_data.N;
    
    fine_mesh = plot_data.xmin:plot_data.L/10000:plot_data.xmin+plot_data.L;
    Nmesh = length(fine_mesh);
    x_incl_mesh = [x; fine_mesh'; zeros(size(fine_mesh))'];
    plot_data.ode_params.Ntr = Nmesh;
    acceltot = eval(strcat(plot_data.ode_params.function, '(x_incl_mesh,plot_data.ode_params)'));
    Emacro = acceltot(N+Nmesh+1:2*N+Nmesh);
    Emicro = acceltot(2*N + Nmesh + 1:end);
    
%     if plot_info.micro
    delx = plot_data.delx;
        plot(mod(fine_mesh(1:1000),delx),Emicro(1:1000),'.','MarkerSize',6)
        hold on
        plot(mod(x(1:N/4),delx),Emacro(1:N/4),'o')
        hold off
        
        for ii = 2:4
            hold on            
            plot( mod(fine_mesh(1000*(ii-1)+1:1000*ii),delx),Emicro((ii-1)*1000+1:1000*ii),'.','MarkerSize',6)        
            plot(mod(x((ii-1)*N/4+1:ii*N/4),delx),Emacro((ii-1)*N/4+1:ii*N/4),'o')
            hold off
        end
%     end

    if plot_info.setx
    xlim([plot_info.xlim(1),plot_info.xlim(2)])
    end
    if plot_info.sety
    ylim([plot_info.ylim(1),plot_info.ylim(2)])
    end
    title(sprintf('%.02f spatial period of micro and macro fields at time = %.02f',plot_data.Nx/4, plot_data.time));
%     legend('E micro 1st period', 'E macro 1st period','E micro 2nd period',...
%         'E macro 2nd period','E micro 3rd period','E macro 3rd period', ...
%         'E micro 4th period','E macro 4th period')
    xlabel('x\omega_p/v_0')
    ylabel('E\epsilon_0 v_0^2/(|e|\omega_p^2)')
    set(gca,'fontsize',plot_data.figure_font)
end

function plot_spectrum(plot_data,plot_info)
    
    Nx = plot_data.Nx;
    L = plot_data.L;
    klist = 2*pi/plot_data.delx/Nx*(-Nx/2:Nx/2-1);
    plot(klist,abs(fftshift(fft(plot_data.E))))
    title(sprintf('Fourier spectrum of E at time = %f',plot_data.time));
    xlabel('k')
    
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
