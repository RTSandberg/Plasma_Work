% script for computing and presenting pre-run diagnostics
function prerun(figure_name, save_movie, subplot_array,plot_data)%(xvec0,vvec0,f0vec,pointsize,xmin,L,figure_font,figure_name,save_movie)
%pre run diagnostics

    figure(1)
    set(gcf,'Position',[7 91 1059 894])

    for plot_info = subplot_array
        subplot(plot_info.m,plot_info.n,plot_info.p)
        eval(strcat(plot_info.plot_feature, '(plot_data,plot_info)' ) )
    end       


    
    if save_movie
        savefig([figure_name 'phase_initial'])
    end
pause(.3)
end

function phase_initial(plot_data,plot_info)
    
    
    scatter(plot_data.xvec0,plot_data.vvec0,plot_data.pointsize, plot_data.f0vec);
    if plot_info.setx
        xlim([plot_info.xlim(1),plot_info.xlim(2)])
    end
    if plot_info.sety
        ylim([plot_info.ylim(1),plot_info.ylim(2)])
    end
    title('Initial phase space particles')
    colorbar()
    xlabel('x\omega_p/v_0')
    ylabel('v/v_0')
    set(gca,'fontsize',plot_data.figure_font)
end

function interpolate_phase_space(plot_data,plot_info)


        xvis = plot_info.xlim(1):plot_info.delxvis:plot_info.xlim(2);
        vvis = plot_info.ylim(1):plot_info.delvvis:plot_info.ylim(2);
        [X,V] = meshgrid(xvis,vvis);
        xperiodic = plot_data.xvec0;
        vperiodic = plot_data.vvec0;
        xperiodic = [xperiodic - plot_data.L; xperiodic; xperiodic+plot_data.L];
        vperiodic = [vperiodic; vperiodic; vperiodic];
        f0vecperiodic = plot_data.f0vec;
        f0vecperiodic = [f0vecperiodic;f0vecperiodic;f0vecperiodic];
        
        finterp = griddata(xperiodic,...
            vperiodic,f0vecperiodic,...
            X,V);
        imagesc(plot_info.xlim,plot_info.ylim,finterp)
        title('Initial phase space')
        xlabel('x\omega_p/v_0')
        ylabel('v/v_0')
        colorbar()
        set(gca,'fontsize', plot_data.figure_font, 'YDir','normal')
end

function plot_micro_E(plot_data,plot_info)
    xvec0 = plot_data.xvec0;
    vvec0 = plot_data.vvec0;
    N = plot_data.N;
    
    fine_mesh = plot_data.xmin:plot_data.L/4000:plot_data.xmin+plot_data.L;
    Nmesh = length(fine_mesh);
    x_incl_mesh = [xvec0; vvec0; fine_mesh'; zeros(size(fine_mesh))'];
    plot_data.ode_params.Ntr = Nmesh;
    acceltot = eval(strcat(plot_data.ode_params.function, '(x_incl_mesh,plot_data.ode_params)'));
    Emacro = acceltot(N+Nmesh+1:2*N+Nmesh);
    Emicro = acceltot(2*N + Nmesh + 1:end);
    
    micro_title = '';
    macro_title = '';
    plot(-20,1)
    plot_list = [];
    if plot_info.micro
        micro_plot = plot(fine_mesh,Emicro,'.','MarkerSize',6,'DisplayName','Micro field');
        micro_title = 'micro';
        plot_list = micro_plot;
    end
    if plot_info.macro
        hold on
        macro_plot = plot(xvec0(1:N),Emacro,'o','DisplayName','Macro field');
        hold off
        macro_title = 'macro';
        plot_list = [plot_list, macro_plot];
    end
    if plot_info.macro && plot_info.micro
        micro_title = 'micro and';
    end
        
%     end
    if plot_info.setx
    xlim([plot_info.xlim(1),plot_info.xlim(2)])
    end
    if plot_info.sety
    ylim([plot_info.ylim(1),plot_info.ylim(2)])
    end
    title(sprintf('Initial %s%s field',micro_title, macro_title));
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
        macro_plot = plot(xvec0(1:N),phimacro,'o','DisplayName','Macro potential');
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

