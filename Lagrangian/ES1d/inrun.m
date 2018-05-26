function inrun(plot_in_run,save_movie, subplot_array,plot_data)
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



    if plot_in_run
        
        figure(1)
        
        for plot_info = subplot_array
            subplot(plot_info.m,plot_info.n,plot_info.p)
            eval(strcat(plot_info.plot_feature, '(plot_data,plot_info)' ) )
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
    if plot_info.setx
    xlim([plot_info.xlim(1),plot_info.xlim(2)])
    end
    if plot_info.sety
    ylim([plot_info.ylim(1),plot_info.ylim(2)])
    end
    title(sprintf('Phase space at time = %.02f',plot_data.time));
    colorbar()
    xlabel('x')
    ylabel('v')
    set(gca,'fontsize',plot_data.figure_font)
end

function plot_interpolated_phase_space(plot_data,plot_info)
        edistribution = xvweight(plot_data.x(1:plot_data.N), ...
            plot_data.x(plot_data.N+1:2*plot_data.N), plot_data.f0vec, ...
            plot_data.xmesh, plot_data.delx,plot_data.xmin,plot_info.delvvis, plot_info.ylim(2));

        imagesc([plot_info.xlim(1),plot_info.xlim(2)],[plot_info.ylim(1),plot_info.ylim(2)],edistribution)
        title('phase space')
        xlabel('x\omega_p/v_0')
        ylabel('v/v_0')
        colorbar()
        set(gca,'fontsize', plot_data.figure_font, 'YDir','normal')
end

function alt_plot_interpolated_phase_space(plot_data,plot_info)

        xvis = plot_info.xlim(1):plot_info.delxvis:plot_info.xlim(2);
        vvis = plot_info.ylim(1):plot_info.delvvis:plot_info.ylim(2);
        [X,V] = meshgrid(xvis,vvis);
        finterp = griddata(plot_data.x(1:plot_data.N),...
            plot_data.x(plot_data.N+1:end),plot_data.f0vec,...
            X,V);
        imagesc(plot_info.xlim,plot_info.ylim,finterp)
        title('phase space')
        xlabel('x\omega_p/v_0')
        ylabel('v/v_0')
        colorbar()
        set(gca,'fontsize', plot_data.figure_font, 'YDir','normal')
end


function plot_Edp(plot_data,plot_info)
        plot(plot_data.xmesh,plot_data.density,plot_data.xmesh, plot_data.E)%,xmesh,phi)
        title(sprintf('Time = %f',plot_data.time));
    if plot_info.setx
    xlim([plot_info.xlim(1),plot_info.xlim(2)])
    end
    if plot_info.sety
    ylim([plot_info.ylim(1),plot_info.ylim(2)])
    end
        legend('density','field')%,'potential')
        xlabel('x \omega_p/v_0')
        
        set(gca,'fontsize', plot_data.figure_font)
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
%     xlabel('x')
%     ylabel('v')
%     set(gca,'fontsize',figure_font)
