% script for computing and presenting pre-run diagnostics
function prerun()%(xvec0,vvec0,f0vec,pointsize,xmin,L,figure_font,figure_name,save_movie)
%pre run diagnostics

    load('output_data')
    figure(6)
    pointsize ;
    scatter(xvec0,vvec0,pointsize, f0vec);
    xlim([xmin,xmin+L])
    title('Initial phase space')
    colorbar()
    xlabel('x')
    ylabel('v')
    set(gca,'fontsize',figure_font)
    
    if save_movie
        print([figure_name 'phase_initial'],'-dpng')
    end
pause(.3)