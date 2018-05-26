% function make_input_deck()
% generate input deck
% this script contains a list of quantities and saves them to a struct
% compatible with the 1d electrostatic Lagrangian Vlasov equation solver

clear
close all


run_day = 'May_25';
run_name = 'convergence_warm_equil';
input_deck = ['./input_decks/' run_name '_input.mat'];


tf = 3;
delt = .1; 
% 
xmin = 0; 
L = 2*pi; 
Nx =20;
delx = L/Nx;
% 
delv = .0002; 
vmax = .001; 
% 
% %f0vec comes one of three ways: (1) precomputed in a file, (2) by
% function to be evaluated on x,v, or (3) by selecting physical features
% most likely we'll use features
save(input_deck)

spots = struct('x',{1},'v',{0},'N',{.5} );
spots(2) = struct('x',{3},'v',{0},'N',{.5} );
spots(3) = struct('x',{4},'v',{0},'N',{.5} );
spots = struct([]);

beams = [];
beams = struct('v',{0}, 'vth', {0.003}, 'amplitude',{.001}, 'perturb', {'n'},...
        'wavelength',{L}, 'locationphase', {0}, 'n0',{1});
f0vec = make_f0_features(input_deck,spots,beams);

m = 1;
q = -1;
% 
figure_font = 22; 
pointsize =10;
% 
method_params = struct('method','rk2','delt', delt, 'periodic',1,'xmin',0,'period',L,'a',1);
ode_params = struct();
%

%what diagnostics?
% prerun - phase space initialization, initial field not supported right
% now

% in run - phase space particles, weighted phase space, macro E, density, potential
%           micro E, density, potential
%   choose 1 - 4 options out of 
%   1) plot_phase_space_part, plots phase space coordinates
%   2) plot_interpolated_phase_space
%   3) alt_plot_interpolated_phase_space : use MATLAB interpolant
%   4) plot_Edp, plots E, density, and potential
plot_in_run =1; 
save_movie = 0; 
plot_rows = 2; plot_col = 1;
p=1;
subplot_array(p) = struct('m', plot_rows, 'n', plot_col, 'p', p, ...
    'plot_feature', 'plot_phase_space_part','setx',1,'xlim',[xmin,xmin+L],...
    'delxvis',delx,'sety',1,'ylim',[-110*vmax,110*vmax],'delvvis',110*delv);
p=p+1;
subplot_array(p) = struct('m', plot_rows, 'n', plot_col, 'p', p, ...
    'plot_feature', 'alt_plot_interpolated_phase_space','setx',1,'xlim',[xmin,xmin+L],...
    'delxvis',delx,'sety',1,'ylim',[-110*vmax,110*vmax],'delvvis',11*delv);
% p=p+1;
% subplot_array(p) = struct('m', plot_rows, 'n', plot_col, 'p', p, ...
%     'plot_feature', 'plot_Edp','setx',1,'xlim',[xmin,xmin+L],...
%     'delxvis',delx,'sety',0,'ylim',[-.001,.001],'delvvis',.001);
% post run - E, density, potential; v vs x; x vs t; 2 particle case:
% analytic; cold case: analyticplot_initial =1;
plot_dephi =0; 
plot_part= 0; 
plot_two = 0; 
plot_phase= 0; 
% 
% 
save(input_deck)
mytryLagrangeVlasov2(input_deck)