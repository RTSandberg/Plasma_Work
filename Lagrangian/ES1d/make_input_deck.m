function make_input_deck(run_day, run_name, input_deck,xmin,L,Nx,delv,vmin,vmax,delt,tf)
% generate input deck
% this script contains a list of quantities and saves them to a struct
% compatible with the 1d electrostatic Lagrangian Vlasov equation solver





save_movie = 0; 

% tf = 100;
% delt = 1; 
% 
% xmin = 0; 
% L = 2*pi; 
% Nx =32;
delx = L/Nx;
% 
% delv = .2;
% vmin = 0;
% vmax = .2; 
% 
% %f0vec comes one of three ways: (1) precomputed in a file, (2) by
% function to be evaluated on x,v, or (3) by selecting physical features
% most likely we'll use features
save(input_deck)

shear.shearQ = 0; shear.slope = .2*delx / delv;

spots = struct('x',{1},'v',{0},'N',{.5} );
spots(2) = struct('x',{3},'v',{0},'N',{.5} );
spots(3) = struct('x',{4},'v',{0},'N',{.5} );
spots = struct([]);

beams = [];
beams = struct('v',{.1}, 'vth', {.00}, 'amplitude',{.02}, 'perturb', {'n'},...
        'wavelength',{L}, 'locationphase', {0}, 'n0',{1});
% beams(2) = struct('v',{-.1}, 'vth', {0.0000}, 'amplitude',{.01}, 'perturb', {'n'},...
%         'wavelength',{L}, 'locationphase', {0}, 'n0',{.5});
f0vec = make_f0_features(input_deck,spots,beams,shear);

m = 1;
q = -1;
% 
figure_font = 22; 
pointsize =10;
% 
method_params = struct('method','rk2','delt', delt, 'periodic',1,'xmin',0,'period',L,'a',1);
ode_params = struct('smooth',0);
%

%what diagnostics?
% prerun - phase space initialization, initial field not supported right
%   choose 1 - 2 options out of 
%   1) phase_initial, plots phase space coordinates
%   2) interpolate_phase_space
%   3) plot_micro_E
num_prerun=0; prerun_subplot_array  = struct([]);
num_prerun=num_prerun+1;
prerun_subplot_array = [prerun_subplot_array, struct('p', num_prerun, ...,
    'plot_feature', 'phase_initial',...
    'setx',1,'xlim',[xmin,xmin+L],'delxvis',delx,'sety',1,...
    'ylim',[vmin,vmax],'delvvis',1*delv,'micro',0,'macro',0)];
% num_prerun=num_prerun+1;
% prerun_subplot_array = [prerun_subplot_array, struct('p', num_prerun, ...,
%     'plot_feature', 'interpolate_phase_space',...
%     'setx',1,'xlim',[xmin,xmin+L],'delxvis',delx,'sety',1,...
%     'ylim',[vmin,vmax],'delvvis',.01*delv,'micro',0,'macro',0)];
num_prerun=num_prerun+1;
prerun_subplot_array = [prerun_subplot_array, struct('p', num_prerun, ...,
    'plot_feature', 'plot_micro_E',...
    'setx',1,'xlim',[xmin,xmin+L],'delxvis',delx,'sety',1,...
    'ylim',[-.03,.03],'delvvis',.001,'micro',0,'macro',1)];
plot_rows = num_prerun; plot_cols = 1;
for ii = 1:num_prerun
    prerun_subplot_array(ii).m = plot_rows;
    prerun_subplot_array(ii).n = plot_cols;
end


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
plot_in_run =0;
diagnostic_increment = 50;

num_inrun=0; inrun_subplot_array  = struct([]);
num_inrun=num_inrun+1;
inrun_subplot_array = [inrun_subplot_array, struct('p', num_inrun, ...,
    'plot_feature', 'plot_phase_space_part',...
    'setx',1,'xlim',[xmin,xmin+L],'delxvis',delx,'sety',1,...
    'ylim',[vmin,vmax],'delvvis',1*delv,'micro',0,'macro',1)];
% num_inrun=num_inrun+1;
% subplot_array = [subplot_array, struct('p', num_inrun, ...,
%     'plot_feature', 'aperiodicity',...
%     'setx',1,'xlim',[xmin,xmin+delx],'delxvis',delx,'sety',1,...
%     'ylim',[25*vmin,25*vmax],'delvvis',1*delv,'micro',1,'macro',0)];
% num_inrun=num_inrun+1;
% subplot_array = [subplot_array, struct('p', num_inrun, ...,
%     'plot_feature', 'periodic_plot_micro_E',...
%     'setx',1,'xlim',[xmin,xmin+delx],'delxvis',delx,'sety',0,...
%     'ylim',[5*vmin,5*vmax],'delvvis',1*delv,'micro',1,'macro',0)];


% num_inrun=num_inrun+1;
% inrun_subplot_array = [inrun_subplot_array, struct('p', num_inrun, ...
%     'plot_feature', 'alt_plot_interpolated_phase_space',...
%     'setx',1,'xlim',[0,L],'delxvis',delx,'sety',1,...
%     'ylim',[vmin,vmax],'delvvis',.01*delv,'micro',1,'macro',0)];
num_inrun=num_inrun+1; inrun_subplot_array = [inrun_subplot_array, struct('p',...
    num_inrun, 'plot_feature', 'plot_micro_E',...
    'setx',1,'xlim',[xmin,xmin+L],'delxvis',delx,'sety',1,...
    'ylim',[-.025,.025],'delvvis',11*delv,'micro',0,'macro',1)];
% num_inrun=num_inrun+1;
% subplot_array = [subplot_array,  struct('p', num_inrun, ...
%     'plot_feature', 'plot_Edp','setx',1,'xlim',[xmin,xmin+L],...
%     'delxvis',delx,'sety',1,'ylim',[-.03,.03],'delvvis',.01,'micro',1,'macro',0)];

plot_rows = num_inrun; plot_cols = 1;
for ii = 1:num_inrun
    inrun_subplot_array(ii).m = plot_rows;
    inrun_subplot_array(ii).n = plot_cols;
end


% post run - E, density, potential; v vs x; x vs t; 2 particle case:
% analytic; cold case: analyticplot_initial =1;
plot_dephi =1; 
plot_part= 0; 
plot_two = 0; 
plot_phase= 0; 
inter_particle_separation = 1;
normE = 1;
% 
% 
save(input_deck)