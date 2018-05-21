% function make_input_deck()
% generate input deck
% this script contains a list of quantities and saves them to a struct
% compatible with the 1d electrostatic Lagrangian Vlasov equation solver

clear
close all


run_day = 'May_18';
run_name = 'test_deck';
input_deck = ['./input_decks/' run_name '_input.mat'];

plot_initial =1;
plot_in_run =1; 
save_movie = 0; 
plot_dephi =1; 
plot_part= 1; 
plot_two = 0; 
plot_phase= 1; 
% 
tf = 4*pi;
delt = .6; 
% Nt = ceil(tf/delt);
% 
xmin = 0; 
L = 2*pi; 
Nx =8;
delx = L/Nx;
% 
delv = .06; 
vmax = .03; 
% 
% %f0vec comes one of three ways: (1) precomputed in a file, (2) by
% function to be evaluated on x,v, or (3) by selecting physical features
% most likely we'll use features
save(input_deck)

spots = struct('x',{1,3},'v',{0,0},'N',{.5,.5} );
spots = struct([]);

beams = [];
beams = struct('v',{.00}, 'vth', {0.00}, 'amplitude',{.001}, 'perturb', 's',...
        'wavelength',{L}, 'locationphase', {0}, 'n0',{1});
f0vec = make_f0_features(input_deck,spots,beams);

m = 1;
q = 1;
% 
figure_font = 22; 
pointsize =10;
% 
a =1;
% 
save(input_deck)
mytryLagrangeVlasov2(input_deck)