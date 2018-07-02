%%% get_dispersion


%%% generate_convergence_data_LPM
clear
close all

run_day = 'June_12';
topic_name = 'dispersion';
% input_deck = ['./input_decks/' topic_name '_' run_name '_input.mat'];
xmin = 0; 


% collect data
run_name = 'dispersion';
input_deck = ['./input_decks/' topic_name '_' run_name '_input.mat'];

% Llist = [.5,2*pi, 100];
% deltlist = [1, .5, .25, .1, .05, .025, .01];

% collect for v = -.2, L = .5, delt = .1
% v = -.2; vmin = -.3; vmax = -.1; delv = .2;
delv = .02;
vmin = -.2;
vmax = .2; 

L = 2*pi; 
delt = .1; tf = 50;

Nx = 100;

wavelengths = L* [2,1];%,.5,1/3];%,1/4];%,1/5,1/6,1/7,1/10,1/20,1/30,1/40,1/50];
% wavelengths = [L];
numk = length(wavelengths);

omegas = zeros(numk,1);

for jj = 1:numk
    wavelength = wavelengths(jj);
%     delt = deltlist(jj);
    make_input_deck(run_day, topic_name, run_name, input_deck,xmin,L,Nx,delv,vmin,vmax,delt,tf,wavelength)


    mytryLagrangeVlasov2(input_deck)
    analytics=postrun;
    omegas(jj) = analytics(3);
end

    