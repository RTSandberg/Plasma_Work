%%% convergence_single_slice


%%% generate_convergence_data_LPM
clear
close all

% two movies, v0 = 1.5, Nx = 6 and 128, L = 2*pi
run_day = 'June_12';
% input_deck = ['./input_decks/' run_name '_input.mat'];
xmin = 0; 


% collect data
run_name = 'convergence';
input_deck = ['./input_decks/' run_name '_input.mat'];

% Llist = [.5,2*pi, 100];
% deltlist = [1, .5, .25, .1, .05, .025, .01];

% collect for v = 1.5, L = 2pi, delt = .01
v = 1.5; vmin = 1.4; vmax = 1.6; delv = .2;
L = 2*pi; 
delt = .01; tf = 100;

Nxlist = [4, 8, 16, 32, 64, 128, 200, 256];
numNx = length(Nxlist);
% numdelt = length(deltlist);

analytics = zeros(numNx,6);

for jj = 1:numNx

    Nx = Nxlist(jj);
    make_input_deck(run_day, run_name, input_deck,xmin,L,Nx,delv,vmin,vmax,delt,tf)


    mytryLagrangeVlasov2(input_deck)
    analytics(jj,:)=postrun;
end
    