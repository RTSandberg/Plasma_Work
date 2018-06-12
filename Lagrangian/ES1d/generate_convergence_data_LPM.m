

%%% generate_convergence_data_LPM
clear
close all

% two movies, v0 = 1.5, Nx = 6 and 128, L = 2*pi
run_day = 'June_12';
run_name = 'convergence_big_v_small_N';
% input_deck = ['./input_decks/' run_name '_input.mat'];
xmin = 0; 
L = 2*pi; Nx = 6; 
vmin = 1.4; vmax = 1.6; delv = .2; delt = .05;
% make_input_deck(run_day, run_name, input_deck,xmin,L,Nx,delv,vmin,vmax,delt)
% 
% load(input_deck)
% tf = 50;
% save_movie = 1; plot_in_run = 1; 
% beams.v = 1.5;
% save(input_deck)
% 
% mytryLagrangeVlasov2(input_deck)
% postrun;
% 
% run_name = 'convergence_big_v_big_N';
% input_deck = ['./input_decks/' run_name '_input.mat'];
% L = 2*pi; Nx = 128; 
% vmin = 1.4; vmax = 1.6; delv = .2; delt = .05; 
% 
% make_input_deck(run_day, run_name, input_deck,xmin,L,Nx,delv,vmin,vmax,delt)
% 
% load(input_deck)
% tf = 50; 
% save_movie = 1; plot_in_run = 1;
% beams.v = 1.5;
% save(input_deck)
% 
% mytryLagrangeVlasov2(input_deck)
% postrun;

% collect data
run_name = 'convergence';
input_deck = ['./input_decks/' run_name '_input.mat'];

Llist = [.5,2*pi, 100];
Nxlist = [4, 8, 16, 32, 64, 128, 200, 256];
deltlist = [1, .5, .25, .1, .05, .025, .01];
numNx = length(Nxlist);
numdelt = length(deltlist);
% collect for v = 1.5
analyticsv1p5fixt = zeros(2,numNx,6);
analyticsv1p5fixx = zeros(2,numdelt,6);
for ii = 1:2 
    L = Llist(ii+1); %Llist = [.5, 2*pi, 100], .5 is too short for v=1.5
    delt = .01;
    for jj = 1:numNx
        Nx = Nxlist(jj);
        make_input_deck(run_day, run_name, input_deck,xmin,L,Nx,delv,vmin,vmax,delt)
        
        
        mytryLagrangeVlasov2(input_deck)
        analyticsv1p5fixt(ii,jj,:)=postrun;
    end
    
    for jj = 1:numdelt
        delt = deltlist(jj);
        make_input_deck(run_day, run_name, input_deck,xmin,L,Nx,delv,vmin,vmax,delt)
        
        
        mytryLagrangeVlasov2(input_deck)
        analyticsv1p5fixx(ii,jj,:)=postrun;
    end
end

% for v = -.2
analyticsvnp2fixt = zeros(3,numNx,6);
analyticsvnp2fixx = zeros(3,numdelt,6);
for ii = 1:3 
    L = Llist(ii);
    delt = .01;
    for jj = 1:numNx
        Nx = Nxlist(jj);
        make_input_deck(run_day, run_name, input_deck,xmin,L,Nx,delv,vmin,vmax,delt)
        
        mytryLagrangeVlasov2(input_deck)
        analyticsvnp2fixt(ii,jj,:)=postrun;
    end
    for jj = 1:numdelt
        delt = deltlist(jj);
        
        make_input_deck(run_day, run_name, input_deck,xmin,L,Nx,delv,vmin,vmax,delt)
        mytryLagrangeVlasov2(input_deck)
        analyticsvnp2fixx(ii,jj,:)=postrun;
    end
end

figure_name = ['../../output_files/' run_day '/' run_name '/' run_name '_'];
data_name = [figure_name 'analytics'];
save(data_name);