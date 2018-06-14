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

% collect for v = -.2, L = .5, delt = .1
v = -.2; vmin = 1.4; vmax = 1.6; delv = .2;
L = 0.5; 
delt = .1; tf = 100;

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

figure
loglog(L./Nxlist,analytics(:,1),'.','MarkerSize',40)
hold on
for ii = [2,4,5,6]
    loglog(L./Nxlist, analytics(:,ii),'.','MarkerSize',40-2*ii)
end
hold off
legend('|E|_{\infty}','|E|_1','velocity \sigma','|separation|_{\infty}','|separation|_1')
title(sprintf('Cold plasma convergence, L=%.2f, delt=%.2f, v=%.2f',L,delt,v));
xlabel('\Delta x \omega_p/v_0')
set(gca,'fontsize',22)





    