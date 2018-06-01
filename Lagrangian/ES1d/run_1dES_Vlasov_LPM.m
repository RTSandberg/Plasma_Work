%%% run_1dES_Vlasov_LPM
clear
close all

run_day = 'June_1';
run_name = 'convergence';
input_deck = ['./input_decks/' run_name '_input.mat'];

make_input_deck(run_day, run_name, input_deck)
load(input_deck)
tf = 2;
save(input_deck)
mytryLagrangeVlasov2(input_deck)
% 
postrun_analytics = postrun