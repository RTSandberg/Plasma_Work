%%% generating test files
clear
close all

input_deck = './input_decks/test_deck_input';
load(input_deck);


test_name = 'test_multiples';
spots = struct('x', {1,5}, 'v', {0,0} ,'feature_weight', {.5,.5});
% beams = struct('v',{.02}, 'vth', {0.025}, 'amplitude',{.001}, 'perturb', 'x',...
%         'wavelength',{L}, 'locationphase', {0}, 'feature_weight', {1});
beams = struct([]);


f0 = make_f0_features(input_deck,spots,beams);
test_file_name = strcat('./unit_tests/make_f0vec_features_tests/', test_name)
% save( test_file_name,'xmin','L','Nx','delv','vmax','spots','beams','f0')

