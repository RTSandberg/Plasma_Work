%%% LPM analysis script
topic = 'ooptest/diagnostics';
run_day = 'Dec_14_2018';
run_name = 'test1';
sim_root = ['../../../PlasmaResearch/output_s/developing_LPM/developing_features/' topic '/' run_day '/' run_name ];

Nt = 10;
fig = figure;
for ii = 0:Nt
% h5disp([sim_root '/raw_particle_data/plasma/' ...
%     sprintf('raw-plasma-%d.h5',ii)],'/set1');

% load data
data = h5read([sim_root '/raw_particle_data/plasma/' ...
    sprintf('raw-plasma-%d.h5',ii)],'/set1');
x = data(1,:); v = data(2,:);

% plot
plot(x,v)
%save
print(fig,[sim_root '/movie_images/' num2str(ii)],'-dpng');


end