%%% LPM analysis script
topic = 'ooptest/diagnostics';
run_day = 'Dec_14_2018';
run_name = 'test1';
sim_root = ['../../../PlasmaResearch/output_s/developing_LPM/developing_features/' topic '/' run_day '/' run_name ];

L = 2*pi;
Nt = 200;
tf = 4*2*pi;
dt = tf/Nt;
fig = figure;
for ii = 0:Nt
    if mod(ii,1)==0
% % h5disp([sim_root '/raw_particle_data/plasma/' ...
%     sprintf('raw-plasma-%d.h5',ii)],'/set1');

% load data
data = h5read([sim_root '/data/raw_particle_data/plasma/' ...
    sprintf('raw-plasma-%d.h5',ii)],'/set1');
Edata = h5read([sim_root '/data/raw_field_data/raw_E/' ...
    sprintf('raw-E-field-%d.h5',ii)],'/set1');
x = data(1,:); v = data(2,:);
E = Edata(2,:);
x = mod(x,L);

% plot
subplot(2,1,1)
plot(x,v,'*')
xlim([0,2*pi])
ylim([-.05,.05])
xlabel('x')
ylabel('v')

title(sprintf('phase space t=%.02f',ii*dt))
set(gca,'fontsize',20)

subplot(2,1,2)
plot(x,E,'*')
ylim([-.05,.05])
%save
% print(fig,[sim_root '/movie_images/' num2str(ii)],'-dpng');
pause(.01)
    end

end