% the next iteration in Vlasov solvers
% begin 12/13/2018
% 

% eventually, need to read input deck here!

clear
close all

% file managers
topic = 'ooptest/diagnostics';
run_day = 'Dec_14_2018';
run_name = 'test1';
sim_root = ['../../../PlasmaResearch/output_s/developing_LPM/developing_features/' topic '/' run_day '/' run_name ];


% smoothing parameter
delta = .2;

% time
tf = pi/2;
Nt = 10;
delt = tf/Nt;
ndump = 1;          % diagnostic dump frequency

% domain
xmin = 0; 
L=2*pi;

% initialize species
num_species = 1;

plasma_species = LPM_species.empty(num_species,0);
plasma_species(1) = LPM_species('plasma',-1,1,'stream',0,10,L,.2,sim_root,Nt);
% plasma_species(2)= LPM_species('plasma2',-1,1,'stream',0,10,L,.2);

% calculate rhobar
rhobar = 0;
for ii = 1:num_species
    rhobar = rhobar - sum(plasma_species(ii).weights)*plasma_species(ii).q/L;
end

plasma_fields = Plasma_fields(1,delta,rhobar,L,sim_root,Nt);

% compute initial diagnostics
    for ii = 1:num_species
        plasma_species(ii).save_diagnostics(0)
    end
    
for time_count = 1:Nt
    
    % update species
    
    time_step(plasma_species, plasma_fields, num_species, delt);
    % compute diagnostics
    for ii = 1:num_species
        plasma_species(ii).save_diagnostics(time_count)
    end
    % do amr
end

