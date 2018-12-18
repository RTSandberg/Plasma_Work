%%% LPM global field evaluating functions
function time_step(species_array, plasma_fields, num_species, delt)
    % Input:
    %       species array


    % gather -> one x, weights
    [totalx, totalweight] = gather_points(species_array, num_species);
    % calculate E -> Etotal
    E = plasma_fields.calc_E(totalx, totalx, totalweight);
    % scatter E -> 
    % update individual species
    % for each species, take update x with v, update v with section of E
    particle_count = 0;
    for ii = 1:num_species
        Ntemp = species_array(ii).Npoints;
        species_array(ii).update(species_array(ii).v,...
            E(particle_count+1:particle_count + Ntemp), delt);
        particle_count = particle_count + Ntemp;
    end

    function [totalx, totalweight] = gather_points(species_array, num_species)
        % eventually do separately positive and negative charges
        % gather all x, weights
        num_parts = 0;
        for ii = 1:num_species
            num_parts = num_parts + species_array(ii).Npoints;
        end
        totalx = zeros([1,num_parts]);
        totalweight = zeros([1,num_parts]);
        points_counter = 0;
        for ii = 1:num_species
            tempN = species_array(ii).Npoints;
            totalx(points_counter+1:points_counter+tempN) = species_array(ii).x;
            totalweight(points_counter+1:points_counter+tempN) = ...
                species_array(ii).q*species_array(ii).weights;
            points_counter = points_counter+ tempN;
        end
    end


% calc E with total weights

% break back to each species

% update each species

end

function integrate

%forward Euler:
% calculate E
% xdot = v
% break up E,v to send to individual species




% leapfrog:
% calculate E
% update v
% update x

% fe:
% don't stagger initially
% calculate E
% update v, x

% RK2
% calculate E
% calculate v*, x*
% calculate E*
% calculate v**, x**
% combine *, ** to update E,v

end

