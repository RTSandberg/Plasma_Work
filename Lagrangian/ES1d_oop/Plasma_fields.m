classdef Plasma_fields
    properties
        is_periodic
        deltasq
        rhobar
        L
        field_mesh
        E_values
    end
    
    methods
        function obj = Plasma_fields(is_periodic, delta,rhobar,L,sim_root,Nt)
            obj.is_periodic = is_periodic;
            obj.deltasq = delta^2;
            obj.rhobar = rhobar;
            obj.L = L;
        end
        
        function E = calc_E(field_obj,targets, positions, weights)
            
            [E_interact, neutralizing] = calc_E_free_space(field_obj,...
            targets, positions, weights);

            if field_obj.is_periodic
                E = E_interact + neutralizing + field_obj.rhobar * targets;
            end
            
            field_obj.E_values = E;
            field_obj.field_mesh = targets;
        end
        
        function [E_interact, neutralizing] = calc_E_free_space(field_obj,...
                targets, positions, weights)
             E_interact = zeros(size(targets));
                for ii = 1:length(targets)
                    E_interact(ii) = weights*((targets(ii)-positions')./...
                        sqrt((targets(ii)-positions').^2 + field_obj.deltasq));
                end

                E_interact = .5*E_interact;
                neutralizing = 1/field_obj.L * weights * positions';
                
        end
        
    end
end