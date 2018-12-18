% the field class
classdef LPM_species < species
    properties
        deltasq
    end
    methods
        function obj = LPM_species(name,q,m,type,do_amr,Npoints,L,delta,sim_root,Nt)
            obj = obj@species(name,q,m,type,do_amr,Npoints,L,sim_root,Nt);
            obj.deltasq = delta^2;
        end
        function [E_interact,E_neutral] = calc_E(this_species,targets,L)
            
            E_interact = zeros(size(targets));
            for ii = 1:length(targets)

                E_interact(ii) = this_species.weights* ...
                    ((targets(ii)-this_species.x')./...
                    sqrt((targets(ii)-this_species.x').^2 + this_species.deltasq));
            end
                
                E_interact = .5*this_species.q*E_interact;
                E_neutral = this_species.q/L * this_species.weights ...
                    * this_species.x';
        end
        
        function [E_interact,E_neutral] = calc_E_fast(this_species,targets,L)
            %temp = sort(xsvec);
            [~, indx,ind_sorted] = unique(this_species.x);
            indx = indx';
            ind_sorted = ind_sorted'; % convert MATLAB default columns from unique to rows
            %[~,ind_sorted] = sort(xsvec);
            overlap = 0;
            E_interact = zeros(size(indx));

            if length(indx)==this_species.Npoints % if no points are on top of eachother
                f0_sorted = this_species.weights(indx);
                E_interact(1) = sum(f0_sorted(2:end));

            else
               f0_sorted = zeros(size(indx));
               for ii = 1:this_species.Npoints
                   f0_sorted(ind_sorted(ii)) = f0_sorted(ind_sorted(ii)) ...
                       + this_species.weights(ii);
               end
                overlap = 1;
            %     'You have sheet overlap!'
            end

            E_interact(1) = -sum(f0_sorted(2:end));
            for ii = 2:length(indx)
                E_interact(ii) = E_interact(ii-1) + f0_sorted(ii)+f0_sorted(ii-1);
            end
            
                E_interact = .5*this_species.q*E_interact;
                E_neutral = this_species.q/L * this_species.weights ...
                    * this_species.x';
        end
        
        function phi = calc_potential(this_species,targets,L)
                phi = zeros(size(this_species.x));
        
                for ii = 1:this_species.Npoints
                    tempvec = (targets(ii)-this_species.x).^2;
                    phi(ii) = (sqrt(tempvec+this_species.deltasq)-...
                        tempvec/L)*this_species.weights';
                end
                phi = -.5*this_species.q*phi;
        end
        function potential = calc_potential_fast(obj,targets)
        end
    end
    
end