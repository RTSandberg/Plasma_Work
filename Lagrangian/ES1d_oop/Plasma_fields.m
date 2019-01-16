classdef Plasma_fields < handle
    properties
        is_periodic
        deltasq
        rhobar
        L
        field_mesh
        E_values
        sim_root
    end
    
    methods
        function obj = Plasma_fields(is_periodic, delta,rhobar,L,sim_root,Nt,field_mesh)
            obj.is_periodic = is_periodic;
            obj.deltasq = delta^2;
            obj.rhobar = rhobar;
            obj.L = L;
            obj.field_mesh = field_mesh;
            
            %initialize diagnostic
            obj.sim_root = sim_root;
            status = rmdir([ sim_root '/data/raw_field_data/'], 's');
            mkdir([ sim_root '/data/raw_field_data/raw_E'])
%             for ii = 1:Nt+1
%                 filename_E = sprintf([sim_root '/data/raw_field_data/raw_E/raw-E-field-%i.h5'],ii-1);
%                 h5create(filename_E,'/set1',[2,length(field_mesh)])
%             end
        end
        
        function E = calc_E(field_obj,targets, sources, source_charges)
            
            if field_obj.is_periodic
                sources = mod(sources,field_obj.L);
                targets = mod(targets,field_obj.L);
            end
            
            E_interact = calc_E_free_space_sortsum(field_obj,...
            targets, sources, source_charges);

            if field_obj.is_periodic
                E = E_interact + 1/field_obj.L * source_charges * sources'...
                    + field_obj.rhobar * targets;
            end
            
            field_obj.E_values = E;
            field_obj.field_mesh = targets;
        end
        
        function E_interact = calc_E_free_space(this_field,...
                targets, sources, source_charges)
             E_interact = zeros(size(targets));
                for ii = 1:length(targets)
                    E_interact(ii) = source_charges*sign(targets(ii)-...
                        sources');
                end

                E_interact = .5*E_interact;
        end
        
        function E_interact = calc_E_free_space_regular(field_obj,...
                targets, sources, source_charges)
             E_interact = zeros(size(targets));
                for ii = 1:length(targets)
                    E_interact(ii) = source_charges*((targets(ii)-sources')./...
                        sqrt((targets(ii)-sources').^2 + field_obj.deltasq));
                end

                E_interact = .5*E_interact;
        end
        
        function E_interact = calc_E_free_space_sortsum(field_obj,...
                targets, sources, charges)
            
            Ns = length(sources);
            [source_sorted, indx,ind_sorted] = unique(sources);
            
            N5 = length(targets);
            [target_sorted, indx_tar,ind_sorted_tar] = unique(targets);
            
            %[~,ind_sorted] = sort(xsvec);
            overlap = 0;
            E_interact = zeros(size(indx'));

            if length(indx)== Ns % if no points are on top of eachother
                charges_sorted = charges(indx);
%                 E_interact(1) = sum(charges_sorted(2:end));

            else
               charges_sorted = zeros(size(indx));
               for ii = 1:Ns
                   charges_sorted(ind_sorted(ii)) = ...
                       charges_sorted(ind_sorted(ii)) + charges(ii);
               end
                overlap = 1;
            %     'You have sheet overlap!'
            end
            
            % initialize E_interact
            jj = 1;
            preveq = 0;
            E_interact(1) = -sum(charges_sorted);
            while target_sorted(1) > source_sorted(jj)
                E_interact(1) = E_interact(1) + 2*charges_sorted(jj);
                jj = jj+1;
            end
                
            if abs(target_sorted(1) - source_sorted(jj)) < 1e-12
                preveq = 1;
                E_interact(1) = E_interact(1) + charges_sorted(jj);
            end
            
            %main loop
            for ii = 2:length(indx_tar)
                    E_interact(ii) = E_interact(ii-1);
                if preveq
                    E_interact(ii) = E_interact(ii) + charges_sorted(jj);
                    preveq = 0;
                    jj = jj+1;
                end
                
                while jj <= Ns && target_sorted(ii)>source_sorted(jj)
                    E_interact(ii) = E_interact(ii) + charges_sorted(jj);
                    jj = jj+1;
                end
                    
                if jj <= Ns && abs(target_sorted(ii) - source_sorted(jj)) < 1e-12
                    preveq = 1;
                    E_interact(ii) = E_interact(ii) + charges_sorted(jj);
                end
            end
            E_interact = .5*E_interact(ind_sorted);
        end
        
        
        function save_diagnostics(this_field,iteration_number)
                name1 = [this_field.sim_root '/data/raw_field_data/raw_E' ];
%                 name2 = sprintf(['/raw-' this_plasma.name '-%i.csv'],iteration_number);
                nameE = sprintf(['/raw-E-field-%i.csv'],iteration_number);
                
                csvwrite([name1 nameE],[this_field.field_mesh;this_field.E_values]);
        end
    end
end