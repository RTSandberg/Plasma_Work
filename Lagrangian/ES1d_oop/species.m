classdef species < handle
    % Plasma species object
    properties
        name            % name of species
        q               % species charge
        m               % species mass
        type            % is species a stream, a patch, a point, or distributed
        do_amr          % IF species is stream, do point insertion
        
        Npoints
        alphas
        x
        v
        f0
        weights
        
        integration_method
        
        sim_root
        % list of diagnostics to perform
    end
    methods
        function obj = species(name,q,m,type,do_amr,Npoints,L,sim_root,Nt)
            obj.name = name;
            obj.q = q;
            obj.m = m;
            
            obj.Npoints = Npoints;
            delalpha = L/Npoints;
            obj.alphas = 0:delalpha:L-.5*delalpha;
            obj.x = obj.alphas+ .01*cos(2*pi/L*obj.alphas);
            obj.v = 0*ones(size(obj.alphas));
            obj.f0 = ones(size(obj.alphas));
            obj.weights = obj.f0 * delalpha;
            
            %set up diagnostics
            obj.sim_root = sim_root;
            status = rmdir([ sim_root '/data/raw_particle_data/' name ], 's');
            mkdir([ sim_root '/data/raw_particle_data/' name ])
            
            for ii = 1:Nt+1
                filename = sprintf([sim_root '/data/raw_particle_data/' name '/raw-' name '-%i.h5'],ii-1);
                filename_E = sprintf([sim_root '/data/raw_particle_data/' name '/raw-' name '-fields-%i.h5'],ii-1);
                h5create(filename,'/set1',[2,Npoints])
                h5create(filename_E,'/set1',[2,Npoints])
            end
        end
        
%         function 
        
        function update_v(this_plasma, vdot, delt)
            this_plasma.v = this_plasma.v + delt*vdot;
        end
        
        function update_x(this_plasma, xdot, delt)
            this_plasma.x = this_plasma.x + delt*xdot;
        end
        
        function save_diagnostics(this_plasma,iteration_number)
            
%                 filename = ;
%                 h5create(filename,'/set1',[2,Npoints])
                name1 = [this_plasma.sim_root '/data/raw_particle_data/' this_plasma.name];
                name2 = sprintf(['/raw-' this_plasma.name '-%i.h5'],iteration_number);
                nameE = sprintf(['/raw-' this_plasma.name '-fields-%i.h5'],iteration_number);
                
                h5write([name1 name2],'/set1',[this_plasma.x;this_plasma.v]);
                h5write([name1 nameE],'/set1',[this_plasma.x;this_plasma.v]);
                
%                 data = h5read(filename,'/set1');
% xtemp = data(1,:); 
% vtemp = data(2,:);
        end
    end
    methods(Abstract=true)
        calc_E
        calc_potential
    end
end