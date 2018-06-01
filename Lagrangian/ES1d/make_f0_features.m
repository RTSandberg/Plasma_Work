%%% make_f0vec
%%% generate f0vec by features
%%% makes the distribution vector by 
%%% mesh
%%% input: 
%%%          input_deck:    .mat file of parameters for defining initial grid
%%%          beams:         struct array of beams and beam parameters
%%%          spots:         struct array of spots and spot parameters
%%%                     IMPORTANT NOTE: if either beams or spots is empty,
%%%                     must be given as struct([])
%%%
%%%          spot fields:   x:      scalar, x coordinate of spot,
%%%                         v:      scalar, v coordinate of spot
%%%                         N:      scalar, integral of spot = number of
%%%                         particles represented by spot
%%%          beam fields:   v:      scalar, v value of beam
%%%                         vth:    scalar, thermal velocity or spread of beam
%%%                         perturb:character, indicate whether perturbation is desired
%%%                                     'n' : no perturbation desired
%%%                                     's' : sinusoidal perturbation for
%%%                                     standing wave
%%%                                     't' : traveling sine wave, wave
%%%                                     speed determined by wavelength
%%%                                     'x' : pulse perturbation in x
%%%                         amplitude:scalar, gives amplitude of
%%%                         perturbation
%%%                         wavelength:scalar, gives wavelength of sinusoid
%%%                         perturbation.  Currently does nothing for pulse
%%%                         perturbation
%%%                         locationphase: scalar, means nothing for
%%%                         sinusoidal perturbation, gives location of
%%%                         pulse perturbation
%%%                         n0:     scalar, base density of beam
%%%
%%%         
%%% output:
%%%         f0vec:  vector of initial f values, evaluated at points
%%%         determined by x, v mesh

%%%

function f0 = make_f0_features(input_deck,spots,beams,shear)


    load(input_deck)

    delx = L/Nx;
%     vmin = -vmax;
        
    x0 = xmin+.5*delx : delx : xmin + L;
    v0 = vmin + .5*delv : delv : vmax; Nv = length(v0);
    [X0,V0] = meshgrid(x0, v0);
    if shear.shearQ
    % shearing the inital distribution
        X0 = X0 + shear.slope*V0;
        X0 = mod(X0 - xmin,L)+xmin;
    end
    
    f0 = zeros(size(X0));
    
    % currently, only supports cold spots.  Don't know what a warm spot
    % means
    

    for spotindex = 1:numel(spots)
        spot = spots(spotindex);
        
        x1 = spot.x; v1 = spot.v;
        x1 = mod(x1-xmin,L)+xmin;

        jj = ceil((x1-xmin)/delx);
        ii = ceil((v1-vmin)/delv);
        if ii > Nv; ii = Nv; end
        if ii < 1; ii = 1; end
        f0(ii,jj) = 1/delx/delv * spot.N;
    end
    
    for beamindex = 1:numel(beams)
        
        beam = beams(beamindex);
           
        ii = round((beam.v-vmin)/delv);
        if ii > Nv; ii = Nv; end
        if ii < 1; ii = 1; end
        
        n0 = beam.n0;
        
        k = 2*pi/beam.wavelength;
        vth = beam.vth;

        if vth <= 0 && vth > -1
            f0(ii,:) = n0 * 1/delv*ones([1,Nx]);
            
           if beam.perturb == 's' %standing sine       
               
                
                f0(ii,:) = f0(ii,:) +  1./delx/delv * beam.amplitude*sin(k * x0) ;
               
           elseif beam.perturb == 't' %traveling sine
               vp = 1./k;
               amp = beam.amplitude;
               f0 = (n0/delv + 1./delx/delv * amp*sin(k * X0) )...
                    .* ( ( (V0 - vp*amp/delx/n0*sin(k*X0) ) >= -.5*delv)...
                    & ( (V0 - vp*amp/delx/n0*sin(k*X0) ) < .5*delv) );
                
               
           elseif beam.perturb == 'x' %pulse in x
               % currently only support a single cell perturbation
               x1 = mod(beam.locationphase-xmin,L) + xmin;
               jj = ceil((x1-xmin)/delx);
               if jj <= 0; jj = 1; end;
               amp = beam.amplitude;
               
               f0(ii,:) =  ones([1,Nx]);
               f0(ii,jj) = f0(ii,jj) + beam.amplitude;
               f0(ii,:) = n0 * 1/delv/delx/(Nx+amp) * f0(ii,:);
               
           elseif beam.perturb == 'v' %pulse in v 
           end
           
        elseif beam.vth < -1 % uniform beam
            iithick = round(beam.amplitude  / delv);
            
%             f0(ii-iithick:ii+iithick,:) = n0 * 1/delv*ones([2*iithick + 1,Nx]);
            f0 = n0/2/vmax * ones(size(V0));
           if beam.perturb == 's' %standing sine       
                f0 = f0 +  1./delx/delv * beam.amplitude*sin(k * X0) ;
           end
           
        else % warm beams with Gaussian spread
           f0 = n0 * 1/sqrt(pi)/vth*exp(-V0.^2/vth^2) ;
           
           if beam.perturb == 's'
               f0 = 1/sqrt(pi)/vth.*exp(-(V0-beam.v).^2/vth^2)...
                   .* (n0+ beam.amplitude.*sin(k * X0));
           elseif beam.perturb == 't'
               vp = 1./k;
               f0 =  1/sqrt(pi)/vth...
                   *exp(-(V0-beam.v - vp*beam.amplitude/n0/delx*sin(k*X0)).^2/vth^2)...
                   .* (n0 + beam.amplitude/delx.*sin(k * X0));
           elseif beam.perturb == 'x'               
               x1 = mod(beam.locationphase-xmin,L) + xmin;
               jj = ceil((x1-xmin)/delx);
               if jj <= 1; jj = 1; end
               f0(:,jj) = (f0(:,jj) + beam.amplitude);
               f0 = f0/ (1+Nv*delx*delv*beam.amplitude);
           elseif beam.perturb == 'v'
               
           end
               
        end
    end
    
%    sum(sum(f0))*delx*delv % checking normalization
end

function set_grid_params()
%%% from spots and beams struct arrays, determines vmax, delv, delx, Nx
%%% set x0, v0, X0, V0

% for spots
spotxlist = [];
spotvlist = [];
for spot = spots
    spotxlist = [spotxlist spot.x];
    spotvlist = [spotvlist spot.v];
end
xmat = spotxlist'*ones(size(spotxlist));
xdiffs = abs(xmat - xmat');
delx = min(min(xdiffs));
vmat = spotvlist'*ones(size(spotvlist));
vdiffs = abs(vmat - vmat');
delv = min(min(vdiffs));

end
