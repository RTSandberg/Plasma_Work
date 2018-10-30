% density diagnostic
function distribution = xvweight(positions, velocities, f0vec, meshpositions, delx,xmin, delv, vmax)
% take in list of positions and densities
% list of mesh positions and parameters
% uses 0th order weighting to calculate density
% for viewing purposes
% loop through phase space particles
% for each, find the cell containing it and add density

Np = length(positions);
Nx = length(meshpositions);
Nv = floor( (2*vmax) / delv );
distribution = zeros([Nv,Nx]);
vmin = -vmax;
for ii = 1:Np
    xii = positions(ii);
        % xii = xmin + delx*(jj-1) + r, 0<= r < 1
    vii = velocities(ii);
    jj = floor((xii - xmin)/ delx ) + 1;
    kk = floor((vii - vmin)/ delv ) +1;
    if kk < 1
        kk = 1;
    end
    if kk > Nv
        kk = Nv;
    end
    distribution(kk,jj) = distribution(kk,jj) + f0vec(ii);
end  