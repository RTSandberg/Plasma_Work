% density diagnostic
function density = xweight(positions, f0vec, meshpositions, delx,xmin,delv,q)
% take in list of positions and densities
% list of mesh positions and parameters
% uses 0th order weighting to calculate density
% for viewing purposes
% loop through phase space particles
% for each, find the cell containing it and add density

Np = length(positions);
Nm = length(meshpositions);
density = zeros([Nm,1]);
% for ii = 1:Np
%     xii = positions(ii);
%     jj = floor((xii - xmin)/ delx ) + 1;
%         % xii = xmin + delx*(jj-1) + r, 0<= r < 1
%     density(jj) = density(jj) + delv*q*f0vec(ii);
% end       

%1st order weighting
for ii = 1:Np
    xii = positions(ii);
    jj = floor((xii - xmin-.5*delx)/ delx ) + 1;
    xjj = xmin + .5*delx +(jj-1)*delx;
    if jj < 1
        density(end) = density(end) + (xjj+delx-xii)/delx*delv*q*f0vec(ii); 
        density(jj+1) = density(jj+1) + (xii-xjj)/delx*delv*q*f0vec(ii);
    else
        density(jj) = density(jj) + (xjj+delx-xii)/delx*delv*q*f0vec(ii);
        if jj <Nm
            density(jj+1) = density(jj+1) + (xii-xjj)/delx*delv*q*f0vec(ii);
        else
            density(1) = density(1) + (xii-xjj)/delx*delv*q*f0vec(ii);
        end
    end
end
        
        
        