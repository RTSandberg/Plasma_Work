% density diagnostic
function density = xweight(positions, weights, meshpositions, delx,xmin,delv)
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
weights = weights * delv;
order = 2;
%1st order weighting
if order == 1
for ii = 1:Np
    xii = positions(ii);
    jj = floor((xii - xmin-.5*delx)/ delx ) + 1;
    xjj = xmin + .5*delx +(jj-1)*delx;
    if jj < 1
        density(end) = density(end) + (xjj+delx-xii)/delx*weights(ii); 
        density(jj+1) = density(jj+1) + (xii-xjj)/delx*weights(ii);
    else
        density(jj) = density(jj) + (xjj+delx-xii)/delx*weights(ii);
        if jj <Nm
            density(jj+1) = density(jj+1) + (xii-xjj)/delx*weights(ii);
        else
            density(1) = density(1) + (xii-xjj)/delx*weights(ii);
        end
    end
end
elseif order == 2        
        

%2nd order weighting
delxsq = delx^2;
for ii = 1:Np
    xii = positions(ii);
    jj = floor((xii - xmin)/ delx ) + 1;
    xjj = xmin + .5*delx +(jj-1)*delx;
    
    if jj < 1
        density(end-1) = density(end-1) + .5*(1.5 - (xii-xjj+delx)/delx)^2*weights(ii);
        density(end) = density(end) + (.75 - (xii-xjj)^2/delxsq)*weights(ii); 
        density(jj+1) = density(jj+1) + .5*(1.5 - (xjj+delx-xii)/delx)^2*weights(ii);
    elseif jj == 1
        density(end) = density(end) + .5*(1.5 - (xii-xjj+delx)/delx)^2*weights(ii);
        density(1) = density(1) + (.75 - (xii-xjj)^2/delxsq)*weights(ii); 
        density(jj+1) = density(jj+1) + .5*(1.5 - (xjj+delx-xii)/delx)^2*weights(ii);
        
%         density(jj) = density(jj) + (xjj+delx-xii)/delx*delv*weights(ii);

    else

        if jj <Nm

        density(jj-1) = density(jj-1) + .5*(1.5 - (xii-xjj+delx)/delx)^2*weights(ii);
        density(jj) = density(jj) + (.75 - (xii-xjj)^2/delxsq)*weights(ii); 
        density(jj+1) = density(jj+1) + .5*(1.5 - (xjj+delx-xii)/delx)^2*weights(ii);

%             density(jj+1) = density(jj+1) + (xii-xjj)/delx*delv*weights(ii);
        elseif jj == Nm

        density(jj-1) = density(jj-1) + .5*(1.5 - (xii-xjj+delx)/delx)^2*weights(ii);
        density(jj) = density(jj) + (.75 - (xii-xjj)^2/delxsq)*weights(ii); 
        density(1) = density(1) + .5*(1.5 - (xjj+delx-xii)/delx)^2*weights(ii);

%             density(1) = density(1) + (xii-xjj)/delx*delv*weights(ii);
        else

        density(jj-1) = density(jj-1) + .5*(1.5 - (xii-xjj+delx)/delx)^2*weights(ii);
        density(1) = density(1) + (.75 - (xii-xjj)^2/delxsq)*weights(ii); 
        density(2) = density(2) + .5*(1.5 - (xjj+delx-xii)/delx)^2*weights(ii);
        end
    end
end    
    
end