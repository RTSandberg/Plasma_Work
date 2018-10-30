function density = weight(positions, velocities, charges, N_mesh, xmin, delx_mesh)

%% simple test:
% density = weight(3/8,0,1,4,0,1/4) 
% gives [0,1,0,0]'
% density = weight(3/8+.01,0,1,4,0,1/4)
% give [0 0.9600 .0400 0]'
% density = weight(.4,0,1,1,10,0,1)
% gives [.9 0 0 0 0 0 0 0 0 .1]'
density = zeros(N_mesh,1);
current_density = zeros(N_mesh,1);
N_part = size(positions);
order = 1;
if order == 1
    for i=1:N_part
        xi = positions(i);
%        vi = velocities(i);
        j = 1 + floor( (positions(i)-xmin-delx_mesh*.5) / delx_mesh);

            xj = xmin + .5*delx_mesh +(j-1)*delx_mesh;
        if j == 0

            density(N_mesh) = density(N_mesh) + (xj+delx_mesh-xi)/delx_mesh^2*charges(i);
            density(j+1) = density(j+1) + (xi - xj)/delx_mesh^2*charges(i);
%            current_density(N_mesh) = current_density(N_mesh) + vi*(xj+delx_mesh-xi)/delx_mesh*W(i);
%            current_density(j+1) = current_density(j+1) + vi*(xj+delx_mesh-xi)/delx_mesh*W(i);
        else
            density(j) = density(j) + (xj + delx_mesh - xi)/delx_mesh^2*charges(i);
%            current_density(j) = current_density(j) + vi*(xj+delx_mesh-xi)/delx_mesh*W(i);
            if j < N_mesh
                density(j+1) = density(j+1) + (xi - xj)/delx_mesh^2*charges(i);
%                current_density(j) = current_density(j) + vi*(xj+delx_mesh-xi)/delx_mesh*W(i);
            else % j == N_mesh            
                density(1) = density(1) + (xi - xj)/delx_mesh^2*charges(i);
%                current_density(1) = current_density(1) + vi*(xj+delx_mesh-xi)/delx_mesh*W(i);
            end
        end
    end

else
    %2nd order weighting
    delx = delx_mesh;
    weights = charges/delx;
delxsq = delx_mesh^2;
for ii = 1:N_part
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

        if jj <N_mesh

        density(jj-1) = density(jj-1) + .5*(1.5 - (xii-xjj+delx)/delx)^2*weights(ii);
        density(jj) = density(jj) + (.75 - (xii-xjj)^2/delxsq)*weights(ii); 
        density(jj+1) = density(jj+1) + .5*(1.5 - (xjj+delx-xii)/delx)^2*weights(ii);

%             density(jj+1) = density(jj+1) + (xii-xjj)/delx*delv*weights(ii);
        elseif jj == N_mesh

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
end
