function accel = field2accel(field, positions, charges, masses, xmin, delx_mesh,N_mesh)

    N_part = length(positions);
    accel = zeros(N_part,1);

    for iter=1:N_part
        xi = positions(iter);
        j = 1+floor( (positions(iter)-xmin-.5*delx_mesh) / delx_mesh);

        xj = xmin + .5*delx_mesh + (j-1)*delx_mesh;
        if j== 0
            accel(iter) = accel(iter) + (xj+delx_mesh-xi)/delx_mesh * field(end);
        else
            accel(iter) = accel(iter) + (xj+delx_mesh-xi)/delx_mesh * field(j);
        end
        if j < N_mesh
            accel(iter) = accel(iter) + (xi - xj)/delx_mesh * field(j+1);
        else
            accel(iter) = accel(iter) + (xi - xj)/delx_mesh * field(1);
        end
    end

    accel = charges./masses.*accel;
% 
% for i=1:N_part
%     xi = positions(i);
%     j = 1+floor( (positions(i)-xmin-.5*delx_mesh) / delx_mesh);
%     if j == 0
%         j = N_mesh;
%     end
%     xj = xmin + .5*delx_mesh + (j-1)*delx_mesh;
%     force(i) = force(i) + (xj+delx_mesh-xi)/delx_mesh * field(j);
%     if j < N_mesh
%         force(i) = force(i) + (xi - xj)/delx_mesh * field(j+1);
%     else
%         force(1) = force(1) + (xi - xj)/delx_mesh * field(1);
%     end
% end