%%% field_solve
function [E,phibar,dbar]= field_solve(density, delx_mesh)
 
    Nm = length(density);
    dbar = fft(density);
    kvec = (0:Nm-1)';
    Kvec = 1./delx_mesh*sin(2 *pi* kvec/Nm);
    kappavec = 2/delx_mesh * sin(pi*kvec/Nm);
    phibar = zeros(Nm,1);
    phibar(2:end) = dbar(2:end)./kappavec(2:end).^2;
    Ebar = -1i * Kvec.*phibar;
    E = ifft(Ebar);
    E = real(E);