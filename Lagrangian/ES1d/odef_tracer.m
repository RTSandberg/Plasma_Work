function atotvec = odef_tracer(x,xtracer,f0vec,c1,c2,L)
N = length(f0vec); 



Nt = floor(length(xtracer)/2);
xvec = x(1:N); vvec = x(N+1:end);
xtvec = xtracer(1:Nt); vtvec = xtracer(Nt+1:end);
alpha = c1/L * xvec'*f0vec;

[xj,xi] = meshgrid(xvec,xtvec);
Kmat = .5*sign(xi-xj);
avec = c1*Kmat*f0vec + alpha + c2 * xtvec;

atotvec = [vtvec;avec];