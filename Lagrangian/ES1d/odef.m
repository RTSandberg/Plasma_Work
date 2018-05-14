function atotvec = odef(x,f0vec,c1,c2,L)
N = length(f0vec);
xvec = x(1:N); vvec = x(N+1:end);
alpha = c1/L * xvec'*f0vec;

[xj,xi] = meshgrid(xvec,xvec);
Kmat = .5*sign(xi-xj);
avec = c1*Kmat*f0vec + alpha + c2 * xvec;

avec = avec;

atotvec = [vvec;avec];