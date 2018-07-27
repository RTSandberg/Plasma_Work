
function coeffs = linear_growth(x,y,lineara,linearb)
ind = find(x > lineara);
xlinear = x(ind);
y = y(ind);

ind = find(xlinear < linearb);
xlinear = xlinear(ind);
y = y(ind);

X = [ones(size(xlinear)),xlinear];
coeffs = X\y;

figure
plot(xlinear,y,xlinear,coeffs(1) + coeffs(2)*xlinear)
