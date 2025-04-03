clear;
clc;

T = @(x) 624835./(949.5 - 16.5.*x) - 308.061;

K = @(x) 1e3*exp(-5033./T(x));

dx = @(x) 1./(K(x).*(1-x));

tbatch = integral(dx,0,.75);

eq1 = @(x) (3600 + integral(dx,0,x)).*(K(x).*(1-x)) - x;

sol = fsolve(eq1, .75);

tbatch1 = integral(dx,0,sol);
