clc;
clear;

T = @(x) 373 + 711.1.*x;

k1 = @(x) 2e-3*exp((-100e3/8.314).*((1./T(x)) - (1/373)));
k2 = @(x) 3e-5*exp((-150e3/8.314).*((1./T(x)) - (1/373)));

eq1 = @(x) x./(((1-x).*(1.25-x)).^.5) - k1(x)./k2(x);

solution = fsolve(eq1,0.1);

T1 = .95*(373 + 711.1*solution);

xa1 = (T1 - 373)/711.1;

eq2 = @(x) 1./(k1(x).*((1-x).*(1.25-x)).^.5 - k2(x).*x);

tbatch = integral(eq2,0,xa1);
