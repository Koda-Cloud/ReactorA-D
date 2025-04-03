clc;
clear;

P = 3;
k = 1.2e4;
Nao = 1.2;
Nto = 4.8;
a = 0.1;
R = 0.08206;
T = 650;
Da = 0.007;

Pab = @(xa) (Nao.*(1-xa).*P)./(Nto + 2.*Nao.*xa);

ca = @(xa) Pab(xa)./(R*T);

phi = @(xa) ((3.*k.*ca(xa).*a^2)/(2*Da)).^0.5;

eta = @(xa) 1./phi(xa);

r = @(xa) k.*(Pab(xa)./(R*T)).^2;

ra = @(xa) eta(xa).*r(xa);

dxa = @(xa) 1./ra(xa);

sol = Nao*integral(dxa,0,0.85);
