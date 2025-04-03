clc;
clear;

y = true;

eps = 0.412;
r = 0.3;
a = r/3;
Da = 0.004;
P = 2;
R = 0.08206;
nao = 12;

Pa = @(xa) (1-xa)./(1+xa).*P;

k1 = @(T) 2.35e13*exp(-19124./T);
k2 = @(T) 2.35e12*exp(-19124./T);

phi = @(T) a*(k1(T)./Da)^0.5;
eta = @(T) (1./phi(T))*((1./tanh(phi(T).*3)) - 1./(3.*phi(T)));

ca = @(xa,T) Pa(xa)./(R*T);

dxa1 = @(xa) 1./((1-eps).*k1(578).*eta(578).*ca(xa,578));

sol = nao*integral(dxa1,0,0.9);

phi2 = @(T) a.*(k2(T)./Da)^0.5;
eta2 = @(T) (1./phi2(T)).*((1./tanh(phi2(T).*3)) - 1./(3.*phi2(T)));

Ti = 578;

while y == true

dxa2 = @(xa) 1./((1-eps).*k2(Ti).*eta2(Ti).*(ca(xa,Ti)));

eq1 = nao*integral(dxa2,0,0.9);

sol1 = sol - eq1;

if abs(sol1) < 1
    y = false;
else
    y = true;
    T = Ti;
    Ti = T + 1;
end

end
