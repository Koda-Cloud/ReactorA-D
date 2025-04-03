clc;
clear;

global sol2 km1 a P Da rhob rhop nao km2 sol3

r = 0.3;
a = 0.3/3;
P1 = 0.7;
P = 1.5;
Da = 0.007;
nao = 12;
rhob = 0.6;
rhop = 0.85;
km1 = 0.07;
km2 = 1.4;

phi = @(k) a.*(k./Da).^0.5;
eta1 = @(k) (1./phi(k)).*((1./tanh(phi(k).*3)) - 1./(3*phi(k)));

eq1 = @(k) -2.5e-5 + eta1(k).*k.*P1;

sol1 = fsolve(eq1,1);

dxa1 = @(xa) 1./((rhob/rhop).*eta1(sol1).*sol1.*(1-xa).*P);

w1 = nao*integral(dxa1,0,0.9);

bi = @(km) (km.*a)./Da;

eta2 = @(k) eta1(k)./(1 + (phi(k)./bi(km1)).*(phi(k).*eta1(k)));

eq2 = @(k) -2.5e-5 + eta2(k).*k.*((P1)  - 2.5e-5/(km1));

sol2 = fsolve(eq2,1);

[w,y] = ode45(@ex1,[0 1e6],0);

plot(w,y(:,1));

eta3 = @(k) eta1(k)./(1 + (phi(k)./bi(km2)).*(phi(k).*eta1(k)));

eq3 = @(k) -2.5e-5 + eta3(k).*k.*(P1 - 2.5e-5/km2);

sol3 = fsolve(eq3,.0001);

[w1,y1] = ode45(@ex2,[0 1e6],0);

figure
plot(w1,y1(:,1));

function F  = ex1(w,x)
    global sol2 a km1 P Da rhob rhop nao

    xa = x;

    k1 = sol2;
    Pa = (1-xa)*P;
    bi = @(km) (km.*a)./Da;
    phi = @(k) a.*(k./Da).^0.5;
    eta1 = @(k) (1./phi(k)).*((1./tanh(phi(k).*3)) - 1./(3*phi(k)));
    eta2 = eta1(k1)./(1 + (phi(k1)./bi(km1)).*(phi(k1).*eta1(k1)));

    fun = @(rw) eta2*k1*(Pa - rw./km1) - rw;

    F = (rhob/rhop)*fsolve(fun,0.0001)*(1/nao);

end

function F  = ex2(w,x)
    global sol2 a km1 P Da rhob rhop nao km2 sol3

    xa = x;

    k1 = sol3;
    Pa = (1-xa)*P;
    bi = @(km) (km.*a)./Da;
    phi = @(k) a.*(k./Da).^0.5;
    eta1 = @(k) (1./phi(k)).*((1./tanh(phi(k).*3)) - 1./(3*phi(k)));
    eta2 = eta1(k1)./(1 + (phi(k1)./bi(km2)).*(phi(k1).*eta1(k1)));

    fun = @(rw) eta2*k1*(Pa - rw./km2) - rw;

    F = (rhob/rhop)*fsolve(fun,0.0001)*(1/nao);

end
    