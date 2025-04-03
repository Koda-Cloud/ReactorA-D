clc;
clear;

global xatrue xatrue2

Q = 250;

sol = fsolve(@soe,[0.4,300]);

sol2 = fsolve(@soe2,[0.4,300]);

xatrue = (sol(1) + sol2(1))/2;

Ttrue = 300 + 72*xatrue;
cao = 4;

k = @(T) 3e7*exp(-5838./T);
K = @(T) 1.9e-11*exp(9059./T);

r = k(Ttrue)*cao*((1-xatrue) - (cao*xatrue^2)/K(Ttrue));

eq1 = @(T) r - k(T)*cao*((1-xatrue) -  (cao*xatrue^2)/K(T));

sol3 = fsolve(eq1,300);

sol4 = fsolve(@soe3,[0.4,330]);

sol5 = fsolve(@soe4,[0.4,330]);

xatrue2 = (sol4(1) + sol5(1))/2;

Ttrue2 = 327.97 + 72*(xatrue2 - xatrue);

r2 = k(Ttrue2)*cao*((1-xatrue2) - (cao*xatrue2^2)/K(Ttrue2));

eq2 = @(T) r2 - k(T)*cao*((1-xatrue2) - (cao*xatrue2^2)/K(T));

sol6 = fsolve(eq2,300);

sol7 = fsolve(@soe5,[0.4,300]);

sol8 = fsolve(@soe6,[0.4,300]);

xatrue3 = (sol7(1) + sol8(1))/2;

Ttrue3 = 321.1127 + 72*(xatrue3 - xatrue2);

T1 = @(x) 300 + 72.*x;

k1 = @(x) 3e7*exp(-5838./T1(x));
K1 = @(x) 1.9e-11*exp(9059./T1(x));

dxa1 = @(x) 1./(k1(x).*(cao.*(1-x) - (cao.*x).^2/K1(x)));

V1 = Q*cao*integral(dxa1,0,xatrue);

T2 = @(x) 327.97 + 72.*(x - xatrue);

k2 = @(x) 3e7*exp(-5838./T2(x));
K2 = @(x) 1.9e-11*exp(9059./T2(x));

dxa2 = @(x) 1./(k2(x).*(cao.*(1-x) - (cao.*x).^2/K2(x)));

V2 = Q*cao*integral(dxa2,xatrue,xatrue2);

T3 = @(x) 321.1127 + 72.*(x - xatrue2);

k3 = @(x) 3e7*exp(-5838./T3(x));
K3 = @(x) 1.9e-11*exp(9059./T3(x));

dxa3 = @(x) 1./(k3(x).*(cao.*(1-x) - (cao.*x).^2/K3(x)));

V3 = Q*cao*integral(dxa3,xatrue2,xatrue3);

function F = soe(vars)
    
    xa = vars(1);
    T = vars(2);

    cao = 4;

    eq1 = 9059/(log(5.28e10*(xa^2/(1-xa))*cao)) - T;
    eq2 = 300 + 72*xa -T;

    F = [eq1;eq2];
end

function F = soe2(vars)
    
    xa = vars(1);
    T = vars(2);

    cao = 4;

    eq1 = 9059/(log(1.34e11*(xa^2/(1-xa))*cao)) - T;
    eq2 = 300 + 72*xa -T;

    F = [eq1;eq2];
end

function F = soe3(vars)
    global xatrue
    xa = vars(1);
    T = vars(2);

    cao = 4;

    eq1 = 9059/(log(5.28e10*(xa^2/(1-xa))*cao)) - T;
    eq2 = 327.97 + 72*(xa - xatrue) -T;

    F = [eq1;eq2];
end

function F = soe4(vars)
    global xatrue
    xa = vars(1);
    T = vars(2);

    cao = 4;

    eq1 = 9059/(log(1.34e11*(xa^2/(1-xa))*cao)) - T;
    eq2 = 327.97 + 72*(xa - xatrue) -T;

    F = [eq1;eq2];
end

function F = soe5(vars)
    global xatrue2
    xa = vars(1);
    T = vars(2);

    cao = 4;

    eq1 = 9059/(log(5.28e10*(xa^2/(1-xa))*cao)) - T;
    eq2 = 321.1127 + 72*(xa - xatrue2) -T;

    F = [eq1;eq2];
end

function F = soe6(vars)
    global xatrue2
    xa = vars(1);
    T = vars(2);

    cao = 4;

    eq1 = 9059/(log(1.34e11*(xa^2/(1-xa))*cao)) - T;
    eq2 = 321.1127 + 72*(xa - xatrue2) -T;

    F = [eq1;eq2];
end
