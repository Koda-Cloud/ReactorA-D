    clc;
clear;

global xatrue1 xatrue2

Q = 2.5e-4;
cao = 2e3;
cbo = 3e3;

kf = @(T) exp(7.16 - 10100./T);
kr = @(T) exp(30.4 - 16100./T);

sol = fsolve(@soe,[0.01,300]);

sol1 = fsolve(@soe2,[0.01,500]);

xatrue1 = (sol(1) + sol1(1))/2;

Ttrue1 = 300 + 640*xatrue1;

r1 = kf(Ttrue1).*cao.*(1-xatrue1).*(cbo - cao*xatrue1) - kr(Ttrue1).*cao*xatrue1;

eq1 = @(T) r1 - kf(T).*cao.*(1-xatrue1).*(cbo - cao*xatrue1) + kr(T).*cao*xatrue1;

Tf2 = fsolve(eq1,350);

sol2 = fsolve(@soe3,[0.01,300]);

sol3 = fsolve(@soe4,[0.01,500]);

xatrue2 = (sol2(1) + sol3(1))/2;

Ttrue2 = 406.5 + 640*(xatrue2 - xatrue1);

r2 = kf(Ttrue2).*cao.*(1-xatrue2).*(cbo - cao*xatrue2) - kr(Ttrue2).*cao*xatrue2;

eq2 = @(T) r2 - kf(T).*cao.*(1-xatrue2).*(cbo - cao*xatrue2) + kr(T).*cao*xatrue2;

Tf3 = fsolve(eq2,400);

sol4 = fsolve(@soe5,[0.01,300]);

sol5 = fsolve(@soe6,[0.01,500]);

xatrue3 = (sol4(1) + sol5(1))/2;

Ttrue3 = 402.54 + 640*(xatrue3 - xatrue2);

Q1 = -Q*((cao*(1-xatrue1) + (cbo - cao*xatrue1))*25*-(Tf2 - Ttrue1) + cao*xatrue1*50*-(Tf2 - Ttrue1));

Q2 = -Q*((cao*(1-xatrue2) + (cbo - cao*xatrue2))*25*-(Tf3 - Ttrue2) + cao*xatrue2*50*-(Tf3 - Ttrue2));

T1 = @(x) 300 + 640*x;

kf1 = @(x) exp(7.16 - 10100./T1(x));
kr1 = @(x) exp(30.4 - 16100./T1(x));

dxa1 = @(x) 1./(kf1(x).*cao.*(1-x).*(cbo - cao*x) - kr1(x).*cao.*x);

V1 = Q*cao*integral(dxa1,0,xatrue1);

T2 = @(x) 406.5 + 640*(x - xatrue1);

kf2 = @(x) exp(7.16 - 10100./T2(x));
kr2 = @(x) exp(30.4 - 16100./T2(x));

dxa2 = @(x) 1./(kf2(x).*cao.*(1-x).*(cbo - cao*x) - kr2(x).*cao.*x);

V2 = Q*cao*integral(dxa2,xatrue1,xatrue2);

T3 = @(x) 402.54 + 640*(x - xatrue2);

kf3 = @(x) exp(7.16 - 10100./T3(x));
kr3 = @(x) exp(30.4 - 16100./T3(x));

dxa3 = @(x) 1./(kf3(x).*cao.*(1-x).*(cbo - cao*x) - kr3(x).*cao.*x);

V3 = Q*cao*integral(dxa3,xatrue2,xatrue3);

function F = soe(vars)

    xa = vars(1);
    T = vars(2);

    cao = 2e3;
    cbo = 3e3;

    kf = exp(7.16 - 10100/T);
    kr = exp(30.4 - 16100/T);

    eq1 = xa/((1-xa)*(cbo - cao*xa)) - kf/kr;
    eq2 = 300 + 640*xa - T;

    F = [eq1;eq2];

end

function F = soe2(vars)

    xa = vars(1);
    T = vars(2);

    cao = 2e3;
    cbo = 3e3;

    kf = exp(7.16 - 10100/T);
    kr = exp(30.4 - 16100/T);

    eq1 = 10100*kf*(1-xa)*(cbo - cao*xa) - 16100*kr*xa;
    eq2 = 300 + 640*xa - T;

    F = [eq1;eq2];

end

function F = soe3(vars)
global xatrue1
    xa = vars(1);
    T = vars(2);

    cao = 2e3;
    cbo = 3e3;

    kf = exp(7.16 - 10100/T);
    kr = exp(30.4 - 16100/T);

    eq1 = xa/((1-xa)*(cbo - cao*xa)) - kf/kr;
    eq2 = 406.5 + 640*(xa - xatrue1) - T;

    F = [eq1;eq2];

end

function F = soe4(vars)
global xatrue1
    xa = vars(1);
    T = vars(2);

    cao = 2e3;
    cbo = 3e3;

    kf = exp(7.16 - 10100/T);
    kr = exp(30.4 - 16100/T);

    eq1 = 10100*kf*(1-xa)*(cbo - cao*xa) - 16100*kr*xa;
    eq2 = 406.5 + 640*(xa - xatrue1) - T;

    F = [eq1;eq2];

end

function F = soe5(vars)
global xatrue2
    xa = vars(1);
    T = vars(2);

    cao = 2e3;
    cbo = 3e3;

    kf = exp(7.16 - 10100/T);
    kr = exp(30.4 - 16100/T);

    eq1 = xa/((1-xa)*(cbo - cao*xa)) - kf/kr;
    eq2 = 402.54 + 640*(xa - xatrue2) - T;

    F = [eq1;eq2];

end

function F = soe6(vars)
global xatrue2
    xa = vars(1);
    T = vars(2);

    cao = 2e3;
    cbo = 3e3;

    kf = exp(7.16 - 10100/T);
    kr = exp(30.4 - 16100/T);

    eq1 = 10100*kf*(1-xa)*(cbo - cao*xa) - 16100*kr*xa;
    eq2 = 402.54 + 640*(xa - xatrue2) - T;

    F = [eq1;eq2];

end
