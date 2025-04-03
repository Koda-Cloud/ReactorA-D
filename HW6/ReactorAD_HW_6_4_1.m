clc;
clear;

sol = fsolve(@soe,[.7,298]);

T2 = 298 + 72*0.5*sol(1);

xa2 = 0.8;
xa1 = 0.325;
cao = 4;
Q = 250;

k1 = 3e7*exp(-5838/321.58);
k2 = 1.9e-11*exp(9059/321.58);

eq = @(V) xa1 + k1*((1-xa2) - (cao*(xa2^2))/k2)*(V/Q) - xa2;

sol1 = fsolve(eq,10000);

r = 18e3*k1*cao*((1-xa2) - (cao*(xa2^2))/k2)*sol1;



function F = soe(vars)

    xa = vars(1);
    T = vars(2);

    cao = 4;

    eq1 = 9059/(log(5.29e10*(xa^2/(1-xa))*cao)) - T;
    eq2 = 298 + 72*xa - T;

    F = [eq1,eq2];
end
