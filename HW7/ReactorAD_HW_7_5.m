clc;
clear;

sol = fsolve(@soe,[0.01,400]);


function F = soe(vars)

    xa = vars(1);
    T = vars(2);

    K = 5e5*exp((25e3/1.987)*((1/T) - (1/323)));

    eq1 = xa^2/(1-xa)^2 - K;
    eq2 = 300 + 625*xa - T;

    F = [eq1;eq2];

end

