clc;
clear;

tau = integral(@(t) t,1,2);

E = @(t) (0).*(t < 1) + (1).*( t >= 1 & t <= 2) + (0).*(t > 2);

[t,x] = ode45(@soe,[0 1e5],0);

xa = trapz(t,x.*E(t));

[t1,x1] = ode45(@soe2,[0 1e5],0);

xa1 = trapz(t,x1.*E(t));

function F = soe(t,x)

    k = 0.3;
    cao = 2;

    dxa = k*cao^2*(1-x)^3;

    F = dxa;
end

function F = soe2(t,x)

    cao = 2;
    
    T = 305 + 160.*x;
    k = 0.3*exp(-2000*((1/T) - 1/300));
    
    dxa = k*cao^2*(1-x)^3;

    F = dxa;
end
