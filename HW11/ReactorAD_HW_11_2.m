clc;
clear;

sol = fsolve(@soe,[1,6]);

function F = soe(vars)

    r = vars(1);
    k = vars(2);

    ca = 0.24;
    cao = 1.2;

    cas = ca - r/12.5;


    phi = 4.382*(k*(cas^0.5))^0.5;
    eta = (1.0357 + 0.3173*phi + 0.000437*phi^2)/(1 + 0.4172*phi + 0.139*phi^2);
    
    eq1 = r - eta*k*cas^1.5;
    eq2 = cao - ca - eta*k*cas^1.5*5;

    F = [eq1;eq2];

end
