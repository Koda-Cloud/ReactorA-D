clc;
clear;

sol = fsolve(@soe,[0.001,1]);

function F = soe(vars)

    ca1 = vars(1);
    tau = vars(2);

    cao = 1.2;
    ca2 = 0.12;
    k = 5.5;

    phi1 = 14.61*ca1^0.5;
    phi2 = 14.61*ca2^0.5;

    eta1 = (1/phi1)*((1/tanh(3*phi1)) - 1/(3*phi1));
    eta2 = (1/phi2)*((1/tanh(3*phi2)) - 1/(3*phi2));

    eq1 = cao - ca1 - k*eta1*ca1^2*tau;
    eq2 = ca1 - ca2 - k*eta2*ca2^2*tau;

    F = [eq1;eq2];

end
