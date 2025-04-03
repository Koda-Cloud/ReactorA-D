clc;
clear;

global xa1


solution = fsolve(@soe,[0.1, 700]);

xa1 = solution(1);

sol = fsolve(@soe2,[0.1, 372]);



function F = soe(vars)

    xa1 = vars(1);
    T = vars(2);

    V = 0.47;
    Q = 0.416e-3;

    K = 4e6*exp(-7900/T);

    eq1 = K*(1-xa1)*(V/Q) - xa1;
    eq2 = 343 + 44923.3*(1-xa1)*K - T;

    F = [eq1;eq2];
end

function F = soe2(vars)

    global xa1

    xa2 = vars(1);
    T = vars(2);

    V = 0.47;
    Q = 0.416e-3;

    K = 4e6*exp(-7900/T);

    eq1 = xa1 + K*(1-xa2)*(V/Q) - xa2;
    eq2 = 363.86 + 40926.7*(1-xa2)*K - T;

    F = [eq1;eq2];
end
