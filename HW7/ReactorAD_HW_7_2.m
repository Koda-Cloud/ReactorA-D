clc;
clear;

global xa1pre xa2pre

nao = 62;

sol = fsolve(@soe,[0.01,900]);

Ttrue1 = sol(2) - 25;

xa1 = (Ttrue1 - 723)/479.1;

Tf2 = (0.4*Ttrue1 + 0.3*648)/0.7;

na1f = 0.4*nao*(1-xa1) + 0.3*nao;


xa1pre = (0.7*nao - na1f)/(0.7*nao);

sol2 = fsolve(@soe2,[0.01,900]);

Ttrue2 = sol2(2) - 25;

xa2 = (Ttrue2 - 809.32)/479.1 + xa1pre;

Tf3 = (0.7*Ttrue2 + 0.3*648);

na2f = 0.7*nao*(1-xa2) + 0.3*nao;

xa2pre = (nao - na2f)/nao;

sol3 = fsolve(@soe3,[0.01,900]);

Ttrue3 = sol3(2) - 25;

xa3 = (Ttrue3 - 836.56)/479.1 + xa2pre;

na3 = nao*(1-xa3);

sol4 = fsolve(@soe4,[0.01,900]);

xa3loc = (nao - na3)/nao;


function F = soe(vars)

    xa = vars(1);
    T = vars(2);

    eq1 = 11835/log(((7.825e4*xa)/(1-xa))*((1-0.0475*xa)/(0.115 - 0.0475*xa))^0.5) - T;
    eq2 = 723 + 479.1*xa - T;

    F = [eq1;eq2];

end

function F = soe2(vars)
global xa1pre
    xa = vars(1);
    T = vars(2);

    eq1 = 11835/log(((7.825e4*xa)/(1-xa))*((1-0.0475*xa)/(0.115 - 0.0475*xa))^0.5) - T;
    eq2 = 809.32 + 479.1*(xa - xa1pre) - T;

    F = [eq1;eq2];

end


function F = soe3(vars)
global xa2pre
    xa = vars(1);
    T = vars(2);

    eq1 = 11835/log(((7.825e4*xa)/(1-xa))*((1-0.0475*xa)/(0.115 - 0.0475*xa))^0.5) - T;
    eq2 = 836.56 + 479.1*(xa - xa2pre) - T;

    F = [eq1;eq2];

end

function F = soe4(vars)
global xa2pre
    xa = vars(1);
    T = vars(2);

    eq1 = 11835/log(((7.825e4*xa)/(1-xa))*((1-0.0475*xa)/(0.115 - 0.0475*xa))^0.5) - T;
    eq2 = 648 + 479.1*xa - T;

    F = [eq1;eq2];

end
