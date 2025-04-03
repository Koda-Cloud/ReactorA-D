clc;
clear;

eq1 = @(x,T) 5e5*exp((30e3/1.987)*((1./T)-(1/323))) - x.^2./(1-x).^2;
eq2 = @(x,T) 300 + 600.*xa - T;

fimplicit(eq1,[0 0.99 300 1000]);

view([90 -90])

global xa1 xa2

Nao = 10;

solution = fsolve(@soe1,[0.4 300]);

xa1 = 0.9*solution(1);

sol = fsolve(@soe2,[0.4 350]);

xa2 = 0.9*sol(1);

sol2 = fsolve(@soe3, [.4 350]);

xa3 = 0.9*sol2(1);

Na2 = Nao*(1-xa1)*(1-xa2);
Na3 = Nao*(1-xa1)*(1-xa2)*(1-xa3);

xaoverall2 = (Nao - Na2)/Nao;
xaoverall3 = (Nao - Na3)/Nao;

T1 = 300 + 600.*xa1;
T2 = 350 + 600*xa2;
T3 = 350 + 600*xa3;



function F = soe1(vars)

    xa = vars(1);
    T = vars(2);

    eq1 = 300 + 600.*xa - T;
    eq2 = 5e5*exp((30e3/1.987)*((1./T)-(1/323))) - xa.^2./(1-xa).^2;

    F = [eq1,eq2];

end

function F = soe2(vars)

    global xa1

    xa = vars(1);
    T = vars(2);

    eq1 = 350 + 600*(xa) - T;
    eq2 = 5e5*exp((30e3/1.987)*((1./T)-(1/323))) - (xa1 + xa.*(1-xa1)).^2./((1-xa1).*(1-xa)).^2;

    F = [eq1,eq2];

end

function F = soe3(vars)

global xa1 xa2

    xa = vars(1);
    T = vars(2);

    eq1 = 350 + 600*(xa)- T;
    eq2 = 5e5*exp((30e3/1.987)*((1./T)-(1/323))) - (xa1 + xa2*(1-xa1) + (1-xa1)*(1-xa2).*xa).^2./((1-xa1)*(1-xa2).*(1-xa)).^2;

    F = [eq1,eq2];

end
