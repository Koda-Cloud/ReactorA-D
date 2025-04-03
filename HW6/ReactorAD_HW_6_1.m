clc;
clear;

solution = fsolve(@soe,[0.1, 24, 220]);


function F = soe(vars)

    ca1 = vars(1);
    V = vars(2);
    T = vars(3);

    Q = 50;
    Cao = 0.25;
   
    Ca2 = 0.2*Cao;

    K = 60*exp(-2.5 - (450/(T + 460)));

    eq1 = Cao - K*ca1*(V/Q) - ca1;
    eq2 = ca1 - K*Ca2*(V/Q) - Ca2;
    eq3 = 150 + (20000/1875)*K*ca1*V - T;

    F = [eq1;eq2;eq3];
end
