clc;
clear;

global sol

xa = 0.62;
R = 8.314;
P = 1;  
Nto = 2.1;
E1 = 41570;
E2 = 83140;
Cp = 42;
    
Nao = .10*Nto;
    
k1 = @(T) 1e3*exp(-41570./(R.*T));
k2 = @(T) 1e6*exp(-83140./(R.*T));
    
eq1 = @(T) (41570.*k1(T).*(1-xa) - 83140.*k2(T).*xa);

sol = fsolve(eq1,573);

solution = fsolve(@soe,[0.5,650,20]);

% Q = -88.2*(573-sol) + (41570-83140)*((P*Nao)/Nto)*(k1(sol)*(1-xa) - k2(sol)*xa)*V;

eq2 = @(T)  -(E1 - E2)*(P/Nto)*(k1(T)*(1-solution(1)) - k2()*solution(1))*solution(3);



function F = soe(vars)

    global sol

    xa1 = vars(1);
    T = vars(2);
    V = vars(3);

    xa2 = 0.62;
    R = 8.314;
    P = 1;  
    Nto = 2.1;

    k1 = 1e3*exp(-41570./(R.*T));
    k2 = 1e6*exp(-83140./(R.*T));

    k11 = 1e3*exp(-41570./(R.*sol));
    k22 = 1e6*exp(-83140./(R.*sol));
    
    eq1 = (41570.*k1.*(1-xa1) - 83140.*k2.*xa1);
    eq2 = (P/Nto)*(k1*(1-xa1) - k2*xa1)*V - xa1;
    eq3 = xa1 + (P/Nto)*(k11*(1-xa2) - k22*xa2)*5 - xa2;

    F = [eq1,eq2,eq3];
end
