clc;
clear;


xa = 0.62;
R = 8.314;
P = 1;  
Nto = 2.1;
    
Nao = .10*Nto;
    
k1 = @(T) 1e3*exp(-41570./(R.*T));
k2 = @(T) 1e6*exp(-83140./(R.*T));
    
eq1 = @(T) (41570.*k1(T).*(1-xa) - 83140.*k2(T).*xa);

sol = fsolve(eq1,573);

V = xa/((P/Nto)*(k1(sol)*(1-xa) - k2(sol)*xa)); 

Q = -88.2*(573-sol) + (41570-83140)*((P*Nao)/Nto)*(k1(sol)*(1-xa) - k2(sol)*xa)*V;

eq = @(T) 41570*5*((P*0.1*Nto)/Nto)*(k1(618)*0.38 - k2(618)*0.62) + 42*2.1*(T - 618);

sol2 = fsolve(eq, 600);

eq2 = @(T) 41570*7.24*((P*0.1*Nto)/Nto)*(k1(664)*0.52 - k2(664)*0.48) + 42*2.1*(T - 664);

sol3 = fsolve(eq2,600);

sol4 = fsolve(@soe2,[0.5,700,6]);

solution = fsolve(@soe,[650, 22]);





function F = soe(vars)

    T = vars(1);
    V = vars(2);

    xa = 0.62;
    R = 8.314;
    P = 1;  
    Nto = 2.1;

    Nao = .10*Nto;

    k1 = 1e3*exp(-41570./(R.*T));
    k2 = 1e6*exp(-83140./(R.*T));

    eq1 = (P/Nto).*((41570./(R.*T.^2)).*k1.*Nao.*(1-xa) - (83140./(R.*T.^2)).*k2.*Nao*xa);
    eq2 = ((P*V)/(Nto)).*(k1*(1-xa) - k2*xa) - xa;

    F = [eq1,eq2];

end

function F = soe2(vars)

    xa1 = vars(1);
    T = vars(2);
    V = vars(3);

    xa2 = 0.62;
    R = 8.314;
    P = 1;  
    Nto = 2.1;

    k1 = 1e3*exp(-41570./(R.*T));
    k2 = 1e6*exp(-83140./(R.*T));
    k11 = 1e3*exp(-41570./(R.*618));
    k22 = 1e6*exp(-83140./(R.*618));

    eq1 = (P/Nto)*(k1*(1-xa1) - k2*xa1)*5 - xa1;
    eq2 = xa1 + (P/Nto)*(k11*(1-xa2) - k22*xa2)*V - xa2;
    eq3 = 41570*k1*(1-xa1) - 83140*k2*xa1;

    F = [eq1,eq2,eq3];

end
