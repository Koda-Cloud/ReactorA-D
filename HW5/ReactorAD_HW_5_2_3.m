clear;
clc;

[t,y] = ode45(@soe,[0 5000],[0.0185,0,0,650]);
plot(t,y(:,1:3));
legend('na','nb','nd');
ylabel('Molar flowrate (mol/s)');
xlabel('Volume (L)');

%plot(t,y(:,4));

[t,x] = ode45(@soe2,[0 500],[0.00344008,0.0150599,0.0150599]);
figure
plot(t,x(:,1:3));
legend('na','nb','nd');
ylabel('Molar flowrate (mol/s)');
xlabel('Volume (L)');

function F = soe(t,x)

    na = x(1);
    nd = x(2);
    ne = x(3);
    T = x(4);

    nt = na + nd + ne;
    
    R = 1.987;
    
    k = 7.94e9*exp(-27.4e3/(R*T));
    
    R1 = 8.314;
    P = 2;

    dnadv = -k*(na/nt)*(P/(R1*T));
    dnddv = k*(na/nt)*(P/(R1*T));
    dnedv = k*(na/nt)*(P/(R1*T));
    dTdv = ((-15170 + 2.2*(T - 298))/(75.2*na + 37.7*nd + 35.3*ne))*((k*na)/nt)*(P/(R1*T));

    F = [dnadv;dnddv;dnedv;dTdv];
end

function F = soe2(t,x)

    na = x(1);
    nd = x(2);
    ne = x(3);

    nt = na + nd + ne;
    
    R = 1.987;
    T = 650;
    
    k = 7.94e9*exp(-27.4e3/(R*T));
    
    R1 = 8.314;
    P = 2;

    dnadv = -k*(na/nt)*(P/(R1*T));
    dnddv = k*(na/nt)*(P/(R1*T));
    dnedv = k*(na/nt)*(P/(R1*T));

    F = [dnadv;dnddv;dnedv];
end
