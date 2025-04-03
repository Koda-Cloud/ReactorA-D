clc;
clear;

[t, y] = ode45(@soe,[0 50],0.549);

plot(t,y(:,1));

[t1,y1] = ode45(@soe2,[0 50],0.766);

figure
plot(t1,y1(:,1));

function F = soe(t,x)
    
    ca = x;
  
    cao = 1;
    tau = 0.5;

    dca = cao/tau  - ca/tau - 3*exp(-0.2*t)*ca^2;

    F = dca;

end

function F = soe2(t,x)
    
    ca = x;

    cao = 1;
    tau = 0.5;

    phi = 10*exp(-0.003*t)*(ca)^0.5;
    eta = (1.0357 + 0.3173*phi + 0.000437*phi^2)/(1 + 0.4172*phi + 0.139*phi^2);

    dca = cao/tau - ca/tau - eta*3*exp(-0.2*t)*ca^2;

    F = dca;

end
