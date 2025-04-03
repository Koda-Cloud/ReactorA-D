clear;
clc;

[t,y] = ode45(@soe,[0 240],[0.8,0,0,0,320]);
plot(t,y(:,1:4));
legend('ca','cb','cc','cd');
figure
plot(t,y(:,5));


function F = soe(t,x)

    ca = x(1);
    cb = x(2);
    cc = x(3);
    cd = x(4);
    T = x(5);

    k1 = 0.04*exp((-7000/1.987)*((1/T) - (1/320)));
    k2 = 0.03*exp((-6000/1.987)*((1/T) - (1/320)));

    r1 = k1*ca;
    r2 = k2*cb;

    dcadt = -r1;
    dcbdt = 2*r1 - r2;
    dccdt = 2*r2;
    dcddt = r2;
    dTdt = -40*r1 - 60*r2;

    F = [dcadt;dcbdt;dccdt;dcddt;dTdt];
end
