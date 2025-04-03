clc;
clear;
w0 = linspace(0,5,1000);

[w,y] = ode45(@ex1,w0,0);

w1 = w;

plot(w1,y)
grid

function F  = ex1(w,x)

    xa = x;

    P = 10;
    k1 = 0.5;
    k2 = 0.3;
    T = 500;
    
    Pab = (1-xa)/(1.2-(0.5*xa))*P;

    PL = @(rw) Pab - rw./0.9;

    phi = @(rw)  0.5.*(PL(rw)./(0.08206*T)).^0.5;

    eta = @(rw) 1 - 0.95.*phi(rw) - 0.6.*phi(rw).^2;

    fun = @(rw) eta(rw).*k1.*(PL(rw)./(1+(k2.*PL(rw)))).^2 - rw;

    F = fsolve(fun,0.0001);

end
