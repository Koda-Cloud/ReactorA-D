clc;
clear;

w0 = linspace(0,5,100);

[w,y] = ode45(@ex1,w0,0);

w1 = w*8000;

plot(w1,y)
grid

function F  = ex1(w,x)

    xa = x;

    P = 10;
    k1 = 0.5;
    k2 = 0.3;

    Pab = ((0.8*(1-xa))/(1-(0.4*xa)))*P;

    fun = @(rw) k1*((Pab - rw/0.9)/(1 + k2*(Pab - rw/0.9)))^2 - rw;

    F = fsolve(fun,0.0001);

end
