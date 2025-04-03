clc;
clear;

cao = 4;

x = linspace(0,0.325);
xa = linspace(0.00001,.999);
T = 298 + 72*x;

T1 = 321.6;

T3 = 9059./(log(5.29e10.*(xa.^2./(1-xa)).*cao));

plot(x,T);
hold on
plot([0.325 0.8],[T1 T1]);
plot(xa,T3);

axis([0 1 298 600]);
