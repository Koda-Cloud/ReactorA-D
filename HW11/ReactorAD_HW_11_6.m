clc;
clear;

t = [0 2 4 6 8 10 12 14 16 18 20 22 24];
co = [0 1.4 5.6 9.1 9.5 6.9 4.6 2.5 1.7 1 0.60 0.20 0];
c = co*1e-6;

x = trapz(t,c);
E = c./x;

plot(t,E);

tau = trapz(t,t.*E);
var = trapz(t,(t-tau).^2.*E);

conv = zeros(1,length(t));

for i = 1:length(t)

    dxa = @(xa) 1./(1-xa).^2;

    eq1 = @(x) t(1,i)*16.56*0.07 - integral(dxa,0,x);
    
    sol = fsolve(eq1,0.00001);

    conv(1,i) = sol;

end

xmean = trapz(t,conv.*E);
