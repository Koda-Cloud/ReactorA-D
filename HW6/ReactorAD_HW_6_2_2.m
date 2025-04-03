clc;
clear;

T1 = linspace(250,500,100);
V = 0.47;
Q = 0.416e-3;

K = 4e6*exp(-7900./T1);

tau = V/Q;

xa1 = (K.*tau)./(1 + K.*tau);
xa2 = (T1 - 343)./39.76;
xa3 = (T1 - 363.8)./35.67 + .7307;
xa4 = (.7307 + K.*tau)./(1 + K.*tau);

plot(T1,xa1);
hold on
plot(T1,xa2);
plot(T1,xa3);
plot(T1,xa4);
legend('Equilibrium Conversion: CSTR 1','Tempreature Line: CSTR 1','Equilibrium Conversion: CSTR 2','Tempreature Line: CSTR 2');
ylabel('Conversion');
xlabel('Temperature (K)');
axis([250 500 -3 3]);
