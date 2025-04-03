clc;
clear;

y = @(t) (0).*(t > 0 & t < 5) + (3.*exp((-1/4).*(t-5))).*(t > 5);

const = integral(@(t) y(t),5,inf);

E = @(t) y(t)/const;

t = linspace(0,50,101);

plot(t,E(t));

tau = integral(@(t) t.*E(t),5,inf);

[w,u] = ode45(@soe,[0 100],[2 0 0 0]);

figure
plot(w,u(:,1));
hold on
plot(w,u(:,2));
plot(w,u(:,3));
plot(w,u(:,4));

bconst = trapz(w,u(:,2));
EB = u(:,2)/bconst;
cb = trapz(w,EB.*u(:,2));

function F = soe(w,vars)

    ca = vars(1);
    cb = vars(2);
    cc = vars(3);
    cd = vars(4);

    k1 = 0.8;
    k2 = 0.1;
    k3 = 0.15;

    r1 = k1*ca;
    r2 = k2*cb;
    r3 = k2*cc;

    dca = -r1;
    dcb = r1 - r2;
    dcc = r2-r3;
    dcd = r3;

    F = [dca;dcb;dcc;dcd];

end
