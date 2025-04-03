clear;
clc;

[t,x] = ode45(@batch,[0 100],[1 0 0]);
plot(t,x);
legend("ca","cb","cc");
ylabel("Concentration in mol/L");
xlabel("Time in Min");

global tau

yy = zeros(10,3);

for tau = 1:1:100
    xx = fsolve(@cstr_fun, [1 1 1]);
    yy(tau,1:3) = xx;
end


function y = batch(t,x)

    ca = x(1);
    cb = x(2);
    cc = x(3);

    r1 = (0.05*ca)/(1+ca+cc);
    r2 = (0.03*cb)/(1+cb);

    y1 = -r1;
    y2 = r1-r2;
    y3 = r2;
    y = [y1;y2;y3];
end

function y = cstr_fun(x)

    global tau
    
    ca = x(1);
    cb = x(2);
    cc = x(3);
    cao = 1;

    r1 = (0.05*ca)/(1+ca+cc);
    r2 = (0.03*cb)/(1+cb);

    y1 = ca - cao + r1*tau;
    y2 = cb - r1*tau + r2*tau;
    y3 = cc - r2*tau;

    y = [y1,y2,y3];

end
