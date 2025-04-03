clc;
clear;

nebo = 707/3600;
nso = 7.104/3600;
nbo = 0.293/3600;
ntolo = 4.968/3600;
nh2oo = 7777/3600;

[w,y] = ode45(@soe,[0 6e5],[nebo nso nbo ntolo nh2oo 0 0 0 0 0 950]);

figure
plot(w,y(:,1));
hold on
plot(w,y(:,2));
legend('NEB','NS');
ylabel('Molar Flowrate (kmol/s)');
xlabel('Mass of Catalyst (kg)');

figure
plot(w,y(:,11));


function F = soe(w,x)

    neb = x(1);
    ns = x(2);
    nb = x(3);
    ntol = x(4);
    nh20 = x(5);
    nh2 = x(6);
    nc2h4 = x(7);
    nco = x(8);
    nco2 = x(9);
    nch4 = x(10);
    T = x(11);

    P = 1.25;

    ntot = neb + ns + nb + ntol + nh20 + nh2 + nc2h4 + nco + nco2 + nch4;

    peb = (neb/ntot)*P;
    ps = (ns/ntot)*P;
    ph2 = (nh2/ntot)*P;
    ph20 = (nh20/ntot)*P;
    pc2h4 = (nc2h4/ntot)*P;
    pch4 = (nch4/ntot)*P;
    pco = (nco/ntot)*P;

    k1  = 1.51286*exp(-10925/T);
    k2 = 5.6197e5*exp(-25e3/T);
    k3 = 1.3446*exp(-11e3/T);
    k4 = 9.3016e-1*exp(-12.5e3/T);
    k5 = 5.3163e-2*exp(-7.9e3/T);
    k6 = 1.6769e9*exp(-8.85e3/T);
    kp = 8.2e6*exp(-15.2e3/T);

    cpeb = -43.091 + 7.07e-1*T - 4.810e-4*T^2 + 1.301e-7*T^3; 
    cps = -28.243 + 6.158e-1*T - 4.022e-4*T^2 + 9.933e-8*T^3; 
    cpb = -33.911 + 4.743e-1*T - 3.016e-4*T^2 + 7.129e-8*T^3; 
    cptol = -24.350 + 5.124e-1*T - 2.765e-4*T^2 + 4.910e-8*T^3; 
    cpc2h4= 3.805 + 1.566e-1*T - 8.347e-5*T^2 + 1.755e-8*T^3; 
    cpch4 = 19.247 + 5.212e-2*T + 1.197e-5*T^2 - 1.131e-8*T^3;
    cph2o = 32.236 + 1.923e-3*T + 1.055e-5*T^2 - 3.596e-9*T^3; 
    cpco = 10.863 - 1.285e-2*T + 2.789e-5*T^2 - 1.271e-8*T^3; 
    cpco2 = 19.791 + 7.342e-2*T - 5.601e-5*T^2 + 1.715e-8*T^3; 
    cph2 = 27.138 + 9.272e-3*T - 1.381e-5*T^2 + 7.644e-9*T^3;

    ncp = neb*cpeb + ns*cps + nb*cpb + ntol*cptol + nc2h4*cpc2h4 + ...
        nch4*cpch4 + nh20*cph2o + nco*cpco + nco2*cpco2 + nh2*cph2;

    r1 = k1*(peb - (ps*ph2)/kp);
    r2 = k2*peb;
    r3 = k3*peb*ph2;
    r4 = k4*ph20*pc2h4;
    r5 = k5*ph20*pch4;
    r6 = k6*(P/T^3)*ph20*pco;

    dh1 = 120737 + 4.56*T;
    dh2 = 108802 - 7.95*T;
    dh3 = -53171 -13.19*T;
    dh4 = 82054 + 8.83*T;
    dh5 = 211979 + 1658*T;
    dh6 = -45217 + 10.46*T;

    dhr = dh1*r1 + dh2*r2 + dh3*r3 + dh4*r4 + dh5*r5 + dh6*r6;

    deb = -r1 - r2 - r3;
    ds = r1;
    dh2 = r1 - r3 + 2*r4 + 3*r5 + r6;
    db = r2;
    dc2h4 = r2 - 0.5*r4;
    dtol = r3;
    dch4 = r3 - r5;
    dh20 = -r4 - r5 - r6;
    dco = r4 + r5 - r6;
    dco2 = r6;
    dT = -dhr/ncp;

    F = [deb;ds;db;dtol;dh20;dh2;dc2h4;dco;dco2;dch4;dT];

end
