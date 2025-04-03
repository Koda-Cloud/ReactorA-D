%clears window memory and previous figures
clc;
clear;
clf;

%establishes symbols to be used in symbolic math
syms T x

%the heat capacities of all species involved in the problems
cp_SO2 =  7.116 + 9.512e-3*x + 3.511e-6*(x^2);
cp_SO3 =  6.077 + 25.537e-3*x - 0.687e-6*(x^2);
cp_O2 =  6.148 + 3.102e-3*x - 0.923e-6*(x^2);
cp_N2 =  6.457 + 1.39e-3*x -0.069e-6*(x^2);

%the change in enthlpy for the delta H rxn
cp_rxn = -cp_SO2 + cp_SO3  - (.5*cp_O2);

%the change in enthalpy for the specieis going into the first reactor
cp_in = 11*cp_SO2 + 10*cp_O2 + 79*cp_N2;

%the complete delta H rxn for the entire problem
dhrxn =  -23490 + int(cp_rxn,[298 T]);

%output
dhrxn

%calculation for the delta H in for the energy balance to solve for extent
dhin =  int(cp_in,[T 700]);

%the extent of the reaction for the first reactor in terms of temperature
ex =  dhin/dhrxn;

%the equilibrium constant for the reaction with respect to temperature
k = exp(28.25 - 2.07*log(T) - 3.14e-7*(T^2) + 3.64e-3*T + (11511.7/T) - 27.9);

P = 1.5;

%the function that will be used to solve for the adiabatic temperature
%within the first reactor
Keq1 = k - (ex/((11-ex)*(((10 - .5*ex)/(100-.5*ex))^.5)*(P^.5))) ==0;

%numerically solving for the adiabatic temperture in the first reactor
vpasolve(Keq1)

%the cps and corresponding mol intakes for the integral of Delsta H in for
%the second reactor
cp_in2 = 3.11*cp_SO2 + 6.055*cp_O2 + 79*cp_N2 + 7.89*cp_SO3;

%Delta H in for the second reactor
dhin2 = int(cp_in2,[T 700]);

%extent of the reaction for the second reactor
ex2 = dhin2/dhrxn;

%equation to solve for the temperature in the second reactor
Keq2 = k - ((7.89+ex2)/((3.11-ex2)*(((6.055 - .5*ex2)/(96.055-.5*ex2))^.5)*(P^.5))) ==0;

%numerically solving
vpasolve(Keq2)

%the next two lines are evaluating the extent of the first reaction and
%outputting
T = 898.1;

ex1 = eval(subs(ex))

%same as a bove but for the second reactor
T = 765.3;

ex2 = eval(subs(ex2))

%calculating the conversion rates for the first and second reactors.
conv1 = ((11 - (11- ex1))/11)

conv2 = ((11- (11-ex1-ex2))/11)
