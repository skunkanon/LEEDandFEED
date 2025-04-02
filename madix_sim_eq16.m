function [t_span,tmp_span,rate_span] = madix_sim_eq16(beta)
%LOW COVERAGE 


R = 8.314; %gas constant in J/(K*mol)
C_0 = 0.005; %initial surface 
kcal_to_j = 4184; %conversion factor

Eact = 30;
Ea_0 = Eact * kcal_to_j;

v = 8.5 * 10^15;

%coverage-dependent activation energy
y_E = 4 * kcal_to_j; %4kcal/mol into j/mol, eq) 14 
Ea = @(C) Ea_0 - C*y_E;



%defining experimental conditions
%beta = 9; % Kelvin/s
init_tmp = 0; %initial temperature
max_tmp = 500;
t_span = linspace(0,max_tmp/beta,10000); %array of times from start through whole duration

%temperature as function of time 
tmp = @(t) beta*t + init_tmp; %function declaration
tmp_span = tmp(t_span); %independent variable for spectra spanning duration 

%C_0 = 1;
f = 0.9;
g = @(C) (C/C_0)*(C_0 - C + f*C_0);

dCdt = @(t,C) -v * g(C) * exp(-Ea(C)/(R*tmp(t))); % = -N, the rate of desorption
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-19); % Adjust tolerances
[t,C] = ode15s(dCdt, t_span,C_0,options);


dCdt_span = arrayfun(@(t,C) dCdt(t,C), t,C); %IS NEGATIVE, returned multiplied by -1
%Ea_span = arrayfun(@(C) Ea(C),C);
%g_span = arrayfun(@(C) g(C), C);

rate_span = -dCdt_span;

end