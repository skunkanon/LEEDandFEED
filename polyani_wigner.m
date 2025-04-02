function [time_span, tmp_span, rate_span] = polyani_wigner(beta,init_tmp,Ea,gam,preexponent,N_0,max_tmp)
v = preexponent; %pre-exponential
R = 8.314; %gas constant in J/(K*mol)
%N_0 = 0.05; %initial surface 
kcal_to_j = 4184; %conversion factor
Ea_0 = Ea;
y_E = gam;
%Ea_0 = 60 * 10^3; %initial activation E in J/mol

%y_E = 4 * kcal_to_j; %4kcal/mol into j/mol, eq) 14 

%defining experimental conditions
%beta = 20; % Kelvin/s
%init_tmp = 70; %initial temperature
%max_tmp = 1500;
ds_time = 20; %define desorption duration in seconds
time_span = linspace(0,(max_tmp - init_tmp)/beta,1000); %array of times from start through whole duration

%temperature as function of time 
tmp = @(t) beta*t + init_tmp; %function declaration
tmp_span = tmp(time_span); %independent variable for spectra spanning duration 

%coverage-dependent activation energy
Ea = @(N) Ea_0 - N*y_E;

%solve coverage differential equation 
dNdt = @(t, N) -v * N * min(exp(-Ea(N) / (R * tmp(t))),2);

options = odeset('RelTol', 1e-13, 'AbsTol', 1e-19); % Adjust tolerances
[t, N] = ode15s(dNdt, time_span, N_0, options);

%make coords for plot
rate_span = arrayfun(@(t, N) -dNdt(t, N), t, N);
%Ea_span = arrayfun(@(N) Ea(N),N);


[M,I] = max(rate_span);



fprintf('max rate= %d\n', max(rate_span));
fprintf('\n');
fprintf('peak temp= %d\n',tmp_span(I));



%make plots 
%plot(tmp_span,rate_span);
end

%3/19 - Maybe have this function spit back the integral to verify that
%coverage all checks out? 