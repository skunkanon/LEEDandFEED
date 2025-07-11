function [time_span, tmp_span, rate_span,cov_span] = polyani_wigner_comp1(beta,init_tmp,Ea,gam,w,preexponent,N_0,max_tmp)
%First iteration of including the compensation effect in polyani_wigner
%calculations, 6/10
v = preexponent; %pre-exponential
R = 8.314; %gas constant in J/(K*mol)
kcal_to_j = 4184; %conversion factor
Ea_0 = Ea;

y_E = gam;
time_span = linspace(0,(max_tmp - init_tmp)/beta,3000); %array of times from start through whole duration
%temperature as function of time 
tmp = @(t) beta*t + init_tmp; %function declaration
tmp_span = tmp(time_span); %independent variable for spectra spanning duration 
%coverage-dependent activation energy
Ea = @(N) Ea_0 - N*y_E;

Ea = @(N) Ea_0;
%NEW CODE, 6/10
v_theta = @(N) 10^(log10(v) + w * N);
%RESUME OLD CODE 

dNdt = @(t, N) -v_theta(N) * N * exp(-Ea(N) / (R * tmp(t)));

%RESUME OLD CODE 
%solve coverage differential equation 
%dNdt = @(t, N) -v * N * exp(-Ea(N) / (R * tmp(t)));
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-19); % Adjust tolerances
[t, N] = ode15s(dNdt, time_span, N_0, options);
%make coords for plot
rate_span = arrayfun(@(t, N) -dNdt(t, N), t, N);
cov_span = N;
%Ea_span = arrayfun(@(N) Ea(N),N);
[~,I] = max(rate_span);
fprintf('max rate= %d\n', max(rate_span));
fprintf('\n');
fprintf('peak temp= %d\n',tmp_span(I));


%NEW CODE, 6/10

%3/19 - Maybe have this function spit back the integral to verify that
%coverage all checks out? 