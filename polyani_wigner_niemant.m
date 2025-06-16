function [time_span, tmp_span, rate_span,cov_span] = polyani_wigner_niemant(beta,init_tmp,Ea,w,preexponent,N_0,max_tmp, T_c)
v = preexponent; %pre-exponential
R = 8.314; %gas constant in J/(K*mol)
kcal_to_j = 4184; %conversion factor
Ea_0 = Ea;
%y_E = w;
time_span = linspace(0,(max_tmp - init_tmp)/beta,3000); %array of times from start through whole duration
%temperature as function of time 
tmp = @(t) beta*t + init_tmp; %function declaration
tmp_span = tmp(time_span); %independent variable for spectra spanning duration 
%coverage-dependent activation energy
%Ea = @(N) Ea_0 - N*y_E;
%solve coverage differential equation 
dNdt = @(t, N) -v * N * exp(-Ea_0 / (R * tmp(t))) * exp((w * N / R) * (1/tmp(t) - 1/T_c));
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
end

