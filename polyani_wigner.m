function [time_span, tmp_span, rate_span,cov_span] = polyani_wigner(beta,init_tmp,Ea,gam,preexponent,N_0,max_tmp)
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
%solve coverage differential equation 
dNdt = @(t, N) -v * N * exp(-Ea(N) / (R * tmp(t)));
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-19); % Adjust tolerances
[t, N] = ode15s(dNdt, time_span, N_0, options);
%make coords for plot
rate_span = arrayfun(@(t, N) -dNdt(t, N), t, N) ./ beta;
cov_span = N;
%Ea_span = arrayfun(@(N) Ea(N),N);
[~,I] = max(rate_span);
fprintf('max rate= %d\n', max(rate_span));
fprintf('\n');
fprintf('peak temp= %d\n',tmp_span(I));
end

%3/19 - Maybe have this function spit back the integral to verify that
%coverage all checks out? 