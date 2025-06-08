%some constants
v = 10^13; %pre-exponential
R = 8.314; %gas constant in J/(K*mol)
Ea_0 = 6 * 10^4; %initial activation E in J/mol
N_0 = 1; %initial surface % Define the range of k values
k_values = 0.1:0.1:1.0; % k ranges from 0.1 to 1.0 with increments of 0.1


y_E = 16.736; %4cal/mol into j/mol, eq) 14 

%defining experimental conditions
beta = 5; % Kelvin/s
init_tmp = 70; %initial temperature
ds_time = 20; %define desorption duration in seconds
t_span = linspace(0,ds_time,1000); %array of times from start through whole duration

%temperature as function of time 
tmp = @(t) beta*t + init_tmp; %function declaration
tmp_span = tmp(t_span); %independent variable for spectra spanning duration 

%coverage-dependent activation energy
Ea = @(N) Ea_0 - N*y_E;
%defining experimental conditions
beta_4 = 20; % Kelvin/s

%temperature as function of time 
tmp = @(t) beta_4*t + init_tmp; %function declaration
tmp_span = tmp(t_span); %independent variable for spectra spanning duration 

%solve coverage differential equation 
dNdt_4 = @(t,N) -v * N * exp(-Ea(N)/(R*tmp(t)));
[t,N] = ode45(dNdt_4, t_span, N_0);


%make coords for plot
dNdt4_span = arrayfun(@(t, N) dNdt_4(t, N), t, N);
Ea4_span = arrayfun(@(N) Ea(N),N);



%make plots 
plot(tmp_span,N,tmp_span,-dNdt4_span);