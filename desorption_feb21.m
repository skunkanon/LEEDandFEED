%{
2/21 notes:
ultimately should be able to simulate spectra as per falconer-maddix
method with the g(c) and everything; take real data, find activation E and
pre-exponential factor. 

2/23: 

pretty straightforward to get index of max dNdt for peak
temperature. dont know why it starts exploding at higher beta. paper said
35 K/s is max but can't handle even 15. 

%}

%some constants
v = 10^13; %pre-exponential
R = 8.314; %gas constant in J/(K*mol)
N_0 = 0.6; %initial surface 
kcal_to_j = 4184; %conversion factor
Ea_0 = 60 * 10^3; %initial activation E in J/mol

y_E = 4 * kcal_to_j; %4kcal/mol into j/mol, eq) 14 

%defining experimental conditions
beta = 20; % Kelvin/s
init_tmp = 70; %initial temperature
ds_time = 20; %define desorption duration in seconds
t_span = linspace(0,ds_time,100); %array of times from start through whole duration

%temperature as function of time 
tmp = @(t) beta*t + init_tmp; %function declaration
tmp_span = tmp(t_span); %independent variable for spectra spanning duration 

%coverage-dependent activation energy
Ea = @(N) Ea_0 - N*y_E;

%solve coverage differential equation 
dNdt = @(t, N) -v * N * min(exp(-Ea(N) / (R * tmp(t))), 2);

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9); % Adjust tolerances
[t, N] = ode23s(dNdt, t_span, N_0, options);

%make coords for plot
rate_span = arrayfun(@(t, N) -dNdt(t, N), t, N);
%Ea_span = arrayfun(@(N) Ea(N),N);


[M,I] = max(rate_span);



fprintf('max rate= %d\n', max(rate_span));
fprintf('\n');
fprintf('peak temp= %d\n',tmp_span(I));



%make plots 
plot(tmp_span,rate_span);

%plot(t_span,Ea_span);
%plot(tmp_span,-dNdt_span);

%% 2/25 - trying to replicate spectra from given Tp, heating rates - falconer & madix 

Ea_0 = 33.25 * kcal_to_j;
%Ea_0 = 60 * 10^3;
v = 8.5 * 10^15;

%coverage-dependent activation energy
Ea = @(C) Ea_0 - C*y_E;



%defining experimental conditions
beta = 35; % Kelvin/s
init_tmp = 300; %initial temperature
ds_time = 20; %define desorption duration in seconds
t_span = linspace(0,ds_time,10000); %array of times from start through whole duration

%temperature as function of time 
tmp = @(t) beta*t + init_tmp; %function declaration
tmp_span = tmp(t_span); %independent variable for spectra spanning duration 

C_0 = 1;
f = 0.9;
g = @(C) (C/C_0)*(C_0 - C + f*C_0);

dCdt = @(t,C) -v * g(C) * exp(-Ea(C)/(R*tmp(t))); % = -N, the rate of desorption
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9); % Adjust tolerances
[t,C] = ode23s(dCdt, t_span,C_0,options);


dCdt_span = arrayfun(@(t,C) dCdt(t,C), t,C);
Ea_span = arrayfun(@(C) Ea(C),C);
g_span = arrayfun(@(C) g(C), C);

tmp_C_span = tmp(t_span) - 273;
%figure
plot(tmp_C_span, C, 'b', tmp_C_span, -dCdt_span);
legend('coverage', 'rate');
%% 

figure

plot(t, tmp_span, 'r');
legend('temp(time)'); 

figure 
plot(tmp_span, C);



%% different heating rate 'beta'

%defining experimental conditions
beta_2 = 10; % Kelvin/s

%temperature as function of time 
tmp = @(t) beta_2*t + init_tmp; %function declaration
tmp_span = tmp(t_span); %independent variable for spectra spanning duration 

%solve coverage differential equation 
dNdt_2 = @(t,N) -v * N * exp(-Ea(N)/(R*tmp(t)));
[t,N] = ode45(dNdt_2, t_span, N_0);


%make coords for plot
dNdt2_span = arrayfun(@(t, N) dNdt_2(t, N), t, N);
Ea2_span = arrayfun(@(N) Ea(N),N);

%{                                      another instance                         %}

%defining experimental conditions
beta_3 = 15; % Kelvin/s

%temperature as function of time 
tmp = @(t) beta_3*t + init_tmp; %function declaration
tmp_span = tmp(t_span); %independent variable for spectra spanning duration 

%solve coverage differential equation 
dNdt_3 = @(t,N) -v * N * exp(-Ea(N)/(R*tmp(t)));
[t,N] = ode45(dNdt_3, t_span, N_0);


%make coords for plot
dNdt3_span = arrayfun(@(t, N) dNdt_3(t, N), t, N);
Ea3_span = arrayfun(@(N) Ea(N),N);

%{                                      another instance                         %}

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


%{                                      another instance                         %}

%{  

BELOW EXPLODES WHAT

beta_5 = 25; 


%temperature as function of time 
tmp = @(t) beta_5*t + init_tmp; %function declaration
tmp_span = tmp(t_span); %independent variable for spectra spanning duration 

%solve coverage differential equation 
dNdt_5 = @(t,N) -v * N * exp(-Ea(N)/(R*tmp(t)));
[t,N] = ode45(dNdt_5, t_span, N_0);


%make coords for plot
dNdt5_span = arrayfun(@(t, N) dNdt_5(t, N), t, N);
Ea5_span = arrayfun(@(N) Ea(N),N);

%}




plot(tmp_span,-dNdt_span,tmp_span,-dNdt2_span, tmp_span, -dNdt3_span, tmp_span, -dNdt4_span);
%% 

[town,car] = equationfit_test_feb25(1,1);

fprintf('best gen=%d\n',town);
fprintf('RIP in=%d\n',car);
%% 

fig5_N = [0.005, 0.013, 0.018, 0.041];


x_5 = cell(1, length(fig5_N)); 
y_5 = cell(1, length(fig5_N)); 


for i = 1:length(fig5_N)
    [x_5{i}, y_5{i}] = desorption_sim(fig5_N(i));
end


figure;
hold on;
for i = 1:length(fig5_N)
    plot(x_5{i}, y_5{i}, 'DisplayName', sprintf('N = %.3f', fig5_N(i)));
end
hold off;


xlabel('Temperature (K)');
ylabel('-dN/dt');
title('Desorption Rate vs. Temperature for Different N Values');
legend show;
grid on;



fig6_N = [0.026,0.048,0.11,0.2,0.32,0.41,0.66,0.95,1];

x_6 = cell(1, length(fig6_N)); 
y_6 = cell(1, length(fig6_N)); 


for i = 1:length(fig6_N)
    [x_6{i}, y_6{i}] = desorption_sim(fig6_N(i));
end


figure;
hold on;
for i = 1:length(fig6_N)
    plot(x_6{i}, y_6{i}, 'DisplayName', sprintf('N = %.3f', fig6_N(i)));
end
hold off;


xlabel('Temperature (K)');
ylabel('-dN/dt');
title('Desorption Rate vs. Temperature for Different N Values');
legend show;
grid on;


%% EQ 16 IMPLEMENTATION  - VARY BETA ONLY, BUILDING OFF CURVE 'a' IN FIGURE 5
% 2/28 - PEAK TEMP FOR BETA = 9.0 K/s DOESN'T AGREE WITH DESORPTION PLOT
% NOR PAPER TWO SECTIONS UP, IDK WHY BUT WE'LL KEEP GOING. ODD BEHAVIOR FOR
% BETA LOWER THAN ~5 (0.35 GIVES A PEAK TEMP OF 34 CELSIUS, SHOULD BE AROUND 100 FOR FIGURE 1), BUT TAKES BETA
% HIGHER THAN 35 LIKE A CHAMP. 'f' SHOULD STAY CONSTANT I THINK

fig5_beta = linspace(0.2, 35,18); 

e16_x_5 = cell(1, length(fig5_beta)); 
e16_y_5 = cell(1, length(fig5_beta)); 
e16_t_5 = cell(1,length(fig5_beta));
tempindex = zeros(1, length(fig5_beta));
e16_Tp_5 = zeros(1, length(fig5_beta));
int_list = zeros(1,length(fig5_beta));

fprintf("\n");

fprintf("FIGURE 5");
for i = 1:length(fig5_beta)
    [e16_t_5{i},e16_x_5{i}, e16_y_5{i}] = madix_sim_eq16(fig5_beta(i));
    [~, tempindex(i)] = max(e16_y_5{i});
    e16_Tp_5(i) = e16_x_5{i}(tempindex(i));
    fprintf("Peak Temp for B = %.4f\n",fig5_beta(i));
    fprintf("   %.2f C\n\n",e16_Tp_5(i)-273);
end


figure;
hold on;
%{
scale1 = 1.5;

plot(e16_x_5{1},scale1*e16_y_5{1},'DisplayName', sprintf('Beta = %.3f', fig5_beta(1))); %SCALE LOWER BETA TO BE VISIBLE ON THE PLOT, "'Scale Factor =%.3f',scale1," SCREWS UP THE PLOT 

plot(e16_x_5{1},1*e16_y_5{1},'DisplayName', sprintf('Beta = %.3f', fig5_beta(2)));

plot(e16_x_5{3},e16_y_5{3},'DisplayName', sprintf('Beta = %.3f', fig5_beta(3)));

plot(e16_x_5{4},e16_y_5{4},'DisplayName', sprintf('Beta = %.3f', fig5_beta(4)));

hold off;
%}




for i = 1:length(fig5_beta)
    plot(e16_x_5{i}, e16_y_5{i}, 'DisplayName', sprintf('Beta = %.3f', fig5_beta(i)));
    int_list(i) = trapz(e16_t_5{i},e16_y_5{i});
    fprintf('%.4\n', int_list(i));
end
hold off;

for i = 1:length(fig5_beta)


end

xlabel('Temperature (K)');
ylabel('Rate = -dN/dt');
title('Desorption Rate vs. Temperature for Different Beta Values, N = 0.005');
legend show;
grid on;


kcal_to_j = 4184;
R = 8.314;
Ea_0 = 32 * kcal_to_j;

tempx = zeros(1, length(fig5_beta));
tempy = zeros(1, length(fig5_beta));






for i = 1:length(fig5_beta)
    %fprintf("%d\n",log(fig5_beta(i)/e16_Tp_5(i)^2));
    tempx(i) = 1/e16_Tp_5(i);
    tempy(i) = log(fig5_beta(i)/e16_Tp_5(i)^2); 
    
end

%}
p = polyfit(tempx, tempy, 1);

% Extract the slope (m) and intercept (b) of the line
m = p(1); % Slope
b = p(2); % Intercept

% Generate fitted y-values
fitted_y = polyval(p, tempx);

% Plot the original data and the fitted line
figure;
plot(tempx, tempy, 'bo', 'DisplayName', 'Data'); % Original data

hold on;

plot(tempx, fitted_y, 'r-', 'DisplayName', sprintf('Fitted Line: y = %.2fx + %.2f', m, b)); % Fitted line
hold off;

fprintf('Ea = %.2f',(-m*R)/kcal_to_j);
fprintf(' kcal/mol \n');
% Add labels and legend
xlabel('tempx');
ylabel('tempy');
title('Linear Regression');
legend show;
grid on;

%3/8 - Int_list checks out just fine - all returning the same as the
%initial coverage of 0.005. Yet to check if it works for different
%coverages but not much of a reason to think they don't. 
%% 

%STUFF FOR EXCEL, TROUBLESHOOTING ABOVE SECTION, 3/4 

ex_beta = fig5_beta';
ex_Tp = e16_Tp_5';
ex_tempindex = tempindex';
ex_x = tempx';
ex_y = tempy';

%both agree that M = -2009, and Ea = ~4 kcal/mol, way off from 32. 
%3/4 note; beta  10 to 35 returns the right Ea within 0.02 kcal. Ea drops
%with a heating rate lower than that. Eq 16 turns nonlinear for heating
%rates lower than 10 K/s. Something's up with the solver. 

%Yep, got it. Temperature didn't ramp up to 500 K in each case. 

%From 1 kcal/mol, up to 100 (didn't bother testing higher), consistently
%lowballed actual Ea by 0.04-0.05. 

%Changing the frequency factor 'v' doesn't do much. 

%f = 0.9 gives about the closest for 32 kcal.


%% Now going back and fixing replication of Fig 6, 3/4. Slight overestimation at coverages below 0.048. 



fig5_N = [0.005, 0.013, 0.018, 0.041];
fig5_Ea0 = [33.25, 33.16, 33.01, 33.01];


x_5 = cell(1, length(fig5_N)); 
y_5 = cell(1, length(fig5_N)); 
tempindex = zeros(1, length(fig5_N));
Tp_5 = zeros(1, length(fig5_N));

fprintf("\n");

fprintf("FIGURE 5");

fprintf("\n");

for i = 1:length(fig5_N)
    [x_5{i}, y_5{i}] = madix_sim(fig5_N(i),fig5_Ea0(i));
    [~, tempindex(i)] = max(y_5{i});
    Tp_5(i) = x_5{i}(tempindex(i));
    fprintf("Peak Temp for N = %.4f\n",fig5_N(i));
    fprintf("   %.2f C\n\n",Tp_5(i)-273);
end


figure;
hold on;
for i = 1:length(fig5_N)
    plot(x_5{i}, y_5{i}, 'DisplayName', sprintf('N = %.3f', fig5_N(i)));
end
hold off;


xlabel('Temperature (K)');
ylabel('Rate = -dN/dt');
title('Desorption Rate vs. Temperature for Different N Values');
legend show;
grid on;
%% 3/7 - Now re-replicating fig 6 and 4. Peak temp significantly underestimated at coverages above 0.048.


fig6_N = [0.026, 0.048,  0.11, 0.2,  0.32,  0.41,  0.66,  0.95,   1];

fig6_Ea = [33.01, 32.78, 32.7, 32.3, 32.22, 32.15, 32.07, 32.07, 32];

x_6 = cell(1, length(fig6_N)); 
y_6 = cell(1, length(fig6_N)); 
time_6 = cell(1, length(fig6_N)); 
tempindex = zeros(1, length(fig6_N));
Np_6 = zeros(1, length(fig6_N));
Tp_6 = zeros(1, length(fig6_N));
%time_to_peak = zeros(1,length(fig6_N)); %not useful yet, refer to 3/7
%comment 

fprintf("\n");

fprintf("FIGURE 6");

fprintf("\n");

for i = 1:length(fig6_N)
    [time_6{i},x_6{i}, y_6{i}] = madix_sim(fig6_N(i),fig6_Ea(i));
    [Np_6(i), tempindex(i)] = max(y_6{i});
    Tp_6(i) = x_6{i}(tempindex(i));

    fprintf("Peak Temp for N = %.4f\n",fig6_N(i));
    fprintf("   %.2f C\n\n",Tp_6(i)-273);
end


figure;
hold on;
for i = 1:length(fig6_N)
    plot(x_6{i}, y_6{i}, 'DisplayName', sprintf('N = %.3f', fig6_N(i))); %plots the desorption curve family 
    
end
hold off;


xlabel('Temperature (K)');
ylabel('Rate = -dN/dt');
title('Desorption Rate vs. Temperature for Different N Values');
legend show;
grid on;




%{

%3/7 - ignore this N should be max desorption RATE not the coverage at
the point it hits that, but could be useful later 

for i = 1:length(fig6_N)
    time_to_peak(i) = time_6{1,i}(tempindex(i)); %gets the time at which it hits the max desorption
end

N_to_peak = zeros(1,length(fig6_N));

for i = 1:length(fig6_N)
    lol = linspace(0,time_to_peak(i),10000); %totally unoptimized code that returns the integral of the desorption peak up until the peak 
    lmao = trapz(lol,y_6{i});
    N_to_peak(i) = lmao;
    %fprintf('%.4f',lmao);
end

%}

equa16_x = zeros(1,length(fig6_N));
equa16_y = zeros(1,length(fig6_N));

for i = 1:length(fig6_N)
    equa16_x(i) = 1/Tp_6(i);
    equa16_y(i) = log(Np_6(i));
   
end

%equa16_y = fliplr(equa16_y);


p = polyfit(equa16_x, equa16_y, 1); 
y_fit = polyval(p, equa16_x); 



figure;

scatter(equa16_x,equa16_y);
%plot(equa16_x, equa16_y, 'bo', 'DisplayName', 'Data Points'); 
hold on;
%plot(equa16_x, y_fit, 'r-', 'DisplayName', 'Linear Regression'); 
hold off;

xlabel('1 / T_p (1/K)');
ylabel('ln(N_{p})');
title('ln(N_{p}) vs. 1/T_p');
legend show;
grid on;


fprintf('Linreg: y = %.4f * x + %.4f\n', p(1), p(2));

R = 8.314;

kcal_to_j = 4184;


fprintf('Ea = %.4f\n', p(1)*R*kcal_to_j);

% 3/7: Now is a bit more accurate after reducing y_E (activation energy's
% dependence on coverage) by a factor of 10. Now overestimates for N >
% 0.11, and slightly underestimates for N < 0.11. For instance, Tp at
% saturation for y_E = 4 kcal was 118.81, but at y_E at 0.4 kcal was
% 158.01. 