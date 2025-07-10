%Testing out new reduceBGgetspectra() function. Adapted from CLAA_BG_reduce. 6/4
%Run CL_figure5_ACTUAL first.
new_x_0p4_raw = new_fig5_0p4(:,1)'; 
new_y_0p4_raw = new_fig5_0p4(:,2)'; 
new_x_0p8_raw = new_fig5_0p8(:,1)'; 
new_y_0p8_raw = new_fig5_0p8(:,2)'; 
new_x_1p2_raw = new_fig5_1p2(:,1)'; 
new_y_1p2_raw = new_fig5_1p2(:,2)'; 
new_x_1p6_raw = new_fig5_1p6(:,1)'; 
new_y_1p6_raw = new_fig5_1p6(:,2)'; 
new_x_2p0_raw = new_fig5_2p0(:,1)'; 
new_y_2p0_raw = new_fig5_2p0(:,2)'; 
new_x_2p8_raw = new_fig5_2p8(:,1)'; 
new_y_2p8_raw = new_fig5_2p8(:,2)'; 
new_x_4p0_raw = new_fig5_4p0(:,1)'; 
new_y_4p0_raw = new_fig5_4p0(:,2)';
new_x_8p0_raw = new_fig5_8p0(:,1)'; 
new_y_8p0_raw = new_fig5_8p0(:,2)'; 
%%


if true
figure(4); clf;
hold on;
scatter(new_fig5_0p4(:,1)', new_fig5_0p4(:,2)');
scatter(new_fig5_0p8(:,1)', new_fig5_0p8(:,2)');
scatter(new_fig5_1p2(:,1)', new_fig5_1p2(:,2)');
scatter(new_fig5_1p6(:,1)', new_fig5_1p6(:,2)');
scatter(new_fig5_2p0(:,1)', new_fig5_2p0(:,2)');
scatter(new_fig5_2p8(:,1)', new_fig5_2p8(:,2)');
scatter(new_fig5_4p0(:,1), new_fig5_4p0(:,2)');
scatter(new_fig5_8p0(:,1)', new_fig5_8p0(:,2)');
hold off; 
end


%%


figure(3); clf
hold on;
[tempSPAN_actual_0p4 , signalSPAN_actual_0p4, std_0p4] = reduceBGgetspectra(new_x_0p4_raw, new_y_0p4_raw,800,1500,'r');
[tempSPAN_actual_0p8 , signalSPAN_actual_0p8, std_0p8] = reduceBGgetspectra(new_x_0p8_raw, new_y_0p8_raw,800,1500,'g');
[tempSPAN_actual_1p2 , signalSPAN_actual_1p2, std_1p2] = reduceBGgetspectra(new_x_1p2_raw, new_y_1p2_raw,800,1500,'b');
[tempSPAN_actual_1p6 , signalSPAN_actual_1p6, std_1p6] = reduceBGgetspectra(new_x_1p6_raw, new_y_1p6_raw,800,1500,'m');
[tempSPAN_actual_2p0 , signalSPAN_actual_2p0, ~] = reduceBGgetspectra(new_x_2p0_raw, new_y_2p0_raw,800,1400,'k');
[tempSPAN_actual_2p8 , signalSPAN_actual_2p8, ~] = reduceBGgetspectra(new_x_2p8_raw, new_y_2p8_raw,800,1400,'y');
[tempSPAN_actual_4p0 , signalSPAN_actual_4p0, ~] = reduceBGgetspectra(new_x_4p0_raw, new_y_4p0_raw,800,1400,'c');
[tempSPAN_actual_8p0 , signalSPAN_actual_8p0, ~] = reduceBGgetspectra(new_x_8p0_raw, new_y_8p0_raw,800,1400,'r');

%Can add more of the traces, but those left lie at the precipice or past the linear regime. 

hold off;


% Testing out new getcoverageplotcoverage() function. Adapted from 'erley80_test'
fprintf('NEW INSTANCE \n');
temp_init = 300; %K 
beta = 50; %K/s
time = @(x) (x-temp_init)/beta;



N0_0p4 = 0.0908;
N0_0p8 = 0.2004;
N0_1p2 = 0.2877;
N0_1p6 = 0.3856; %5/21 
N0_2p0 = 0.43065;
N0_2p8 = 0.46168;
N0_4p0 = 0.47695;
N0_8p0 = 0.49657;

figure(2); clf;
hold on;
[dNdt_0p4,N_0p4] = getcoverageplotcoverage(time(tempSPAN_actual_0p4), signalSPAN_actual_0p4, N0_0p4,'r');
[dNdt_0p8, N_0p8] = getcoverageplotcoverage(time(tempSPAN_actual_0p8), signalSPAN_actual_0p8, N0_0p8,'g');
[dNdt_1p2, N_1p2] = getcoverageplotcoverage(time(tempSPAN_actual_1p2), signalSPAN_actual_1p2, N0_1p2,'b');
[dNdt_1p6, N_1p6] = getcoverageplotcoverage(time(tempSPAN_actual_1p6), signalSPAN_actual_1p6, N0_1p6, 'm');


hold off;


%% 6/30 - Getting nonlinear regime desorption spectra
[dNdt_2p8, N_2p8] = getcoverageplotcoverage(time(tempSPAN_actual_2p8), signalSPAN_actual_2p8, N0_2p8, 'k');


%% Using new functions to get more accurate values for time, signal, and coverage. 6/9. 
% Revisiting Arrhenius method, 6/26 
CL_theta = 0.15; 
temp_init = 300; %K 
beta = 50; %K/s
time = @(x) (x-temp_init)/beta;
[CL_time_0p8, coverage_0p8_theta_0p15] = getpointgetinterpolation_focusYgetX(time(tempSPAN_actual_0p8), N_0p8, CL_theta, 2);
[~, rate_0p8_theta_0p15] = getpointgetinterpolation_focusXgetY(time(tempSPAN_actual_0p8), dNdt_0p8, CL_time_0p8, 2);

[CL_time_1p2, coverage_1p2_theta_0p15] = getpointgetinterpolation_focusYgetX(time(tempSPAN_actual_1p2), N_1p2, CL_theta, 2);
[~, rate_1p2_theta_0p15] = getpointgetinterpolation_focusXgetY(time(tempSPAN_actual_1p2),dNdt_1p2, CL_time_1p2, 2);


[CL_time_1p6, coverage_1p6_theta_0p15] = getpointgetinterpolation_focusYgetX(time(tempSPAN_actual_1p6), N_1p6, CL_theta, 2);
[~, rate_1p6_theta_0p15] = getpointgetinterpolation_focusXgetY(time(tempSPAN_actual_1p6), dNdt_1p6, CL_time_1p6, 2);
%

%Arrhenius analysis. Still have to fill in the values manually. 6/4

temp = @(t) beta*t + temp_init;

arrh_x = [1/temp(CL_time_0p8), 1/temp(CL_time_1p2), 1/temp(CL_time_1p6)];
%beta = 1;
arrh_y = [log(1 * rate_0p8_theta_0p15 / coverage_0p8_theta_0p15), log(1 * rate_1p2_theta_0p15 / coverage_1p2_theta_0p15) ...
    log(1 * rate_1p6_theta_0p15 / coverage_1p6_theta_0p15)];

%


%{
arrh_x = [ 1/temp(CL_time_1p2), 1/temp(CL_time_1p6)];
%beta = 1;
arrh_y = [log(1 * rate_1p2_theta_0p15 / 0.2) ...
    log(1 * rate_1p6_theta_0p15 / 0.2)];
%}


p = polyfit(arrh_x, arrh_y, 1);
% Extract slope and intercept
slope = p(1);
intercept = p(2);

% Calculate R^2
y_pred = polyval(p, arrh_x);
SS_res = sum((arrh_y - y_pred).^2);  % Sum of squared residuals
SS_tot = sum((arrh_y - mean(arrh_y)).^2);  % Total sum of squares
R_squared = 1 - (SS_res / SS_tot);

% Generate fit line
x_fit = linspace(min(arrh_x), max(arrh_x), 100);
y_fit = polyval(p, x_fit);

% Plot
figure(4); clf;
scatter(arrh_x * 10^3, arrh_y, 'filled');
hold on;
plot(x_fit * 10^3, y_fit, 'r-', 'LineWidth', 1);
xlabel('1/T (1/K)  *10^3');
ylabel('ln(beta * f /s)');
title('Arrhenius Plot');
grid on;
legend('Data', 'Linear Fit');
set(gca, 'YDir','reverse')

% Display the equation
fprintf('FOR %.4f COVERAGE \n', CL_theta);
fprintf('Fit equation: ln(beta * f/s) = %.4f*(1/T) + %.4f\n', slope, intercept); 
fprintf('R^2 = %.4f\n', R_squared);
fprintf('Ea in kJ/mol \n = %.4f\n', slope*(-8.314)/1000); % Actual: 160 +- 20 kJ
fprintf('Pre-exponential factor in s^-1 for 0.8 x 10^19 dosage \n = %.4e\n', exp(intercept));% Actual: 2 x 10^7 s^-1
fprintf('Log10() of pre-exponential \n = %.4e\n', log10(exp(intercept))); % Actual: 7.3 +- 0.5

%OLD CODE, BEFORE FUNCTION SEGREGATION.  6/4
%{ 
%Time trial - with insights from 'erley80test', trying to do Arrhenius
%analysis with actual time integrals 
%Run CLAABGreduce first. 6/3 



%Plotting the background-reduced signals 
n_trace = 8;
cmap = parula(n_trace);

temp_init = 300; %K 
beta = 50; %K/s
time = @(x) (x-temp_init)/beta;

figure(1); clf;
hold on;
plot(time(new_x_0p4), new_y_0p4, 'o-');
plot(time(new_x_0p8), new_y_0p8, 'o-');
plot(time(new_x_1p2), new_y_1p2, 'o-');
plot(time(new_x_1p6), new_y_1p6, 'o-');
plot(time(new_x_2p0), new_y_2p0, 'o-');
plot(time(new_x_2p8), new_y_2p8, 'o-');
plot(time(new_x_4p0), new_y_4p0, 'o-');
plot(time(new_x_8p0), new_y_8p0, 'o-');

set(gca, 'ColorOrder', cmap, 'NextPlot', 'replacechildren');
legend('0.4', '0.8', '1.2', '1.6', '2.0', '2.8', '4.0', '8.0');
xlabel('Time (s)');
ylabel('Spectrometer Signal');
title('Desorption Traces, Signal vs Time')
hold off;

%% Getting outputs of time() of lower coverage desorption spectra

t_0p4 = time(new_x_0p4);
t_0p8 = time(new_x_0p8);
t_1p2 = time(new_x_1p2);
t_1p6 = time(new_x_1p6);



%% 6/4 - Doing what I originally thought was plotted 


%for i = 1:length(signalSPAN_actual_0p8)
yeah = 80;
okay = 122;
for i = yeah:okay
    arrh_x_test(i) = 1/tempSPAN_actual_0p8(i); 
    arrh_y_test(i) = log(beta * signalSPAN_actual_0p8(i) / N_test(i));

end

arrh_x = arrh_x_test(yeah:end);
arrh_y = arrh_y_test(yeah:end);

clear arrh_x_test;
clear arrh_y_test;

%}

%% 6/9 - Cleaning up with new functions. Should all run as one section. Works.
%Scanning raw scanned data
%new_x_0p4_raw = new_fig5_0p4(:,1)'; %ADDED 
%new_y_0p4_raw = new_fig5_0p4(:,2)'; 
new_x_0p8_raw = new_fig5_0p8(:,1)'; 
new_y_0p8_raw = new_fig5_0p8(:,2)'; 
new_x_1p2_raw = new_fig5_1p2(:,1)'; 
new_y_1p2_raw = new_fig5_1p2(:,2)'; 
new_x_1p6_raw = new_fig5_1p6(:,1)'; 
new_y_1p6_raw = new_fig5_1p6(:,2)'; 

%Background reduction
figure(1); clf;
hold on;
%[tempSPAN_actual_0p4 , signalSPAN_actual_0p4] = reduceBGgetspectra(new_x_0p4_raw, new_y_0p4_raw,750,1450,'k');



[tempSPAN_actual_0p8 , signalSPAN_actual_0p8] = reduceBGgetspectra(new_x_0p8_raw, new_y_0p8_raw,750,1450,'g');
[tempSPAN_actual_1p2 , signalSPAN_actual_1p2] = reduceBGgetspectra(new_x_1p2_raw, new_y_1p2_raw,750,1450,'b');
[tempSPAN_actual_1p6 , signalSPAN_actual_1p6] = reduceBGgetspectra(new_x_1p6_raw, new_y_1p6_raw,750,1450,'m');

hold off;

%Initial coverages, from LEED data
N0_0p8 = 0.2004;
N0_1p2 = 0.2877;
N0_1p6 = 0.3856; 

%Getting rates and coverage from the signal
temp_init = 300; %K 
beta = 50; %K/s
time = @(x) (x-temp_init)/beta;
figure(2); clf;
hold on;
[dNdt_0p8, N_0p8] = getcoverageplotcoverage(time(tempSPAN_actual_0p8), signalSPAN_actual_0p8, N0_0p8,'g');
[dNdt_1p2, N_1p2] = getcoverageplotcoverage(time(tempSPAN_actual_1p2), signalSPAN_actual_1p2, N0_1p2,'b');
[dNdt_1p6, N_1p6] = getcoverageplotcoverage(time(tempSPAN_actual_1p6), signalSPAN_actual_1p6, N0_1p6, 'm');
hold off;

%Interpolation scheme to find rates/times 
CL_theta = 0.02; 
[CL_time_0p8, coverage_0p8_theta_0p15] = getpointgetinterpolation_focusYgetX(time(tempSPAN_actual_0p8), N_0p8, CL_theta, 2);
[~, rate_0p8_theta_0p15] = getpointgetinterpolation_focusXgetY(time(tempSPAN_actual_0p8), dNdt_0p8, CL_time_0p8, 2);
[CL_time_1p2, coverage_1p2_theta_0p15] = getpointgetinterpolation_focusYgetX(time(tempSPAN_actual_1p2), N_1p2, CL_theta, 2);
[~, rate_1p2_theta_0p15] = getpointgetinterpolation_focusXgetY(time(tempSPAN_actual_1p2),dNdt_1p2, CL_time_1p2, 2);
[CL_time_1p6, coverage_1p6_theta_0p15] = getpointgetinterpolation_focusYgetX(time(tempSPAN_actual_1p6), N_1p6, CL_theta, 2);
[~, rate_1p6_theta_0p15] = getpointgetinterpolation_focusXgetY(time(tempSPAN_actual_1p6), dNdt_1p6, CL_time_1p6, 2);

%Making points for Arrhenius plot
temp = @(t) beta*t + temp_init;
arrh_x = [1/temp(CL_time_0p8), 1/temp(CL_time_1p2), 1/temp(CL_time_1p6)];
arrh_y = [log(1 * rate_0p8_theta_0p15 / coverage_0p8_theta_0p15), log(1 * rate_1p2_theta_0p15 / coverage_1p2_theta_0p15) ...
    log(1 * rate_1p6_theta_0p15 / coverage_1p6_theta_0p15)];

%Plot and linear regression
p = polyfit(arrh_x, arrh_y, 1);
% Extract slope and intercept
slope = p(1);
intercept = p(2);
% Generate fit line
x_fit = linspace(0, max(arrh_x), 100);
y_fit = polyval(p, x_fit);
% Plot
figure(4); 
cmap = parula(8);
%Error

residuals = arrh_y - (slope * arrh_x + intercept);
residuals_std = std(residuals);
n = length(arrh_x);
arrh_x_mean = mean(arrh_x);
arrh_x_var = sum((arrh_x - arrh_x_mean).^2);
m_std_error = residuals_std/sqrt(arrh_x_var * (n-1));
c_std_error = residuals_std*sqrt(sum(arrh_x.^2)/(n*arrh_x_var));

arrh_y_mean = mean(arrh_y);
ss_total = sum((arrh_y - arrh_y_mean).^2);
ss_residual = sum(residuals.^2);
% Calculate the residuals and standard deviation for the Arrhenius analysis
r_squared = 1 - (ss_residual/ss_total);
r = sqrt(r_squared);
fprintf('Coverage = \n %.2f \n', CL_theta);
fprintf('R_squared = \n %.10f\n', r_squared);
%Actual plot
hold on;
plot(arrh_x * 10^3, arrh_y,'o','LineStyle','none' ,'Color','k');
plot(x_fit * 10^3, y_fit, 'LineWidth', 1, 'Color', 'k' ) ;
xlabel('1/T (1/K)  *10^3');
ylabel('ln(beta * f /s)');
title('Arrhenius Plot');
grid on;
legend('Data', 'Linear Fit');
set(gca, 'YDir','reverse')
hold off;
% Display the equation
fprintf('Fit equation: \n ln(f/s) = %.4f*(1/T) + %.4f\n', slope, intercept); 
fprintf('Ea in kJ/mol \n = %.10f\n', slope*(-8.314)/1000); % Actual: 160 +- 20 kJ
fprintf('Pre-exponential factor in s^-1  \n = %.10e\n', exp(intercept));% Actual: 2 x 10^7 s^-1
fprintf('Log10() of pre-exponential \n = %.10e\n', log10(exp(intercept))); % Actual: 7.3 +- 0.5


%% Making more plots at different coverages to ascertain a linear coverage-activation energy relationship. 