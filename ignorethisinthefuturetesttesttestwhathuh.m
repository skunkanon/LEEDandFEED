%% 6/8 - testing out linear coverage - Ea relationship

coverage_desired = 0.07;
[time_01, coverage_01] = getpointgetinterpolation_focusYgetX(t_Pt_01, c_Pt_01, coverage_desired, 0.0002); %point #1
[time_01_focusXgetY, rate_01] = getpointgetinterpolation_focusXgetY(t_Pt_01,y_Pt_01, time_01, 0.0002);
%6/7- Works.

[time_015, coverage_015] = getpointgetinterpolation_focusYgetX(t_Pt_015, c_Pt_015, coverage_desired, 0.0002); %point #2
[time_015_focusXgetY, rate_015] = getpointgetinterpolation_focusXgetY(t_Pt_015,y_Pt_015, time_015, 0.0002);
[time_02, coverage_02] = getpointgetinterpolation_focusYgetX(t_Pt_02, c_Pt_02, coverage_desired, 0.0002); %point #3
[time_02_focusXgetY, rate_02] = getpointgetinterpolation_focusXgetY(t_Pt_02,y_Pt_02, time_02, 0.0002);
[time_025, coverage_025] = getpointgetinterpolation_focusYgetX(t_Pt_025, c_Pt_025, coverage_desired, 0.0002); %point #4
[time_025_focusXgetY, rate_025] = getpointgetinterpolation_focusXgetY(t_Pt_025,y_Pt_025, time_025, 0.0002);
[time_03, coverage_03] = getpointgetinterpolation_focusYgetX(t_Pt_03, c_Pt_03, coverage_desired, 0.0002); %point #5
[time_03_focusXgetY, rate_03] = getpointgetinterpolation_focusXgetY(t_Pt_03,y_Pt_03, time_03, 0.0002);


%Arrhenius analysis
beta = 8; %K/s 
temp = @(t) beta*t + init_K; 
arrh_x = [1/temp(time_01), 1/temp(time_015), 1/temp(time_02), 1/temp(time_025), 1/temp(time_03)];
beta = 1;
arrh_y = [log(beta * rate_01 / coverage_01), log(beta* rate_015 / coverage_015), log(beta * rate_02 / coverage_02),...
    log(beta * rate_025 / coverage_025), log(beta*rate_03 / coverage_03)];


%%




p = polyfit(arrh_x, arrh_y, 1);
% Extract slope and intercept
slope = p(1);
intercept = p(2);

% Generate fit line
x_fit = linspace(0, max(arrh_x), 100);
y_fit = polyval(p, x_fit);


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



% Plot
figure(4); 
scatter(arrh_x * 10^3, arrh_y, 'filled');
hold on;
plot(x_fit * 10^3, y_fit, 'r-', 'LineWidth', 1);
xlabel('1/T (1/K)  *10^3');
ylabel('ln(beta * f /s)');
title('Arrhenius Plot');
grid on;
legend('Data', 'Linear Fit');
set(gca, 'YDir','reverse')
% Display the equationwoody 
fprintf('Fit equation: ln(beta * f/s) = %.4f*(1/T) + %.4f\n', slope, intercept); 
fprintf('Ea in kcal/mol  \n = %.10f\n', (1/kcal_to_J)*slope*(-8.314)); % Actual: 160 +- 20 kJ
fprintf('Pre-exponential factor in s^-1 \n = %.10e\n', exp(intercept));% Actual: 2 x 10^7 s^-1
fprintf('Log10() of pre-exponential \n = %.10e\n', log10(exp(intercept))); % Actual: 7.3 +- 0.5 
%%

coverages = [c_0p01, c_0p03, c_0p05, c_0p07];
energies = [Ea_0p01, Ea_0p03, Ea_0p05, Ea_0p07];



p = polyfit(coverages, energies, 1);
% Extract slope and intercept
slope = p(1);
intercept = p(2);

% Generate fit line
x_fit = linspace(min(coverages), max(coverages), 100);
y_fit = polyval(p, x_fit);

% Plot
figure(4); clf;
scatter(coverages, energies, 'filled');
hold on;
plot(x_fit, y_fit, 'r-', 'LineWidth', 1);
xlabel('Coverage');
ylabel('Activation Energy (kcal/J)');
title('Arrhenius Plot');
grid on;
legend('Data', 'Linear Fit');
set(gca, 'YDir','reverse')
% Display the equation
fprintf('Fit equation: ln(beta * f/s) = %.4f*(1/T) + %.4f\n', slope, intercept); 
%fprintf('Ea in kcal/mol  \n = %.10f\n', (1/kcal_to_J)*slope*(-8.314)); % Actual: 160 +- 20 kJ
%fprintf('Pre-exponential factor in s^-1 \n = %.10e\n', exp(intercept));% Actual: 2 x 10^7 s^-1
%fprintf('Log10() of pre-exponential \n = %.10e\n', log10(exp(intercept))); % Actual: 7.3 +- 0.5 
%% what am i doing
%{
xdata = ...
 [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = ...
 [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5];

fun = @(x,xdata)x(1)*exp(x(2)*xdata);
x0 = [100,-1];
x = lsqcurvefit(fun,x0,xdata,ydata);

%}
xdata = linspace(0,3);
ydata = exp(-1.3*xdata) + 0.05*randn(size(xdata)); %exponential + noise
lb = [0,-2]; %lower bound
ub = [3/4,-1]; %upper bound
fun = @(x,xdata)x(1)*exp(x(2)*xdata);
x0 = [1/2,-2]; %guess 
x = lsqcurvefit(fun,x0,xdata,ydata,lb,ub);
%%
% Initial guess for parameters [v, Ea_0]
params0 = [1e13 -1000 , 47.5*4184 + 100];  % Adjust these based on your system

% Lower and upper bounds for parameters
lb = [1e10, 40*4184];    % Lower bounds
ub  = [1e15, 60*4184];      % Upper bounds

% Your experimental data
t_data = t_Pt_02 + (2*rand()-1)/100;  % Your time points
rate_data = y_Pt_02() + (2*rand()-1)/100;  % Your experimental rate data

% Fit using lsqcurvefit
[params_fit, resnorm] = lsqcurvefit(@(params,t) polyani_wigner_fit(t,params,init_K,beta), ...
    params0, t_data, rate_data, lb, ub);

% Extract fitted parameters
v_fit = params_fit(1);
Ea_0_fit = params_fit(2);

%% 6/8 - lsqcurvefit on polani_wigner()

kcal_to_J = 4184;
N_01 = 0.1;
N_015 = 0.15;
N_02 = 0.2;
N_025 = 0.25; 
N_03 = 0.3; 
N_035 = 0.35;


init_K = 300; %Kelvin 
pre_exp = 10^13; %s^-1 
Ea_Pt = 47.5; %kcal/mol
gam = 4 *kcal_to_J; %kcal/mol to j/mol 
[t_Pt_01, x_Pt_01, y_Pt_01,c_Pt_01] = polyani_wigner(8,init_K, Ea_Pt * kcal_to_J,gam,pre_exp,N_01,1050); %Choosing 0.1 = N_0 as the starting value, same as Cox '81
[t_Pt_015, x_Pt_015, y_Pt_015,c_Pt_015] = polyani_wigner(8,init_K, Ea_Pt * kcal_to_J,gam,pre_exp,N_015,1050);
[t_Pt_02, x_Pt_02, y_Pt_02,c_Pt_02] = polyani_wigner(8,init_K, Ea_Pt * kcal_to_J,gam,pre_exp,N_02,1050);
[t_Pt_025, x_Pt_025, y_Pt_025,c_Pt_025] = polyani_wigner(8,init_K, Ea_Pt * kcal_to_J,gam,pre_exp,N_025,1050);
[t_Pt_03, x_Pt_03, y_Pt_03,c_Pt_03] = polyani_wigner(8,init_K, Ea_Pt * kcal_to_J,gam,pre_exp,N_03,1050);
[t_Pt_035, x_Pt_035, y_Pt_035,c_Pt_035] = polyani_wigner(8,init_K, Ea_Pt * kcal_to_J,gam,pre_exp,N_035,1050);

fit_polyani_wigner(t_Pt_01, y_Pt_01, [30*4184, 3*4184, 1e12] , 8 , 300, 1050, 0.1);

%{
c_Pt_01 = c_Pt_01 + (2*rand()-1)/1000;
y_Pt_01 = y_Pt_01 + (2*rand()-1)/1000;
c_Pt_015 = c_Pt_015 +(2*rand()-1)/1000;
y_Pt_015 = y_Pt_015 + (2*rand()-1)/1000;
y_Pt_02 = y_Pt_02 + (2*rand()-1)/1000;
c_Pt_02 = c_Pt_02 + (2*rand()-1)/1000;

figure(6); clf(6); 
cmap = parula(6);
hold on;
plot(t_Pt_01, y_Pt_01, t_Pt_01, c_Pt_01, 'Color', cmap(1,:));
plot(t_Pt_015, y_Pt_015,t_Pt_015, c_Pt_015,'Color', cmap(2,:));
plot(t_Pt_02, y_Pt_02, t_Pt_02, c_Pt_02, 'Color', cmap(3,:));
plot(t_Pt_025, y_Pt_025,t_Pt_025, c_Pt_025,'Color', cmap(4,:));
plot(t_Pt_03, y_Pt_03, t_Pt_03, c_Pt_03,'Color', cmap(5,:));
plot(t_Pt_035, y_Pt_035,t_Pt_035, c_Pt_035, 'Color', cmap(6,:));
legend('0.1', '0.15','0.2', '0.25', '0.3', '0.35');
hold off;
 
%}

%%

x_coverage = [0.18, 0.16, 0.14, 0.12, 0.10, 0.08, 0.06, 0.04, 0.02];
y_energy = [177.33, 203.7, 155.917, 117.74, 64.93, 50.97, -9.59 , -72.52, -171.97];

%%

x_coverage = [0.18, 0.16, 0.14, 0.12, 0.10, 0.08, 0.06, 0.04, 0.02];
y_energy = [4.74e7, 6.8e8, 4.69e6, 9.34e4, 4.7e2, 2.09, 3.27e-01, 9.01e-04, 1.12e-07];
y_energy = log10(y_energy);
%%
%Plot and linear regression
p = polyfit(x_coverage, y_energy, 1);
% Extract slope and intercept
slope = p(1);
intercept = p(2);
% Generate fit line
x_fit = linspace(min(x_coverage), max(x_coverage), 100);
y_fit = polyval(p, x_fit);
% Plot
figure(6); clf;
%cmap = parula(8);
%Error

residuals = y_energy - (slope * x_coverage + intercept);
residuals_std = std(residuals);
n = length(x_coverage);
arrh_x_mean = mean(x_coverage);
arrh_x_var = sum((x_coverage - arrh_x_mean).^2);
m_std_error = residuals_std/sqrt(arrh_x_var * (n-1));
c_std_error = residuals_std*sqrt(sum(x_coverage.^2)/(n*arrh_x_var));
fprintf('Standard error in intercept = \n %.5f (kJ/mol) \n', c_std_error);
arrh_y_mean = mean(y_energy);
ss_total = sum((y_energy - arrh_y_mean).^2);
ss_residual = sum(residuals.^2);
% Calculate the residuals and standard deviation for the Arrhenius analysis
r_squared = 1 - (ss_residual/ss_total);
r = sqrt(r_squared);
%fprintf('Coverage = \n %.2f \n', CL_theta);
fprintf('R_squared = \n %.10f\n', r_squared);
%Actual plot
hold on;
plot(x_coverage, y_energy,'o','LineStyle','none' ,'Color','k');
plot(x_fit, y_fit, 'LineWidth', 1, 'Color', 'k' ) ;
xlabel('Coverage');
ylabel('Activation Energy kJ/mol');
%title('Arrhenius Plot');
grid on;
legend('Data', 'Linear Fit');
%set(gca, 'YDir','reverse')
hold off;
% Display the equation
fprintf('Fit equation: \n  = %.2f(kJ/mol)*Coverage + %.2f(kJ/mol)\n', slope, intercept); 
fprintf('Ea_0 in kJ/mol \n = %.2f ± %.2f\n', intercept, c_std_error); 
fprintf('Ea (0.2 coverage) in kJ/mol \n = %.2f ± %.2f\n', intercept + slope*0.2, c_std_error); % Actual: 160 +- 20 kJ
%fprintf('Pre-exponential factor in s^-1  \n = %.10e\n', exp(intercept));% Actual: 2 x 10^7 s^-1
%fprintf('Log10() of pre-exponential \n = %.10e\n', log10(exp(intercept))); % Actual: 7.3 +- 0.5


%%

residuals = y_energy - (slope * x_coverage + intercept);
residuals_std = std(residuals);
n = length(x_coverage);
arrh_x_mean = mean(x_coverage);
arrh_x_var = sum((x_coverage - arrh_x_mean).^2);
m_std_error = residuals_std/sqrt(arrh_x_var * (n-1));
c_std_error = residuals_std*sqrt(sum(x_coverage.^2)/(n*arrh_x_var));
fprintf('Standard error in intercept = \n %.5f \n', c_std_error);
arrh_y_mean = mean(y_energy);
ss_total = sum((y_energy - arrh_y_mean).^2);
ss_residual = sum(residuals.^2);
% Calculate the residuals and standard deviation for the Arrhenius analysis
r_squared = 1 - (ss_residual/ss_total);
r = sqrt(r_squared);
%fprintf('Coverage = \n %.2f \n', CL_theta);
fprintf('R_squared = \n %.10f\n', r_squared);
%Actual plot
figure(10); clf;
hold on;
plot(x_coverage, y_energy,'o','LineStyle','none' ,'Color','k');
plot(x_fit, y_fit, 'LineWidth', 1, 'Color', 'k' ) ;
xlabel('Coverage');
ylabel('Activation Energy kJ/mol');
%title('Arrhenius Plot');
grid on;
legend('Data', 'Linear Fit');
%set(gca, 'YDir','reverse')
hold off;
% Display the equation
fprintf('Fit equation: \n  log10(v) = %.2f*Coverage + %.2f\n', slope, intercept); 
fprintf('log10(v) at 0 coverage \n = %.2f ± %.2f\n', intercept, c_std_error); 
fprintf('log10(v) at 0.2 coverage \n = %.2f ± %.2f\n', intercept + slope*0.2, c_std_error); % Actual: 160 +- 20 kJ
%fprintf('Pre-exponential factor in s^-1  \n = %.10e\n', exp(intercept));% Actual: 2 x 10^7 s^-1
%fprintf('Log10() of pre-exponential \n = %.10e\n', log10(exp(intercept))); % Actual: 7.3 +- 0.5


%%
kcal_to_J = 4184;
Ea_Pt = 47.5; 
init_K = 300;
pre_exp = 10^13;
w = 5 ; %kJ/mol
gam=4 * kcal_to_J;
[t_Pt_01, x_Pt_01, y_Pt_01,c_Pt_01] = polyani_wigner_comp1(8,init_K, Ea_Pt * kcal_to_J,gam,w,pre_exp,N_01,1050); %Choosing 0.1 = N_0 as the starting value, same as Cox '81
[t_Pt_015, x_Pt_015, y_Pt_015,c_Pt_015] = polyani_wigner_comp1(8,init_K, Ea_Pt * kcal_to_J,gam,w,pre_exp,N_015,1050);
[t_Pt_02, x_Pt_02, y_Pt_02,c_Pt_02] = polyani_wigner_comp1(8,init_K, Ea_Pt * kcal_to_J,gam,w,pre_exp,N_02,1050);
[t_Pt_025, x_Pt_025, y_Pt_025,c_Pt_025] = polyani_wigner_comp1(8,init_K, Ea_Pt * kcal_to_J,gam,w,pre_exp,N_025,1050);
[t_Pt_03, x_Pt_03, y_Pt_03,c_Pt_03] = polyani_wigner_comp1(8,init_K, Ea_Pt * kcal_to_J,gam,w,pre_exp,N_03,1050);
[t_Pt_035, x_Pt_035, y_Pt_035,c_Pt_035] = polyani_wigner_comp1(8,init_K, Ea_Pt * kcal_to_J,gam, w,pre_exp,N_035,1050);

%{
c_Pt_01 = c_Pt_01 + (2*rand()-1)/1000;
y_Pt_01 = y_Pt_01 + (2*rand()-1)/1000;
c_Pt_015 = c_Pt_015 +(2*rand()-1)/1000;
y_Pt_015 = y_Pt_015 + (2*rand()-1)/1000;
y_Pt_02 = y_Pt_02 + (2*rand()-1)/1000;
c_Pt_02 = c_Pt_02 + (2*rand()-1)/1000;
%}
figure(6); clf(6); 
cmap = parula(6);
hold on;
plot(t_Pt_01, y_Pt_01, t_Pt_01, c_Pt_01, 'Color', cmap(1,:));
plot(t_Pt_015, y_Pt_015,t_Pt_015, c_Pt_015,'Color', cmap(2,:));
plot(t_Pt_02, y_Pt_02, t_Pt_02, c_Pt_02, 'Color', cmap(3,:));
plot(t_Pt_025, y_Pt_025,t_Pt_025, c_Pt_025,'Color', cmap(4,:));
plot(t_Pt_03, y_Pt_03, t_Pt_03, c_Pt_03,'Color', cmap(5,:));
plot(t_Pt_035, y_Pt_035,t_Pt_035, c_Pt_035, 'Color', cmap(6,:));
legend('0.1', '0.15','0.2', '0.25', '0.3', '0.35');
hold off;


%% 6/11 - fit_polyani_wigner_comp1 

init_params = [170 * 1000, 0 , 0 , 10^7.4]; %Ea, gam, w, pre-exponent 
figure(3); clf;
temp_init = 300;
beta = 50; 
time = @(x) (x-temp_init)/beta;
hold on;

dNdt_0p8(1:18) = 0;
dNdt_0p8(98:end) = 0;

plot(time(tempSPAN_actual_0p8) ,dNdt_0p8, 'r' );
[timesim_0p8, ~, rate_sim_0p8, ~] = polyani_wigner_comp1(50, 300,init_params(1), init_params(2), init_params(3), init_params(4), N0_0p8, 1500);
plot(timesim_0p8, rate_sim_0p8, 'b');
hold off;



%%
fit_polyani_wigner_comp1(time(tempSPAN_actual_0p8), dNdt_0p8, init_params, 50, 300,1500, N0_0p8);
%% 6/11 - fit_polyani_wigner_comp1 

init_params = [170 * 1000, 0 , 0 , 10^7.4]; %Ea, gam, w, pre-exponent 
figure(3); clf;
temp_init = 300;
beta = 50; 
time = @(x) (x-temp_init)/beta;
hold on;

%dNdt_0p8(1:18) = 0;
%dNdt_0p8(98:end) = 0;

plot(time(tempSPAN_actual_1p6) ,dNdt_1p6, 'r' );
[timesim_1p6, ~, rate_sim_1p6, ~] = polyani_wigner_comp1(50, 300,init_params(1), init_params(2), init_params(3), init_params(4), N0_1p6, 1500);
plot(timesim_1p6, rate_sim_1p6, 'b');
hold off;



%%
fit_polyani_wigner_comp1(time(tempSPAN_actual_1p6), dNdt_1p6, init_params, 50, 300,1500, N0_1p6);
%% 6/11 - fit_polyani_wigner_comp1 

init_params = [43.07214353274571 * 4184, -6 * 1000 , 1.0  , 4.35349687969534e+07]; %Ea, gam, w, pre-exponent 
figure(3); clf;
temp_init = 300;
beta = 50; 
time = @(x) (x-temp_init)/beta;
hold on;

%dNdt_1p2(1:18) = 0;
%dNdt_0p8(98:end) = 0;

plot(time(tempSPAN_actual_1p2) ,dNdt_1p2, 'r' );
[timesim_1p2, ~, rate_sim_1p2, ~] = polyani_wigner_comp1(50, 300,init_params(1), init_params(2), init_params(3), init_params(4), N0_1p2, 1500);
plot(timesim_1p2, rate_sim_1p2, 'b');
hold off;


%%
fit_polyani_wigner_comp1(time(tempSPAN_actual_1p2), dNdt_1p2, init_params, 50, 300,1500, N0_1p2);


%%
%6/16- Testing out 'polyani_wigner_niemant', fig4

Ea_0 = 300 * 1000; %J -> kJ/mol
v_0 = 10^13; 
beta = 1; %K/s
w = 25 * 1000; %J -> kJ/mol
theta_0 = linspace(0.1, 1, 10);
Tc = [600, 1000, 2000, 10^20]; 

%{
[time_test, tmp_test, rate_test, cov_test] = polyani_wigner_niemant(beta, 300, Ea_0, w, v_0, theta_0(1),1200, Tc(1));
subplot(2,2,1); cla;
%figure(1); clf(1);
hold on;
for i=1:length(theta_0)
    [time, ~, rate, ~] = polyani_wigner_niemant(beta, 300, Ea_0,w,v_0, theta_0(i), 1200, Tc(1));
    plot(time, rate, 'DisplayName', sprintf('Coverage = %.2f', theta_0(i)));
end
ylim([0 0.025]);
title(sprintf('w = %.1f, T_c = %.1f K', w, Tc(1)))
hold off;

subplot(2,2,2); cla;
%clf;
%figure(2); clf(2);
hold on;
for i=1:length(theta_0)
    [time, ~, rate, ~] = polyani_wigner_niemant(beta, 300, Ea_0,w,v_0, theta_0(i), 1200, Tc(2));
    plot(time, rate, 'DisplayName', sprintf('Coverage = %.2f', theta_0(i)));
end
ylim([0 0.025]);
title(sprintf('w = %.1f, T_c = %.1f K', w, Tc(2)))
hold off;

%}


%%
figure(2); clf;
for i=1:length(Tc)
    % Left column (original figure 2)
    subplot(4,2, 2*i-1); cla;  % This will give indices 1,3,5,7
    colors = parula(length(theta_0));
    for j=1:length(theta_0)
        hold on;
        [time, ~, rate, ~] = polyani_wigner_niemant(beta, 300, Ea_0,w,v_0, theta_0(j), 1200, Tc(i));
        plot(time, rate, 'Color', colors(j,:),'DisplayName', sprintf('Coverage = %.2f', theta_0(j)));
        ylim([0 0.025]);
        title(sprintf('w = %.1f kJ/mol, T_c = %.1f K', w/1000, Tc(i)));
        legend('show', 'Location', 'best');
        hold off;
    end
end

w = -w; 

% Right column (original figure 3)
for i=1:length(Tc)
    subplot(4,2, 2*i); cla;  % This will give indices 2,4,6,8
    colors = parula(length(theta_0));
    for j=1:length(theta_0)
        hold on;
        [time, ~, rate, ~] = polyani_wigner_niemant(beta, 300, Ea_0,w,v_0, theta_0(j), 1200, Tc(i));
        plot(time, rate, 'Color', colors(j,:),'DisplayName', sprintf('Coverage = %.2f', theta_0(j)));
        ylim([0 0.025]);
        title(sprintf('w = %.1f kJ/mol, T_c = %.1f K', w/1000, Tc(i)));
        legend('show', 'Location', 'best');
        hold off;
    end
end

sgtitle('Figure 4, Niemantsverdriet. Attractive (left) and repulsive (right) pairwise interactions', 'FontSize', 14)


%%

gfftrgft