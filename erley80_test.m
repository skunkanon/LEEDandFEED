% 3/19 - Figure 3, Chlorine-saturated Pd(111) and Pt(111) Desorption 
fprintf('############ NEW INSTANCE ############ \n');

Ea_Pd = 60.5; %kcal/mol
Ea_Pt = 47.5; %kcal/mol
kcal_to_J = 4184;

[t_Pd,x_Pd,y_Pd] = polyani_wigner(8,550, Ea_Pd * kcal_to_J,0,10^13,1,1050);
[t_Pt,x_Pt,y_Pt] = polyani_wigner(8,550, Ea_Pt * kcal_to_J,0,10^13,1,1050);

%%
fprintf('\n');

%% 3/19 - Making Plots

figure;
hold on;
plot(x_Pd, y_Pd);
plot(x_Pt,y_Pt);

[M, I] = max(y_Pd);
plot(x_Pd(I), M, 'o', 'MarkerSize',10);





hold off;

%% 5/30 - Making multiple desorption curves at different coverages
N_01 = 0.1;
N_015 = 0.15;
N_02 = 0.2;
N_025 = 0.25; 
N_03 = 0.3; 
N_035 = 0.35;

[t_Pt_01, x_Pt_01, y_Pt_01,~] = polyani_wigner(8,550, Ea_Pt * kcal_to_J,0,10^13,N_01,1050); %Choosing 0.1 = N_0 as the starting value, same as Cox '81
[t_Pt_015, x_Pt_015, y_Pt_015,~] = polyani_wigner(8,550, Ea_Pt * kcal_to_J,0,10^13,N_015,1050);
[t_Pt_02, x_Pt_02, y_Pt_02,~] = polyani_wigner(8,550, Ea_Pt * kcal_to_J,0,10^13,N_02,1050);
[t_Pt_025, x_Pt_025, y_Pt_025,~] = polyani_wigner(8,550, Ea_Pt * kcal_to_J,0,10^13,N_025,1050);
[t_Pt_03, x_Pt_03, y_Pt_03,~] = polyani_wigner(8,550, Ea_Pt * kcal_to_J,0,10^13,N_03,1050);
[t_Pt_035, x_Pt_035, y_Pt_035,~] = polyani_wigner(8,550, Ea_Pt * kcal_to_J,0,10^13,N_035,1050);

figure(1); clf(1); 
cmap = parula(6);
hold on;
plot(t_Pt_01, y_Pt_01, 'Color', cmap(1,:));
plot(t_Pt_015, y_Pt_015, 'Color', cmap(2,:));
plot(t_Pt_02, y_Pt_02, 'Color', cmap(3,:));
plot(t_Pt_025, y_Pt_025, 'Color', cmap(4,:));
plot(t_Pt_03, y_Pt_03, 'Color', cmap(5,:));
plot(t_Pt_035, y_Pt_035, 'Color', cmap(6,:));
legend('0.1', '0.15','0.2', '0.25', '0.3', '0.35');
hold off;

%% 6/1 - More trapz investigation. 

x_test = [0 ,0.1, 0.2, 0.3];
y_test = [2, 2, 2, 2]; %Area should be  (0.3 * 2 = 0.6.)

%trapz(x_test, y_test returns 0.6.) 
%trapz(y_test) returns 6, assuming spacing of 1. 
%trapz(x_Pt_01, y_Pt_01) returns.. 0.8. Not the initial coverage. 


x_test2 = [0 ,0.15, 0.2, 0.3]; %Uneven spacing. 
y_test2 = [2, 2, 2, 2]; %Area should be  (0.3 * 2 = 0.6.)
%trapz(x,y) still returns 0.6. 
%% 6/1 - Plotting the derivative of the coverage curve, to see if it matches with the generated rate curve from 'polyani_wigner' 
figure(3); clf;
hold on;
%plot(x_Pt_01, Ns_01, 'o-'); 
plot(x_Pt_01, y_Pt_01);
plot(x_Pt_01, c_Pt_01);
hold off;
%6/1 (later in the day) - It does. 

%% 5/31 - In this ideal state, making sure simple integration works. Reading off of 'desorption_feb21.m' and coxlambert code. 
% 6/1 - Should simplify for-loop in CL81AA. 
fprintf('NEW INSTance \n');
Ns_01 = zeros(1,length(y_Pt_01));

for i = 1:length(y_Pt_01)
    Ns_01(i) = trapz(y_Pt_01(i:end)); %gets the area underneath the curve past the signal (new_y) at index 'i'
end

[max_Ns_01, ~] = max(Ns_01); %gets maximum of areas for below normalization 
Ns_01 = Ns_01 * (N_01 / max_Ns_01); %trapz assumes a spacing of 1, so if the actual spacing is lower the reported area's going to be higher
 %more on that, it makes the maximum for new_N_ (the first value)be equal to the initial coverage, N0_

fprintf('NEW INSTANCE \n')
figure(2); clf; 
hold on;
%legend('0.4', '0.8', '1.2', '1.6', '2.0', '2.8', '4.0', '8.0');

plot(x_Pt_01 , Ns_01, 'Color', cmap(1,:), 'DisplayName', '0.1 Normalized Converage');
plot(x_Pt_01, y_Pt_01 ,'--','Color', cmap(1,:), 'DisplayName','0.1 Normalized Signal');

set(gca, 'ColorOrder', cmap, 'NextPlot', 'replacechildren');

hold off;


%


Ns_02 = zeros(1,length(y_Pt_02));

for i = 1:length(y_Pt_02)
    Ns_02(i) = trapz(y_Pt_02(i:end)); %gets the area underneath the curve past the signal (new_y) at index 'i'
end

[max_Ns_02, ~] = max(Ns_02); %gets maximum of areas for below normalization 
Ns_02 = Ns_02 * (N_02 / max_Ns_02); %trapz assumes a spacing of 1, so if the actual spacing is lower the reported area's going to be higher
 %more on that, it makes the maximum for new_N_ (the first value)be equal to the initial coverage, N0_

fprintf('NEW INSTANCE \n')
figure(2); 
hold on;
%legend('0.4', '0.8', '1.2', '1.6', '2.0', '2.8', '4.0', '8.0');

plot(x_Pt_02 , Ns_02, 'Color', cmap(1,:), 'DisplayName', '0.1 Normalized Converage');
plot(x_Pt_02, y_Pt_02 ,'--','Color', cmap(1,:), 'DisplayName','0.1 Normalized Signal');

set(gca, 'ColorOrder', cmap, 'NextPlot', 'replacechildren');

hold off;

%

Ns_03 = zeros(1,length(y_Pt_03));

for i = 1:length(y_Pt_03)
    Ns_03(i) = trapz(y_Pt_03(i:end)); %gets the area underneath the curve past the signal (new_y) at index 'i'
end

[max_Ns_03, ~] = max(Ns_03); %gets maximum of areas for below normalization 
Ns_03 = Ns_03 * (N_03 / max_Ns_03); %trapz assumes a spacing of 1, so if the actual spacing is lower the reported area's going to be higher
 %more on that, it makes the maximum for new_N_ (the first value)be equal to the initial coverage, N0_

fprintf('NEW INSTANCE \n')
figure(2); 
hold on;
%legend('0.4', '0.8', '1.2', '1.6', '2.0', '2.8', '4.0', '8.0');

plot(x_Pt_03 , Ns_03, 'Color', cmap(1,:), 'DisplayName', '0.1 Normalized Converage');
plot(x_Pt_03, y_Pt_03 ,'--','Color', cmap(1,:), 'DisplayName','0.1 Normalized Signal');

set(gca, 'ColorOrder', cmap, 'NextPlot', 'replacechildren');

hold off;

%% 6/1 - Arrhenius analysis on these three points for 0.06 coverage. 


arrh_x = [1/755.205, 1/755.225, 1/782.232];
arrh_y = [log(0.02166/0.06),log(0.0108/0.06), log(0.032/0.06)];
%And they're not linear at all. 

%%
p = polyfit(arrh_x, arrh_y, 1);

% Extract slope and intercept
slope = p(1);
intercept = p(2);

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
fprintf('Fit equation: ln(beta * f/s) = %.4f*(1/T) + %.4f\n', slope, intercept); 
fprintf('Ea in kJ/mol for 0.8 x 10^19 dosage \n = %.4f\n', slope*(-8.314)/1000); % Actual: 160 +- 20 kJ
fprintf('Pre-exponential factor in s^-1 for 0.8 x 10^19 dosage \n = %.4e\n', exp(intercept));% Actual: 2 x 10^7 s^-1
fprintf('Log10() of pre-exponential \n = %.4e\n', log10(exp(intercept))); % Actual: 7.3 +- 0.5 



%% 6/1 (PM) - The integral WRT time of the raw simulated signal (the desorption rate) is the initial coverage. The one WRT temperature isn't. 
%The integral WRT temperature is exactly 8 times the actual coverage, with
%a heating rate of 8 K/s (0.8). Raising that to 35 K/s, it's 3.5, 35 times
%the actual coverage. 
%Let's do Arrhenius analysis from the raw differential equation. 

%Last parameter on the polanyi_wigner() output is the coverage, from
%solving the Madix equation. 'y' is its derivative. 
init_K = 300; %Kelvin 
pre_exp = 10^13; %s^-1 
Ea_Pt = 47.5; %kcal/mol

[t_Pt_01, x_Pt_01, y_Pt_01,c_Pt_01] = polyani_wigner(8,init_K, Ea_Pt * kcal_to_J,0,pre_exp,N_01,1050); %Choosing 0.1 = N_0 as the starting value, same as Cox '81
[t_Pt_015, x_Pt_015, y_Pt_015,c_Pt_015] = polyani_wigner(8,init_K, Ea_Pt * kcal_to_J,0,pre_exp,N_015,1050);
[t_Pt_02, x_Pt_02, y_Pt_02,c_Pt_02] = polyani_wigner(8,init_K, Ea_Pt * kcal_to_J,0,pre_exp,N_02,1050);
[t_Pt_025, x_Pt_025, y_Pt_025,c_Pt_025] = polyani_wigner(8,init_K, Ea_Pt * kcal_to_J,0,pre_exp,N_025,1050);
[t_Pt_03, x_Pt_03, y_Pt_03,c_Pt_03] = polyani_wigner(8,init_K, Ea_Pt * kcal_to_J,0,pre_exp,N_03,1050);
[t_Pt_035, x_Pt_035, y_Pt_035,c_Pt_035] = polyani_wigner(8,init_K, Ea_Pt * kcal_to_J,0,pre_exp,N_035,1050);

%{
c_Pt_01 = c_Pt_01 + (2*rand()-1)/1000;
y_Pt_01 = y_Pt_01 + (2*rand()-1)/1000;
c_Pt_015 = c_Pt_015 +(2*rand()-1)/1000;
y_Pt_015 = y_Pt_015 + (2*rand()-1)/1000;
y_Pt_02 = y_Pt_02 + (2*rand()-1)/1000;
c_Pt_02 = c_Pt_02 + (2*rand()-1)/1000;
%}
figure(6); 
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
 
%% 5/31 Transforming the signals to make them 'realistic', as in making their integrals proportional (not equal) to the initial coverage

signalboot = 1.4;
s_Pt_01 = y_Pt_01 * signalboot; 
s_Pt_015 = y_Pt_015 * signalboot; 
s_Pt_02 = y_Pt_02 * signalboot; 
s_Pt_025 = y_Pt_025 * signalboot; 
s_Pt_030 = y_Pt_03 * signalboot; 
%% 5/31 - Now integrating with code from 'coxlambertactualACTUAL' 
%5/29 - Signals must be normalized to 1; divide all of the signal by the
%max value of the top initial coverage, in [the previous, Cox and Lambert] case the maximum of
%new_y_8p0.
%6/2 Or not. 
%[SIM_max_new_y_all , ~] = max(s_Pt_030); 
fprintf('NEW INSTance \n');
SB_Ns_01 = zeros(1,length(s_Pt_01));

for i = 1:length(s_Pt_01) -1 
    SB_Ns_01(i) = trapz(t_Pt_01(i:end) , s_Pt_01(i:end)); %gets the area underneath the curve past the signal (new_y) at index 'i'
end


fprintf('NEW INSTance \n');

[max_Ns_01, ~] = max(SB_Ns_01); %gets maximum of areas for below normalization 
normed_SB_Ns_01 = SB_Ns_01 * (N_01 / max_Ns_01); %trapz assumes a spacing of 1, so if the actual spacing is lower the reported area's going to be higher
 %more on that, it makes the maximum for new_N_ (the first value)be equal to the initial coverage, N0_
normed_s_Pt_01 = s_Pt_01 * (N_01/max_Ns_01);
 %Works. 

figure(5); 
hold on;
plot(t_Pt_01, y_Pt_01, 'k--','LineWidth',5); %Raw rate
plot(t_Pt_01, c_Pt_01, 'k--', 'LineWidth',5); %Raw rate coverage 
plot(t_Pt_01, s_Pt_01,'r', t_Pt_01, SB_Ns_01, 'r'); %Transformed raw rate and coverage
plot(t_Pt_01, normed_s_Pt_01, 'm--', t_Pt_01, normed_SB_Ns_01, 'm--'); %Normalized transformed

hold off;

%Works.



%% 6/7 - Test out new interpolation functions (focusXgetY, focusYgetX) on simulated data. May make the Arrhenius analysis its own function. 

coverage_desired = 0.05;
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





p = polyfit(arrh_x, arrh_y, 1);
% Extract slope and intercept
slope = p(1);
intercept = p(2);

% Generate fit line
x_fit = linspace(min(arrh_x), max(arrh_x), 100);
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
fprintf('Fit equation: ln(beta * f/s) = %.4f*(1/T) + %.4f\n', slope, intercept); 
fprintf('Ea in kcal/mol  \n = %.10f\n', (1/kcal_to_J)*slope*(-8.314)); % Actual: 160 +- 20 kJ
fprintf('Pre-exponential factor in s^-1 \n = %.10e\n', exp(intercept));% Actual: 2 x 10^7 s^-1
fprintf('Log10() of pre-exponential \n = %.10e\n', log10(exp(intercept))); % Actual: 7.3 +- 0.5 

%% 6/8 - Now testing a linear energy-coverage dependence. 


init_K = 300; %Kelvin 
pre_exp = 10^13; %s^-1 
kcal_to_J = 4184; 
Ea_Pt = 47.5; %kcal/mol
gam = 4 *kcal_to_J; %kcal/mol to j/mol 
[t_Pt_01, x_Pt_01, y_Pt_01,c_Pt_01] = polyani_wigner(8,init_K, Ea_Pt * kcal_to_J,gam,pre_exp,N_01,1050); %Choosing 0.1 = N_0 as the starting value, same as Cox '81
[t_Pt_015, x_Pt_015, y_Pt_015,c_Pt_015] = polyani_wigner(8,init_K, Ea_Pt * kcal_to_J,gam,pre_exp,N_015,1050);
[t_Pt_02, x_Pt_02, y_Pt_02,c_Pt_02] = polyani_wigner(8,init_K, Ea_Pt * kcal_to_J,gam,pre_exp,N_02,1050);
[t_Pt_025, x_Pt_025, y_Pt_025,c_Pt_025] = polyani_wigner(8,init_K, Ea_Pt * kcal_to_J,gam,pre_exp,N_025,1050);
[t_Pt_03, x_Pt_03, y_Pt_03,c_Pt_03] = polyani_wigner(8,init_K, Ea_Pt * kcal_to_J,gam,pre_exp,N_03,1050);
[t_Pt_035, x_Pt_035, y_Pt_035,c_Pt_035] = polyani_wigner(8,init_K, Ea_Pt * kcal_to_J,gam,pre_exp,N_035,1050);

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
 





%% 3/19 - Make function that plots desorption and labels them nicely? Making desorption spectra 




% as cells, not linear arrays
x_Erley = cell(1,2);
y_Erley = cell(1,2);

[~,x_Erley{1}, y_Erley{1}] = polyani_wigner(8,550, Ea_Pd * kcal_to_J,0,10^13,0.03,1050);
[~,x_Erley{2}, y_Erley{2}] = polyani_wigner(8,550, Ea_Pt * kcal_to_J,0,10^13,0.01,1050);




%% Testing out the plot() function 


labels = ["Saturated Pd (111)", "Saturated Pt (111)"];
tight = 'Cl Desorption at 8 K/s';


polyani_plot(x_Erley,y_Erley, labels, tight);

%{
%% 6/2 - Run this through above section 
%0.09 Coverage 
beta = 8; %K/s 


temp = @(t) beta*t + init_K; 


arrh_x = [1/temp(21.2), 1/temp(25.65), 1/temp(26.9645),1/temp(27.6527), 1/temp(28.1532),1/temp(28.52)];
arrh_y = [log(beta* 0.0034/0.09),log(beta* 0.0162/0.09), log(beta*0.0249/0.09),log(beta* 0.0314/0.09), log(beta * 0.0366/0.09), log(beta* 0.0409/0.09)];
 
% It works. 
%%
%0.05 Coverage





arrh_x = [1/temp(26.5265), 1/temp(27.9029), 1/temp(28.5911),1/temp(29.029), 1/temp(29.3418),1/temp(29.5921)];
arrh_y = [log(beta* 0.012/0.05),log(beta* 0.0186/0.05), log(beta*0.0232/0.05),log(beta* 0.0267/0.05), log(beta * 0.0296/0.05), log(beta* 0.03209/0.05)];

%% 6/3 - 0.05 Coverage to the 6th decimal point 

arrh_x = [ 1/temp(26.55649) , 1/temp(27.89285211)];
arrh_y = [log( beta * 0.012106936/ 0.05), log(beta * 0.018677 / 0.05 )];

%Gets energy but still overshoots pre-exponential significantly, by a
%factor of 8, although it's a bit closer. 
%% 6/3 - Reducing initial temperature, also 0.05 coverage 
%Remember to change temp() parameters back to 550 if needed 

beta = 8 ;
arrh_x = [ 1/temp(57.80) , 1/temp(59.1216), 1/temp(59.83), 1/temp(60.27)];
arrh_y = [log( beta * 0.0121/ 0.05), log(beta * 0.018 / 0.05 ), log(beta * 0.0232/ 0.05) ,log(beta * 0.026799/ 0.05)];
%Nothing changes.

%% 6/3 - Raising pre-exponential to 10^14 to see if Arrhenius plot reflects that change well 


arrh_x = [ 1/temp(51.611) , 1/temp(52.7676)];
arrh_y = [log( beta * 0.0137/ 0.05), log(beta * 0.0213 / 0.05 )];
%% 6/3 - Now 10^15


arrh_x = [ 1/temp(46.17) , 1/temp(47.2032)];
arrh_y = [log( beta * 0.01556/ 0.05), log(beta * 0.02416 / 0.05 )];

%% 6/2 Testing out integration function. 
Ns_01 = zeros(1,length(y_Pt_01));
for i = 1:length(y_Pt_01)-1 %necessary to subtract by one for some reason but max() returns the initial coverage
    Ns_01(i) = trapz(t_Pt_01(i:end),y_Pt_01(i:end));

end
%figure;
%plot(t_Pt_01, y_Pt_01,'k', t_Pt_01, Ns_01,'o-'); Works. 


%% 6/7 - Interpolation scheme test

x_span = t_Pt_01(:).';
y_span = c_Pt_01(:).';
desired_y = 0.05;

[~,ind_center] = min(abs(y_span - desired_y));
tolerance = 0.5;
avg_spacing = abs(x_span(1) - x_span(end)) / size(x_span,2);
ind_range = round(tolerance/avg_spacing) / 2;
x_span_focus = x_span(ind_center - ind_range : ind_center + ind_range);
y_span_focus =  y_span(ind_center - ind_range :  ind_center + ind_range);


figure(7); clf;
hold on;
plot(x_span(ind_center - ind_range : ind_center + ind_range) , y_span(ind_center - ind_range :  ind_center + ind_range), 'b--','LineWidth',3);
hold off;


x_span_focus_interp = linspace(min(x_span_focus), max(x_span_focus),2000);
y_span_focus_interp = interp1(x_span_focus,y_span_focus, x_span_focus_interp, 'pchip');

figure(7);
hold on;
plot(x_span_focus_interp, y_span_focus_interp, 'r','LineWidth',2);
plot(x_span, y_span, 'k--')
hold off;

[~,ind_actual] = min(abs(y_span_focus_interp -desired_y));
x_actual = x_span_focus_interp(ind_actual);
y_actual = y_span_focus_interp(ind_actual);

%6/7 - Works. 
%}


%whats up cursor

%{
%% Using new functions to get more accurate values for time, signal, and coverage. 6/9. 

CL_theta = 0.2; 
temp_init = 300; %K 
beta = 50; %K/s
time = @(x) (x-temp_init)/beta;
%[CL_time_0p8, coverage_0p8_theta_0p15] = getpointgetinterpolation_focusYgetX(time(tempSPAN_actual_0p8), N_0p8, CL_theta, 0.0002);
%[~, rate_0p8_theta_0p15] = getpointgetinterpolation_focusXgetY(time(tempSPAN_actual_0p8), signalSPAN_actual_0p8, CL_time_0p8, 0.0002);

[CL_time_1p2, coverage_1p2_theta_0p15] = getpointgetinterpolation_focusYgetX(time(tempSPAN_actual_1p2), N_1p2, CL_theta, 0.0002);
[~, rate_1p2_theta_0p15] = getpointgetinterpolation_focusXgetY(time(tempSPAN_actual_1p2), signalSPAN_actual_1p2, CL_time_1p2, 0.0002);


[CL_time_1p6, coverage_1p6_theta_0p15] = getpointgetinterpolation_focusYgetX(time(tempSPAN_actual_1p6), N_1p6, CL_theta, 0.0002);
[~, rate_1p6_theta_0p15] = getpointgetinterpolation_focusXgetY(time(tempSPAN_actual_1p6), signalSPAN_actual_1p6, CL_time_1p6, 0.0002);

%{
%% Arrhenius analysis. Still have to fill in the values manually. 6/4

temp = @(t) beta*t + temp_init;
arrh_x = [1/temp(CL_time_0p8), 1/temp(CL_time_1p2), 1/temp(CL_time_1p6)];
%beta = 1;
arrh_y = [log(1 * rate_0p8_theta_0p15 / coverage_0p8_theta_0p15), log(1 * rate_1p2_theta_0p15 / coverage_1p2_theta_0p15) ...
    log(1 * rate_1p6_theta_0p15 / coverage_1p6_theta_0p15)];

%%
%}


arrh_x = [ 1/temp(CL_time_1p2), 1/temp(CL_time_1p6)];
%beta = 1;
arrh_y = [log(1 * rate_1p2_theta_0p15 / coverage_1p2_theta_0p15) ...
    log(1 * rate_1p6_theta_0p15 / coverage_1p6_theta_0p15)];

p = polyfit(arrh_x, arrh_y, 1);
% Extract slope and intercept
slope = p(1);
intercept = p(2);
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
fprintf('Fit equation: ln(beta * f/s) = %.4f*(1/T) + %.4f\n', slope, intercept); 
fprintf('Ea in kJ/mol for 0.8 x 10^19 dosage \n = %.4f\n', slope*(-8.314)/1000); % Actual: 160 +- 20 kJ
fprintf('Pre-exponential factor in s^-1 for 0.8 x 10^19 dosage \n = %.4e\n', exp(intercept));% Actual: 2 x 10^7 s^-1
fprintf('Log10() of pre-exponential \n = %.4e\n', log10(exp(intercept))); % Actual: 7.3 +- 0.5

%}