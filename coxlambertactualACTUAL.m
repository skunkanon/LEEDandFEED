%5/21 - Splitting up coxlambertACTUAL, redoing scans, making it more
%organized. 
%eventually make theoretical curves with calculated data (Ea and
%pre-exponential), calculate error 



%% Display background-reduced desorption traces:
%Run coxlambertactualactual_BGreduce beforehand:
n_trace = 8;
cmap = parula(n_trace);

figure(10); clf;
hold on;
plot(new_x_0p4, new_y_0p4, 'o-');
plot(new_x_0p8, new_y_0p8, 'o-');
plot(new_x_1p2, new_y_1p2, 'o-');
plot(new_x_1p6, new_y_1p6, 'o-');
plot(new_x_2p0, new_y_2p0, 'o-');
plot(new_x_2p8, new_y_2p8, 'o-');
plot(new_x_4p0, new_y_4p0, 'o-');
plot(new_x_8p0, new_y_8p0, 'o-');
set(gca, 'ColorOrder', cmap, 'NextPlot', 'replacechildren');
legend('0.4', '0.8', '1.2', '1.6', '2.0', '2.8', '4.0', '8.0');
xlabel('Temperature (K)');
ylabel('Spectrometer Signal');
title('Desorption Traces, Signal vs Temperature')
hold off;

%% Take above background-reduced resorption traces, and change temperature to time
temp_init = 300; %K 
beta = 50; %K/s
time = @(x) (x-temp_init)/beta;

figure(2); clf;
hold on;
%plot(time(new_x_0p4), new_y_0p4, 'o-');
plot(time(new_x_0p8), new_y_0p8, 'o-');
plot(time(new_x_1p2), new_y_1p2, 'o-');
plot(time(new_x_1p6), new_y_1p6, 'o-');
%plot(time(new_x_2p0), new_y_2p0, 'o-');
%plot(time(new_x_2p8), new_y_2p8, 'o-');
%plot(time(new_x_4p0), new_y_4p0, 'o-');
%plot(time(new_x_8p0), new_y_8p0, 'o-');

set(gca, 'ColorOrder', cmap, 'NextPlot', 'replacechildren');
legend('0.8 * 10^{19}', '1.2 * 10^{19}', '1.6 * 10^{19}');
xlabel('Time (s)');
ylabel('Spectrometer Signal');
title('Desorption Traces, Signal vs Time')
hold off;


%% Determine coverage - dosage relationship. Copied from the first 'actual'.
%Run 'coxlambertscan_fig2cov.m' first. 


% Step 1: Remove duplicate x-values before interpolation
[unique_x, idx] = unique(normcoverage(:,1)); 
unique_y = normcoverage(idx,2);

% Step 2: Sort the unique values
SNC = sortrows([unique_x, unique_y], 1);

% Step 3: Create finer interpolation grid 
x_interp = linspace(min(SNC(:,1)), max(SNC(:,1)), 2000)';

% Step 4: Interpolate (with extrapolation disabled for safety)
y_interp = interp1(SNC(:,1), SNC(:,2), x_interp, 'pchip');

% Step 5: Combine interpolated data
NC_interp = [x_interp, y_interp];


%Adjusting scale based off of LEED
scaleFactor = (1/3)/(0.589311184524498);



% Step 6: Plotting
figure(4); clf;
hold on;
scatter(normcoverage(:,1), normcoverage(:,2) * scaleFactor, 'b'); % Original data
%plot(SNC(:,1), SNC(:,2), 'r-', 'LineWidth', 2); % Sorted unique data
plot(x_interp, y_interp * scaleFactor, 'r', 'LineWidth', 1.5); % Interpolated curve
hold off;

xlabel('Dosage (molecules/cm^2*10^{19})');
ylabel('Normalized Coverage');
title('Dosage vs Coverage');
%legend('Original Data', 'Interpolated Original Data', 'Location', 'best','FontSize',16);
set(groot, 'DefaultAxesFontSize', 16)
xlim([0 6])
grid on;
%% Now Arrhenius analysis, getting the areas to the right of each signal. 
%6/2 - ignore this integration method this sucks now 

N0_0p4 = 0.0908;
N0_0p8 = 0.2004;
N0_1p2 = 0.2877;
N0_1p6 = 0.3856; %5/21 
N0_2p0 = 0.43065;
N0_2p8 = 0.46168;
N0_4p0 = 0.47695;
N0_8p0 = 0.49657;

%5/29 - Signals must be normalized to 1; divide all of the signal by the
%max value of the top initial coverage, in this case the maximum of
%new_y_8p0.
cmap = parula(5);
[max_new_y, ~] = max(new_y_8p0);

new_N_0p4 = zeros(1,length(new_x_0p4)); %defines empty array, will be filled in with the areas (the remaining coverage) 

for i = 1:length(new_x_0p4)
    new_N_0p4(i) = trapz(new_y_0p4(1:length(new_x_0p4) + 1 - i)); %gets the area underneath the curve past the signal (new_y) at index 'i'
end
[new_max_0p4, ~] = max(new_N_0p4); %gets maximum of areas for below normalization 
new_N_0p4 = new_N_0p4 * (N0_0p4 / new_max_0p4); %trapz assumes a spacing of 1, so if the actual spacing is lower the reported area's going to be higher
 %more on that, it makes the maximum for new_N_ (the first value)be equal to the initial coverage, N0_

fprintf('NEW INSTANCE \n')
figure(5); clf;
hold on;
%legend('0.4', '0.8', '1.2', '1.6', '2.0', '2.8', '4.0', '8.0');
[max_new_0p4, ~] = max(new_y_0p4);
%plot(new_x_0p4 , new_N_0p4, 'o-', 'Color', cmap(1,:), 'DisplayName', '0.4 Remaining Coverage');
%plot(new_x_0p8, new_y_0p8 .* (N0_0p8/max_new_0p8));
%plot(new_x_0p4, new_y_0p4 .* (N0_0p4/max_new_y),'Color',cmap(1,:), 'DisplayName','0.4 Normalized Signal');

%set(gca, 'ColorOrder', cmap, 'NextPlot', 'replacechildren');

hold off;
legend; 

new_N_0p8 = zeros(1,length(new_x_0p8));
for i = 1:length(new_x_0p8)
    new_N_0p8(i) = trapz(new_y_0p8(1:length(new_x_0p8) + 1 - i)); %gets the area underneath the curve past the signal (new_y) at index 'i'
end
[new_max_0p8, ~] = max(new_N_0p8); %gets maximum of areas for below normalization 
new_N_0p8 = new_N_0p8 * (N0_0p8 / new_max_0p8); %trapz assumes a spacing of 1, so if the actual spacing is lower the reported area's going to be higher
 %more on that, it makes the maximum for new_N_ (the first value)be equal to the initial coverage, N0_

fprintf('NEW INSTANCE \n')
figure(5);
hold on;
%legend('0.4', '0.8', '1.2', '1.6', '2.0', '2.8', '4.0', '8.0');
[max_new_0p8, ~] = max(new_y_0p8);
%plot(new_x_0p8 , new_N_0p8, 'o-', 'Color', cmap(2,:),'DisplayName', '0.8 Remaining Coverage');
plot(new_x_0p8, new_y_0p8 .* (N0_0p8/max_new_y),'Color',cmap(1,:), 'DisplayName','0.8 Normalized Signal');

%plot(new_x_0p8, new_y_0p8 .* (N0_0p8/max_new_0p8));

set(gca, 'ColorOrder', cmap, 'NextPlot', 'replacechildren');

hold off;
new_N_1p2 = zeros(1,length(new_x_1p2));
for i = 1:length(new_x_1p2)
    new_N_1p2(i) = trapz(new_y_1p2(1:length(new_x_1p2) + 1 - i)); %gets the area underneath the curve past the signal (new_y) at index 'i'
end
[new_max_1p2, ~] = max(new_N_1p2); %gets maximum of areas for below normalization 
new_N_1p2 = new_N_1p2 * (N0_1p2 / new_max_1p2); %trapz assumes a spacing of 1, so if the actual spacing is lower the reported area's going to be higher
 %more on that, it makes the maximum for new_N_ (the first value)be equal to the initial coverage, N0_

fprintf('NEW INSTANCE \n')
figure(5); 
hold on;
%legend('0.4', '0.8', '1.2', '1.6', '2.0', '2.8', '4.0', '8.0');
[max_new_1p2, ~] = max(new_y_1p2);
%plot(new_x_1p2 , new_N_1p2, 'o-', 'Color',cmap(3,:), 'DisplayName','1.2 Remaining Coverage');
plot(new_x_1p2, new_y_1p2 .* (N0_1p2/max_new_y),'Color', cmap(2,:), 'DisplayName', '1.2 Normalized Signal');

%plot(new_x_1p2, new_y_1p2 .* (N0_1p2/max_new_1p2));

set(gca, 'ColorOrder', cmap, 'NextPlot', 'replacechildren');

hold off;

new_N_1p6 = zeros(1,length(new_x_1p6));

for i = 1:length(new_x_1p6)
    new_N_1p6(i) = trapz(new_y_1p6(1:length(new_x_1p6) + 1 - i)); %gets the area underneath the curve past the signal (new_y) at index 'i'
end
[new_max_1p6, ~] = max(new_N_1p6); %gets maximum of areas for below normalization 
new_N_1p6 = new_N_1p6 * (N0_1p6 / new_max_1p6); %trapz assumes a spacing of 1, so if the actual spacing is lower the reported area's going to be higher
 %more on that, it makes the maximum for new_N_ (the first value)be equal to the initial coverage, N0_

fprintf('NEW INSTANCE \n')
figure(5); 
hold on;
%legend('0.4', '0.8', '1.2', '1.6', '2.0', '2.8', '4.0', '8.0');

%plot(new_x_1p6 , new_N_1p6, 'o-', 'Color', cmap(4,:), 'DisplayName', '1.6 Normalized Converage');
plot(new_x_1p6, new_y_1p6 .* (N0_1p6/max_new_y),'Color', cmap(3,:), 'DisplayName','1.6 Normalized Signal');

set(gca, 'ColorOrder', cmap, 'NextPlot', 'replacechildren');

hold off;

new_N_2p0 = zeros(1,length(new_x_2p0));

for i = 1:length(new_x_2p0)
    new_N_2p0(i) = trapz(new_y_2p0(1:length(new_x_2p0) + 1 - i)); %gets the area underneath the curve past the signal (new_y) at index 'i'
end
[new_max_2p0, ~] = max(new_N_2p0); %gets maximum of areas for below normalization 
new_N_2p0 = new_N_2p0 * (N0_2p0 / new_max_2p0); %trapz assumes a spacing of 1, so if the actual spacing is lower the reported area's going to be higher
 %more on that, it makes the maximum for new_N_ (the first value)be equal to the initial coverage, N0_




fprintf('NEW INSTANCE \n')
figure(5); 
hold on;
%legend('0.4', '0.8', '1.2', '1.6', '2.0', '2.8', '4.0', '8.0');

%plot(new_x_2p0 , new_N_2p0, 'o-', 'Color', cmap(5,:), 'DisplayName', '2.0 Normalized Coverage');
%plot(new_x_2p0, new_y_2p0 .* (N0_2p0/max_new_y), 'Color', cmap(5,:), 'DisplayName', '2.0 Normalized Signal');

set(gca, 'ColorOrder', cmap, 'NextPlot', 'replacechildren');
%}
hold off;


%% 5/28 - Logarithmic plots. 

%arrh_x = [1/895,1/1100 , 1/1121, 1/1187];
%arrh_y = [,log(0.106/0.2), log(0.255/0.2), log(0.304/0.2)];

%Removing first point because of its overdependence on where the peak start
%is defined, and last point because it strays into the nonlinear region of coverage vs dosage. 
arrh_x = [1/955, 1/1148 , 1/1186];
arrh_y = [log(0.00227/0.2),log(0.16/0.2), log(0.26/0.2)];

p = polyfit(arrh_x, arrh_y, 1);

% Extract slope and intercept
slope = p(1);
intercept = p(2);

% Generate fit line
x_fit = linspace(min(arrh_x), max(arrh_x), 100);
y_fit = polyval(p, x_fit);

% Plot
figure(11); clf;
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

%% 5/29 - Logarithmic plots but with time as the independent variable. Need a term for heating rate. 
beta = 50; %K/s

arrh_x = [1/1148 , 1/1186];
arrh_y = [log(beta * 0.16/0.2), log(beta * 0.26/0.2)];

p = polyfit(arrh_x, arrh_y, 1);

% Extract slope and intercept
slope = p(1);
intercept = p(2);

% Generate fit line
x_fit = linspace(min(arrh_x), max(arrh_x), 100);
y_fit = polyval(p, x_fit);

% Plot
figure(11); clf;
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
