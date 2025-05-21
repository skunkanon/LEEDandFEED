%5/21 - Splitting up coxlambertACTUAL, redoing scans, making it more
%organized. 
%eventually make theoretical curves with calculated data (Ea and
%pre-exponential), calculate error 



%% Display background-reduced desorption traces:
%Run coxlambertactualactual_BGreduce beforehand:
n_trace = 8;
cmap = parula(n_trace);

figure(1); clf;
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

new_x_0p4_time = time(new_x_0p4);


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
figure(3); clf;
hold on;
scatter(normcoverage(:,1), normcoverage(:,2) * scaleFactor, 'b'); % Original data
%plot(SNC(:,1), SNC(:,2), 'r-', 'LineWidth', 2); % Sorted unique data
plot(x_interp, y_interp * scaleFactor, 'r', 'LineWidth', 1.5); % Interpolated curve
hold off;

xlabel('Dosage (molecules/cm^2  * 10^{19})');
ylabel('Normalized Coverage');
title('Dosage vs Coverage');
legend('Original Data', 'Interpolated Original Data', 'Location', 'best');
grid on;