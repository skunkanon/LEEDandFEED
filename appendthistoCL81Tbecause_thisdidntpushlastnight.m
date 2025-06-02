%% 4-3 - Deriving coverage from LEED; just testing if the digitization worked 

scatter(normcoverage(:, 1), normcoverage(:, 2)); %regular plot function loops back in on itself 

%% very crude method for calibrating the axes lol. bring the last number in x_interp at line 15 to 60,000 and find where it's  1.4 

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

% Step 6: Plotting
figure;
hold on;
scatter(normcoverage(:,1), normcoverage(:,2), 20, 'b', 'filled'); % Original data
plot(SNC(:,1), SNC(:,2), 'r-', 'LineWidth', 2); % Sorted unique data
plot(x_interp, y_interp, 'g--', 'LineWidth', 1.5); % Interpolated curve
hold off;

xlabel('Normalized Coverage');
ylabel('LEED Intensity');
title('Coverage vs. LEED Intensity (Interpolated)');
legend('Original Data', 'Sorted Unique', 'Interpolated', 'Location', 'best');
grid on;

%% very crude method for calibrating the axes lol. bring the last number in x_interp at line 15 to 60,000 and find where it's 1.4 

% Step 1: Remove duplicate x-values before interpolation
[unique_x, idx] = unique(TDS_raw(:,1)); 
unique_y = TDS_raw(idx,2);

% Step 2: Sort the unique values
SNC = sortrows([unique_x, unique_y], 1);

% Step 3: Create finer interpolation grid 
x_interp = linspace(min(SNC(:,1)), max(SNC(:,1)), 60000)';

% Step 4: Interpolate (with extrapolation disabled for safety)
y_interp = interp1(SNC(:,1), SNC(:,2), x_interp, 'pchip');

% Step 5: Combine interpolated data
NC_interp = [x_interp, y_interp];

% Step 6: Plotting
figure;
hold on;
scatter(TDS_raw(:,1), TDS_raw(:,2), 20, 'b', 'filled'); % Original data
plot(SNC(:,1), SNC(:,2), 'r-', 'LineWidth', 2); % Sorted unique data
plot(x_interp, y_interp, 'g--', 'LineWidth', 1.5); % Interpolated curve
hold off;

xlabel('Normalized Coverage');
ylabel('LEED Intensity');
title('Coverage vs. LEED Intensity (Interpolated)');
legend('Original Data', 'Sorted Unique', 'Interpolated', 'Location', 'best');
grid on;

%% 


% 1.4 x 10^19 corresponds to 1/3 coverage. Interpolated gives 0.5893. Need
%to scale everything down likewise. 

scaleFactor = (1/3)/(0.589311184524498);

%coverageActual = [x_interp .* 10^19, y_interp .* scaleFactor];

scaleFactor_TDS = (1/3)/(0.445454545454545);
scatter(TDS_raw(:,1), TDS_raw(:,2));
figure; 
hold on; 
TDS_x = x_interp .* 10^19;
TDS_y = y_interp .* scaleFactor;

scatter(TDS_x, TDS_y);


% Formatting
xlabel('Dose (molecules/m²)');
ylabel('Coverage (ML)');
title('TDS Coverage vs. Dose');
%legend('Location', 'southeast', 'FontSize', 10);


xlim([0, max(TDS_x)]);
 
%plot(coverageActual(:,1), coverageActual(:,2));
hold off;
%%
scaleFactor_TDS = (1/3)/(0.445454545454545); 

% Generate TDS data (scaled to ML)
TDS_x = x_interp .* 1e19; % Dose in molecules/m²
TDS_y = y_interp .* scaleFactor_TDS; % Coverage in ML

% Find x-values for key coverages
target_y = [1/3, 9/16]; % 0.333 and 0.5625 ML
x_at_1over3 = interp1(TDS_y, TDS_x, target_y(1), 'linear', 'extrap');
x_at_9over16 = interp1(TDS_y, TDS_x, target_y(2), 'linear', 'extrap');

% Plot
figure; 
plot(TDS_x, TDS_y, 'b-', 'LineWidth', 2);
hold on;

% Highlight critical doses
scatter(x_at_1over3, target_y(1), 100, 'r', 'filled', 'DisplayName', '1/3 ML');
scatter(x_at_9over16, target_y(2), 100, 'g', 'filled', 'DisplayName', '9/16 ML');

% Annotate x-values
text(x_at_1over3, target_y(1), ...
    sprintf('Dose: %.2f x10^{19}', x_at_1over3/1e19), ...
    'VerticalAlignment', 'bottom', 'FontSize', 10);
text(x_at_9over16, target_y(2), ...
    sprintf('Dose: %.2f x10^{19}', x_at_9over16/1e19), ...
    'VerticalAlignment', 'bottom', 'FontSize', 10);

% Formatting
xlabel('Dose (molecules/m²)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Coverage (ML)', 'FontSize', 12, 'FontWeight', 'bold');
title('TDS Coverage vs. Dose', 'FontSize', 14);
legend('Location', 'southeast', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 10, 'LineWidth', 1.5);
xlim([0, max(TDS_x)]);
hold off;


%% Getting coverage
%Fig 7: ln(beta * f / s) vs 1/T

fprintf('\n NEW INSTANCE \n');
%Getting peak temp

[Tp, Tp_index] = deal(zeros(1,numel(data_cell)));

%Getting peak temps 
for i = 1:numel(data_cell) 
    x = data_cell{i}(:, 1);  % First column (x-values)
    y = data_cell{i}(:, 2);  % Second column (y-values)
    [~, Tp_index(i)] = max(y);
    Tp(i) = x(Tp_index(i));
end

beta = 50; %K/s


%% 4-8 FIG FIVE 

