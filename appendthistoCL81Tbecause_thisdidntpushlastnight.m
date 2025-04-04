%% 4-3 - Deriving coverage from LEED; just testing if the digitization worked 

scatter(normcoverage(:, 1), normcoverage(:, 2)); %regular plot function loops back in on itself 

%% very crude method for calibrating the axes lol. bring the last number in x_interp at line 15 to 60,000 and find where it's 1.4 

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


%% 
% 1.4 x 10^19 corresponds to 1/3 coverage. Interpolated gives 0.5893. Need
%to scale everything down likewise. 

scaleFactor = (1/3)/(0.589311184524498);

coverageActual = [x_interp .* 10^19, y_interp .* scaleFactor];

scaleFactor_TDS = (1/3)/(0.445454545454545);

figure; 
hold on; 
TDS_x = TDS_raw(:,1) * 10^19;
TDS_y = TDS_raw(:,2) * scaleFactor_TDS;

plot(TDS_x, TDS_y, 'bo');

 
plot(coverageActual(:,1), coverageActual(:,2));
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


