%% 7/31 2 AM 
fprintf('NEW INSTANCE \n')
figure(1); clf;
hold on;

runs=1;

for i = 1:runs
setA_MC = MONTECARLO_HEX_JUL25(0.25,0,0, 5);
scatter(setA_MC(:,1), -gradient(setA_MC(:,2), setA_MC(:,1)), 20, 'filled', 'r', 'DisplayName', 'Coverage = 0.25');
end
[~ ,p25_tmp, p25_rate, ~] = polyani_wigner(5, 300, 31.6*4184, 0 , 10^13, 0.25,600 );
plot(p25_tmp,p25_rate, 'r');

for i = 1:runs
setA_MC = MONTECARLO_HEX_JUL25(0.5,0,0,5);
scatter(setA_MC(:,1), -gradient(setA_MC(:,2), setA_MC(:,1)), 20, 'filled', 'g', 'DisplayName', 'Coverage = 0.50');

end
[~ ,p50_tmp, p50_rate, ~] = polyani_wigner(5, 300, 31.6*4184, 10 , 10^13, 0.5,600 );
plot(p50_tmp,p50_rate, 'g');

for i = 1:runs
setA_MC = MONTECARLO_HEX_JUL25(0.75,0,0,5);
scatter(setA_MC(:,1), -gradient(setA_MC(:,2), setA_MC(:,1)), 20, 'filled', 'b', 'DisplayName', 'Coverage = 0.75');

end
[~ ,p75_tmp, p75_rate, ~] = polyani_wigner(5, 300, 31.6*4184, 0 , 10^13, 0.75,600 );
plot(p75_tmp,p75_rate, 'b');

for i = 1:runs
setA_MC = MONTECARLO_HEX_JUL25(1,0,0,5);
scatter(setA_MC(:,1), -gradient(setA_MC(:,2), setA_MC(:,1)), 20, 'filled', 'm', 'DisplayName', 'Coverage = 1');

end
[~ ,p0_tmp, p0_rate, ~] = polyani_wigner(5, 300, 31.6*4184, 0 , 10^13, 1,600 );
plot(p0_tmp,p0_rate, 'm');


hold off;
xlabel('Temperature (K)');
ylabel('Desorption rate');
title('HEX LATTICE FIRST ORDER ZERO INTERACTION ENERGY JULY 30');
legend show;
grid on;

%general shape checks out but the peak temperatures are off 7/30. could
%revisit. their 2nd order paper has the triangular lattice, but first order
%doesn't. 

%7/30 - changed to 1st order 
%7/31 - added debug statements, can only really replicate the simple
%first-order desorption from polanyi-wigner without energy dependence. meng
%and weinberg can get BENT 
%7/31, 3 AM: making it diffuse 10x more with zero interaction energy does
%nothing

%% 7/31 3 AM 
fprintf('NEW INSTANCE \n')

figure(2); clf; hold on;

runs = 1;

for i = 1:runs
setA_MC_2 = MONTECARLO_HEX_JUL25(0.5,2,0.1, 3);
scatter(setA_MC_2(:,1), -gradient(setA_MC_2(:,2), setA_MC_2(:,1)), 20, 'filled', 'r', 'DisplayName', 'Coverage = 0.25');
end



for i = 1:runs
setA_MC_2 = MONTECARLO_HEX_JUL25(0.5,2,0.4, 3);
scatter(setA_MC_2(:,1), -gradient(setA_MC_2(:,2), setA_MC_2(:,1)), 20, 'filled', 'b', 'DisplayName', 'Coverage = 0.25');
end



for i = 1:runs
setA_MC_2 = MONTECARLO_HEX_JUL25(0.5,2,0.8, 3);
scatter(setA_MC_2(:,1), -gradient(setA_MC_2(:,2), setA_MC_2(:,1)), 20, 'filled', 'g', 'DisplayName', 'Coverage = 0.25');
end


for i = 1:runs
setA_MC_2 = MONTECARLO_HEX_JUL25(0.5,2,1.2, 3);
scatter(setA_MC_2(:,1), -gradient(setA_MC_2(:,2), setA_MC_2(:,1)), 20, 'filled', 'm', 'DisplayName', 'Coverage = 0.25');
end




hold off;

%7/31 4 AM - '94 publication: lower peak definitely grows as coverage
%increases. also would be more efficient to average way smaller lattice
%sizes, probably but dont want to code that yet 


%% 7/31 3 PM- AVERAGING RESULTS FROM MULTIPLE RUNS WITH SMALLER LATTICE SIZE 

figure(3); clf; hold on; 

runs = 20; % Increase number of runs
eps_nn = 2; 
eps_nnn = 0;
RATIO = 3; 
theta = 0.5;


% Initialize arrays to store all results
all_temperatures = [];
all_gradients = [];

for i = 1:runs
    [setA_MC_2, ~, ~] = MONTECARLO_HEX_JUL25(theta,eps_nn,eps_nnn, RATIO);
    
    % Calculate gradient for this run
    temperatures = setA_MC_2(:,1);
    coverage = setA_MC_2(:,2);
    gradient_vals = -gradient(coverage, temperatures);
    
    % Store results
    all_temperatures = [all_temperatures; temperatures];
    all_gradients = [all_gradients; gradient_vals];
end

% Calculate mean across all runs
unique_temps = unique(all_temperatures);
mean_gradients = zeros(size(unique_temps));

for temp_idx = 1:length(unique_temps)
    temp = unique_temps(temp_idx);
    temp_indices = find(all_temperatures == temp);
    mean_gradients(temp_idx) = mean(all_gradients(temp_indices));
end

% Plot averaged result
scatter(unique_temps, mean_gradients, 20, 'filled', 'r');

hold off;


%% 8/4 - TESTING ENERGY VS COVERAGE PLOT
fprintf('NEW INSTANCE \n');


runs = 20;
eps_nn = 2; 
eps_nnn = 0;
RATIO = 3; 
theta = 0.5;

all_coverages = [];
all_Eds = [];

for i = 1:runs
    % Get coverage and Ed arrays from the function
    [~, coverage_array, Ed_array] = MONTECARLO_HEX_JUL25(theta, eps_nn, eps_nnn, RATIO);
    
    % Store results
    all_coverages = [all_coverages; coverage_array(:)];
    all_Eds = [all_Eds; Ed_array(:)];
end

% Average coverage for each unique Ed value
unique_Eds = unique(all_Eds(~isnan(all_Eds)));
mean_coverages = zeros(size(unique_Eds));

for idx = 1:length(unique_Eds)
    ed_val = unique_Eds(idx);
    indices = find(all_Eds == ed_val);
    mean_coverages(idx) = mean(all_coverages(indices));
end

% Plot averaged result
figure(4); clf; hold on;
scatter(mean_coverages, unique_Eds, 5, 'b', 'filled');
ylabel('Average Activation Energy, E_d (kcal/mol)');
xlabel('Coverage, \theta');
title('Coverage vs Average Activation Energy (Averaged)');
grid on;
hold off;
%% 8/4 - PLOTTING BOTH RATE v TEMP & COVERAGEvEA 

% Combined Monte Carlo TPD Analysis: Rate vs Temp and Ed vs Coverage

figure(5); clf;

runs = 50;
eps_nn = 0; 
eps_nnn = 2;
RATIO = 3; 
theta = 0.5;

% Initialize arrays to store all results
all_temperatures = [];
all_gradients = [];
all_coverages = [];
all_Eds = [];

for i = 1:runs
    [setA_MC_2, coverage_array, Ed_array] = MONTECARLO_HEX_JUL25(theta, eps_nn, eps_nnn, RATIO);
    
    % --- For Rate vs Temp ---
    temperatures = setA_MC_2(:,1);
    coverage = setA_MC_2(:,2);
    gradient_vals = -gradient(coverage, temperatures);
    all_temperatures = [all_temperatures; temperatures];
    all_gradients = [all_gradients; gradient_vals];
    
    % --- For Ed vs Coverage ---
    all_coverages = [all_coverages; coverage_array(:)];
    all_Eds = [all_Eds; Ed_array(:)];
end

% --- Process Rate vs Temp ---
unique_temps = unique(all_temperatures);
mean_gradients = zeros(size(unique_temps));
for temp_idx = 1:length(unique_temps)
    temp = unique_temps(temp_idx);
    temp_indices = find(all_temperatures == temp);
    mean_gradients(temp_idx) = mean(all_gradients(temp_indices));
end

% --- Process Ed vs Coverage ---
unique_coverages = unique(all_coverages(~isnan(all_coverages)));
mean_Eds = zeros(size(unique_coverages));
for idx = 1:length(unique_coverages)
    cov_val = unique_coverages(idx);
    indices = find(all_coverages == cov_val);
    mean_Eds(idx) = mean(all_Eds(indices));
end

% --- Plot as subplots ---
subplot(1,2,1); % Left plot: Rate vs Temp
scatter(unique_temps, mean_gradients, 20, 'filled', 'r');
xlabel('Temperature (K)');
ylabel('-d\theta/dT');
title('Averaged Rate vs Temperature');
grid on;

subplot(1,2,2); % Right plot: Ed vs Coverage
scatter(unique_coverages, mean_Eds, 5, 'b', 'filled');
xlabel('Coverage, \theta');
ylabel('Average Activation Energy, E_d (kcal/mol)');
title('Average Activation Energy vs Coverage (Averaged)');
grid on;

sgtitle('Monte Carlo TPD Simulation Results');
%%
