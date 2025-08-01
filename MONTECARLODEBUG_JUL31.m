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

runs = 40; % Increase number of runs
eps_nn = 2; 
eps_nnn = 0;
RATIO = 3; 
theta = 0.5;


% Initialize arrays to store all results
all_temperatures = [];
all_gradients = [];

for i = 1:runs
    setA_MC_2 = MONTECARLO_HEX_JUL25(theta,eps_nn,eps_nnn, RATIO);
    
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