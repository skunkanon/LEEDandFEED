% TPD Simulation via Kinetic Monte Carlo (Hexagonal Lattice, Set A)
% Based on Meng & Weinberg (1994)

%clear; clc;

%% --- Parameters (Set A, Hexagonal Lattice) ---
L = 40;                      % Lattice size (L x L)
kB = 0.001987;               % Boltzmann constant [kcal/mol/K]
Ed0 = 31.6;                  % Desorption barrier for isolated molecule [kcal/mol]
eps_nn = 2.20;               % NN interaction energy [kcal/mol]
nu = 1e13;                   % Pre-exponential factor [1/s]
beta = 5;                    % Heating rate [K/s]
T = 300;                     % Initial temperature [K]
T_max = 600;                 % Max temperature [K]
t = 0;                       % Initial time [s]
theta0 = 0.5;                % Initial coverage

%% --- Initialization ---
lattice = zeros(L);          % 0 = empty, 1 = occupied
N_sites = L^2;
N_ads = round(theta0 * N_sites);
ind_all = randperm(N_sites, N_ads);
lattice(ind_all) = 1;

%% --- Neighbor indexing (hexagonal lattice, periodic boundaries) ---
neighbors = @(i,j) [mod([i-2, i-2, i-1, i+1, i, i]-1, L)+1;
                    mod([j, j+1, j-1, j+1, j-1, j+1]-1, L)+1];

%% --- Simulation Setup ---
dT = 0.05; % Finer increments for smooth spectrum % Temperature increment per iteration [K]
T_range = T:dT:T_max;
TPD_data = zeros(length(T_range),2);
idx = 1;

for T = T_range
    % Surface relaxation (random hops)
    for sweep = 1:5000
        i = randi(L); j = randi(L);
        if lattice(i,j) == 1
            nei = neighbors(i,j);
            dir = randi(6);
            ni = nei(1,dir); nj = nei(2,dir);
            if lattice(ni, nj) == 0
                lattice(ni, nj) = 1;
                lattice(i, j) = 0;
            end
        end
    end

    % Compute desorption rates
    [rows, cols] = find(lattice == 1);
    N_ads = numel(rows);
    rates = zeros(N_ads, 1);

    for k = 1:N_ads
        i = rows(k); j = cols(k);
        nei = neighbors(i,j);
        nn_count = sum(diag(lattice(nei(1,:), nei(2,:))));
        Ed = Ed0 - nn_count * eps_nn;
        rates(k) = nu * exp(-Ed / (kB * T));
    end

    R_tot = sum(rates);
    dt = dT / beta;
    N_desorb = min(poissrnd(R_tot * dt), N_ads);

    if N_desorb > 0
        cum_rates = cumsum(rates) / R_tot;
        for des = 1:N_desorb
            r = rand();
            des_index = find(cum_rates >= r, 1);
            lattice(rows(des_index), cols(des_index)) = 0;
            rates(des_index) = 0; % avoid choosing the same molecule again
            cum_rates = cumsum(rates) / sum(rates);
        end
    end

    % Record TPD data
    coverage = sum(lattice(:)) / N_sites;
    TPD_data(idx,:) = [T, coverage];
    idx = idx + 1;

    if coverage <= 0
        break; % End simulation if surface is empty
    end
end

% Remove unused rows if simulation ends early
TPD_data(idx:end,:) = [];

%% --- Plot Results ---
temp = TPD_data(:,1);
coverage = TPD_data(:,2);
rate = -gradient(coverage, temp);

windowSize = 10; % Smoothing window
rate_smooth = movmean(rate, windowSize);

figure;
scatter(temp, rate_smooth, 20, 'filled');
xlabel('Temperature (K)');
ylabel('Desorption rate (arb. units)');
title('TPD Spectrum (Hexagonal, Set A, Scatterplot)');
grid on;
