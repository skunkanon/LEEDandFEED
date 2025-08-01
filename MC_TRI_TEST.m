function setA_MC = MC_TRI_TEST(theta0)


L = 200;                      % Lattice size (L x L)
kB = 0.001987;               % Boltzmann constant [kcal/mol/K]
Ed0 = 30.0;                  % Desorption barrier for isolated molecule [kcal/mol]
eps_nn = 0.00;               % NN interaction energy [kcal/mol] 
nu = 1e15;                   % Pre-exponential factor [1/s]
beta = 1;                    % Heating rate [K/s] 
T = 300;                     % Initial temperature [K]
T_max = 600;                 % Max temperature [K]
t = 0;                       % Initial time [s]


% --- Initialization ---
lattice = zeros(L);          % 0 = empty, 1 = occupied
N_sites = L^2;
N_ads = round(theta0 * N_sites);
ind_all = randperm(N_sites, N_ads);
lattice(ind_all) = 1;

% --- Neighbor indexing (triangular lattice, periodic boundaries) ---
% For triangular lattice, need different offsets for even/odd rows
even_neigh_offsets = [0 -1; 0 1; -1 -1; -1 0; 1 -1; 1 0];
odd_neigh_offsets  = [0 -1; 0 1; -1 0; -1 1; 1 0; 1 1];

% --- Simulation Setup ---
dT = 1; % Temperature increment per iteration [K]
T_range = T:dT:T_max;
setA_MC = zeros(length(T_range),2);
idx = 1;

for T = T_range
    % Surface relaxation (random hops)
    for sweep = 1:1000
        i = randi(L); j = randi(L);
        if lattice(i,j) == 1
            % Get neighbor offsets based on row parity
            if mod(i, 2) == 0
                neigh_offsets = even_neigh_offsets;
            else
                neigh_offsets = odd_neigh_offsets;
            end
            
            dir = randi(6);
            ni = i + neigh_offsets(dir,1);
            nj = j + neigh_offsets(dir,2);
            
            % Periodic boundary conditions
            if ni < 1, ni = L; elseif ni > L, ni = 1; end
            if nj < 1, nj = L; elseif nj > L, nj = 1; end
            
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
        
        % Count neighbors using correct triangular indexing
        if mod(i, 2) == 0
            neigh_offsets = even_neigh_offsets;
        else
            neigh_offsets = odd_neigh_offsets;
        end
        
        nn_count = 0;
        for n = 1:6
            ni = i + neigh_offsets(n,1);
            nj = j + neigh_offsets(n,2);
            
            % Periodic boundary conditions
            if ni < 1, ni = L; elseif ni > L, ni = 1; end
            if nj < 1, nj = L; elseif nj > L, nj = 1; end
            
            if lattice(ni, nj) == 1
                nn_count = nn_count + 1;
            end
        end
        
        Ed = Ed0 - nn_count * eps_nn;
        rates(k) = nu * exp(-Ed / (kB * T));
    end

    R_tot = sum(rates);

    % Calculate number of desorptions in this temperature step
    dt = dT / beta;
    N_desorb = min(poissrnd(R_tot * dt), N_ads);

    if N_desorb > 0
        cum_rates = cumsum(rates) / R_tot;
        for des = 1:N_desorb
            r = rand();
            des_index = find(cum_rates >= r, 1);
            lattice(rows(des_index), cols(des_index)) = 0;
            rates(des_index) = 0;
            cum_rates = cumsum(rates) / sum(rates);
        end
    end

    % Record TPD data
    coverage = sum(lattice(:)) / N_sites;
    setA_MC(idx,:) = [T, coverage];
    idx = idx + 1;

    if coverage <= 0
        break; % End simulation if surface is empty
    end
end

% Remove unused rows if simulation ends early
setA_MC(idx:end,:) = [];

end 