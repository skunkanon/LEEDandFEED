function [setA_hex, coverage_array, Ed_array] = MONTECARLO_HEX_JUL25(theta0, eps_nn, eps_nnn, RATIO)


L = 100;
kB = 0.001987;         % [kcal/mol/K]
Ed0 = 31.6;            % Isolated desorption barrier [kcal/mol]
%eps_nn = 0;          % NN interaction energy [kcal/mol]
%eps_nnn = 0;           % NNN interaction energy [kcal/mol]
nu = 1e13;             % Pre-exponential factor [1/s]
beta = 5;              % Heating rate [K/s]
T = 200; T_max = 700; dT = 1;

% --- Initialize lattice ---
lattice = zeros(L);
N_sites = L^2;
N_ads = round(theta0 * N_sites);
lattice(randperm(N_sites, N_ads)) = 1;

% Debug output for initialization
fprintf('Initialization: L=%d, N_sites=%d, N_ads=%d, theta0=%.3f\n', L, N_sites, N_ads, theta0);

% --- Hexagonal lattice offsets ---
even_nn = [0 -1; 0 1; -1 -1; -1 0; 1 -1; 1 0];
odd_nn  = [0 -1; 0 1; -1 0; -1 1; 1 0; 1 1];
even_nnn = [-2 0; -1 -2; -1 2; 0 -2; 0 2; 1 -2; 1 2; 2 0];
odd_nnn  = [-2 0; -1 -1; -1 1; 0 -2; 0 2; 1 -1; 1 1; 2 0];

% --- Output array ---
T_range = T:dT:T_max;
setA_hex = zeros(length(T_range), 2);
coverage_array = zeros(length(T_range), 1); % Add: coverage at each T
Ed_array = zeros(length(T_range), 1);       % Add: avg Ed at each T
idx = 1;

for T = T_range
    % --- Surface diffusion with Metropolis bias ---
    for sweep = 1:(RATIO*N_ads) %RATIO ORIGINALLY 0.5
        i = randi(L); j = randi(L);
        if lattice(i,j) == 1
            if mod(i,2)==0
                offsets = even_nn;
            else
                offsets = odd_nn;
            end

            dir = randi(6);
            ni = mod(i + offsets(dir,1)-1, L) + 1;
            nj = mod(j + offsets(dir,2)-1, L) + 1;
            if lattice(ni,nj) == 0
                Ed_initial = Ed0 - eps_nn * count_nn(i, j, lattice, L, even_nn, odd_nn);
                Ed_final   = Ed0 - eps_nn * count_nn(ni, nj, lattice, L, even_nn, odd_nn);
                deltaE = Ed_final - Ed_initial;
                p_move = min(1, exp(-deltaE / (kB * T)));
                if rand() < p_move
                    lattice(ni, nj) = 1;
                    lattice(i, j) = 0;
                end
            end
        end
    end

    % --- Desorption rates ---
    [rows, cols] = find(lattice == 1);
    N_ads = numel(rows);
    if N_ads == 0
        fprintf('All molecules desorbed at T=%.1f K\n', T);
        break;
    end
    rates = zeros(N_ads,1);
    Eds = zeros(N_ads,1); % Add: store Ed for each adsorbate

    for k = 1:N_ads
        i = rows(k); j = cols(k);
        if mod(i,2)==0
            nn_off = even_nn;
            nnn_off = even_nnn;
        else
            nn_off = odd_nn;
            nnn_off = odd_nnn;
        end

        % Count nearest neighbors
        nn_count = 0;
        for n = 1:6
            ni = mod(i + nn_off(n,1)-1, L) + 1;
            nj = mod(j + nn_off(n,2)-1, L) + 1;
            if lattice(ni,nj) == 1
                nn_count = nn_count + 1;
            end
        end

        % Count next-nearest neighbors
        nnn_count = 0;
        for n = 1:8
            ni = mod(i + nnn_off(n,1)-1, L) + 1;
            nj = mod(j + nnn_off(n,2)-1, L) + 1;
            if lattice(ni,nj) == 1
                nnn_count = nnn_count + 1;
            end
        end

        % Effective desorption energy
        Ed = Ed0 - eps_nn * nn_count - eps_nnn * nnn_count;
        Eds(k) = Ed; % Add: store Ed
        rates(k) = nu * exp(-Ed / (kB * T));
    end

    % --- Desorption events ---
    dt = dT / beta;
    R_tot = sum(rates);
    

    N_desorb = min(poissrnd(R_tot * dt), N_ads);

    if N_desorb > 0
        % Create cumulative rate array
        cum_rates = cumsum(rates) / R_tot;
        
        for d = 1:N_desorb
            r = rand();
            idx_sel = find(cum_rates >= r, 1);
            
            if isempty(idx_sel)
                idx_sel = N_ads; % Fallback
            end
            
            % Remove molecule
            lattice(rows(idx_sel), cols(idx_sel)) = 0;
            
            % Update arrays
            rows(idx_sel) = [];
            cols(idx_sel) = [];
            rates(idx_sel) = [];
            N_ads = N_ads - 1;
            
            % Recalculate cumulative rates if there are remaining molecules
            if N_ads > 0
                R_tot = sum(rates);
                if R_tot > 0
                    cum_rates = cumsum(rates) / R_tot;
                else
                    break;
                end
            else
                break;
            end
        end
    end



    %{
    % --- Post-desorption equilibration step ---
    for sweep = 1:(RATIO * 3 *N_ads) %ORIGINALLY 1.5
        i = randi(L); j = randi(L);
        if lattice(i,j) == 1
            if mod(i,2)==0
                offsets = even_nn;
            else
                offsets = odd_nn;
            end

            dir = randi(6);
            ni = mod(i + offsets(dir,1)-1, L) + 1;
            nj = mod(j + offsets(dir,2)-1, L) + 1;
            if lattice(ni,nj) == 0
                Ed_initial = Ed0 - eps_nn * count_nn(i, j, lattice, L, even_nn, odd_nn);
                Ed_final   = Ed0 - eps_nn * count_nn(ni, nj, lattice, L, even_nn, odd_nn);
                deltaE = Ed_final - Ed_initial;
                p_move = min(1, exp(-deltaE / (kB * T)));
                if rand() < p_move
                    lattice(ni, nj) = 1;
                    lattice(i, j) = 0;
                end
            end
        end
    end
    %} 
    % --- Record coverage vs temperature ---
    coverage = sum(lattice(:)) / N_sites;
    setA_hex(idx,:) = [T, coverage];
    coverage_array(idx) = coverage; % Add: store coverage
    if N_ads > 0
        Ed_array(idx) = mean(Eds); % Add: store mean Ed
    else
        Ed_array(idx) = NaN;
    end
    
    % Debug output every 50 temperature steps
    if mod(idx, 50) == 0 || coverage < 0.1
        fprintf('T=%.1f K, Coverage=%.4f, N_ads=%d\n', T, coverage, N_ads);
    end
    
    idx = idx + 1;

    if coverage <= 0
        fprintf('Complete desorption at T=%.1f K\n', T);
        break;
    end
end

% Trim unused rows from output array
setA_hex(idx:end,:) = [];
coverage_array(idx:end) = [];
Ed_array(idx:end) = [];

% Optionally, output as a cell array
coverage_Ed_cell = {coverage_array, Ed_array};

% Final debug output
fprintf('Simulation complete. Final coverage: %.4f\n', setA_hex(end,2));

end

% --- Helper function ---
function nn_count = count_nn(i, j, lattice, L, even_nn, odd_nn)
    if mod(i,2)==0
        offsets = even_nn;
    else
        offsets = odd_nn;
    end

    nn_count = 0;
    for n = 1:6
        ni = mod(i + offsets(n,1) - 1, L) + 1;
        nj = mod(j + offsets(n,2) - 1, L) + 1;
        if lattice(ni, nj) == 1
            nn_count = nn_count + 1;
        end
    end
end
