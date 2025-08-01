% Monte Carlo simulation of first-order desorption on a 40x40 hexagonal lattice (parameter set E)
L = 40;
Beta = 5;               % Heating rate (K/s)
k0 = 1e13;              % Pre-exponential factor (s^-1)
Ed0_kcal = 31.6;        % Ed^0 in kcal/mol (isolated desorption barrier)
E_nn_kcal = 1.99;       % nearest-neighbor interaction (kcal/mol)
E_nnn_kcal = 0.39;      % next-nearest-neighbor interaction (kcal/mol)
theta0 = 0.5;





% Convert energies to units of K (divide energy in J by Boltzmann's constant)
Ed0_K   = Ed0_kcal * 4184 / 8.314; 
E_nn_K  = E_nn_kcal * 4184 / 8.314;
E_nnn_K = E_nnn_kcal * 4184 / 8.314;
% Precompute neighbor and next-neighbor indices for each site (periodic boundary conditions)
N_sites = L * L;
neighbors = zeros(N_sites, 6);
nextneighbors = zeros(N_sites, 12);
% Axial coordinate neighbor differences for hex lattice
neighbor_diffs = [1 0; -1 0; 0 1; 0 -1; 1 -1; -1 1];
second_diffs = [2 0; -2 0; 0 2; 0 -2; 1 1; -1 -1; 2 -1; -2 1; 1 -2; -1 2; 2 -2; -2 2];
for i = 1:L
    for j = 1:L
        idx = sub2ind([L L], i, j);
        % Neighbors (with wrapping)
        for k = 1:size(neighbor_diffs,1)
            di = neighbor_diffs(k,1); dj = neighbor_diffs(k,2);
            ni = mod(i-1+di, L) + 1;
            nj = mod(j-1+dj, L) + 1;
            neighbors(idx, k) = sub2ind([L L], ni, nj);
        end
        % Next-nearest neighbors
        for k = 1:size(second_diffs,1)
            di = second_diffs(k,1); dj = second_diffs(k,2);
            ni = mod(i-1+di, L) + 1;
            nj = mod(j-1+dj, L) + 1;
            nextneighbors(idx, k) = sub2ind([L L], ni, nj);
        end
    end
end

% Prepare figure for TPD curve
figure; hold on;
colors = lines(1);
N0 = round(theta0 * N_sites);
% Initialize lattice with N0 adsorbates at random sites
occ_flat = zeros(N_sites, 1);
perm = randperm(N_sites);
occ_flat(perm(1:N0)) = 1;
% Initial surface relaxation via fast diffusion (Metropolis Monte Carlo)
T_current = 300;  % assume initial temperature 300 K for diffusion equilibration
num_relax_attempts = 100 * N0;   % number of diffusion hop attempts for initial relaxation
for attempt = 1:num_relax_attempts
    occ_indices = find(occ_flat == 1);
    if isempty(occ_indices), break; end
    p = occ_indices(randi(length(occ_indices)));       % choose a random occupied site index
    % choose a random neighbor site of p
    neigh_list = neighbors(p, :);
    q = neigh_list(randi(6));
    if occ_flat(q) == 1
        continue;  % neighbor occupied, failed hop
    end
    % Compute desorption energy (Ed) of particle at p (initial) and if moved to q (final)
    % Initial Ed at p:
    Nn_p  = sum(occ_flat(neighbors(p,:)));        % occupied nearest neighbors of p
    Nnn_p = sum(occ_flat(nextneighbors(p,:)));    % occupied next-nearest neighbors of p
    Ed_p = Ed0_K - Nn_p * E_nn_K - Nnn_p * E_nnn_K;
    % Final Ed at q (after moving): remove p, add at q, then count neighbors of q
    occ_flat(p) = 0;
    Nn_q = sum(occ_flat(neighbors(q,:)));
    Nnn_q = sum(occ_flat(nextneighbors(q,:)));
    Ed_q = Ed0_K - Nn_q * E_nn_K - Nnn_q * E_nnn_K;
    % Metropolis acceptance for diffusion
    if Ed_q >= Ed_p || rand() < exp((Ed_q - Ed_p) / T_current)
        % Accept the hop: particle moves from p to q
        occ_flat(q) = 1;
        % (Note: p is already set to 0 above)
    else
        % Reject: restore p
        occ_flat(p) = 1;
    end
end

% Temperature-programmed desorption simulation (KMC)
T_current = 300;
t = 0;
occ_count = sum(occ_flat);
max_events = occ_count;
T_record = zeros(max_events, 1);
rate_record = zeros(max_events, 1);
event = 0;
while occ_count > 0
    % Compute desorption rates for all occupied molecules
    Nn_all  = sum(occ_flat(neighbors), 2);
    Nnn_all = sum(occ_flat(nextneighbors), 2);
    Ed_all = Ed0_K - Nn_all * E_nn_K - Nnn_all * E_nnn_K;
    rates_all = k0 * exp(-Ed_all / T_current);
    rates_all(occ_flat == 0) = 0;   % zero rate for empty sites
    total_rate = sum(rates_all);
    if total_rate == 0, break; end
    dt = 1 / total_rate;
    t = t + dt;
    T_current = 300 + Beta * t;
    if T_current > 600, break; end
    r = rand() * total_rate;
    cum_rates = cumsum(rates_all);
    idx = find(cum_rates >= r, 1);
    occ_flat(idx) = 0;
    occ_count = occ_count - 1;
    event = event + 1;
    T_record(event) = T_current;
    rate_record(event) = total_rate / N_sites;
    % Relaxation
    num_relax_attempts = 10 * occ_count;
    for attempt = 1:num_relax_attempts
        if occ_count == 0, break; end
        occ_indices = find(occ_flat == 1);
        p = occ_indices(randi(length(occ_indices)));
        q = neighbors(p, randi(6));
        if occ_flat(q) == 1, continue; end
        Nn_p  = sum(occ_flat(neighbors(p,:)));
        Nnn_p = sum(occ_flat(nextneighbors(p,:)));
        Ed_p = Ed0_K - Nn_p * E_nn_K - Nnn_p * E_nnn_K;
        occ_flat(p) = 0;
        Nn_q  = sum(occ_flat(neighbors(q,:)));
        Nnn_q = sum(occ_flat(nextneighbors(q,:)));
        Ed_q = Ed0_K - Nn_q * E_nn_K - Nnn_q * E_nnn_K;
        if Ed_q >= Ed_p || rand() < exp((Ed_q - Ed_p) / T_current)
            occ_flat(q) = 1;
        else
            occ_flat(p) = 1;
        end
    end
end
T_record = T_record(1:event);
rate_record = rate_record(1:event);
plot(T_record, rate_record, 'LineWidth', 1.5, 'Color', colors(1,:), 'DisplayName', sprintf("\\theta_0 = %.1f", theta0));
hold off;
xlabel('Temperature (K)');
ylabel('Desorption rate (per site s^{-1})');
legend('Location','best');
xlim([300 550]);


