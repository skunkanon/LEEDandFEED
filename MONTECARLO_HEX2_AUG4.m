function [TPD_data, energyStats, cov_events, Ed_events] = MONTECARLO_HEX2_AUG4(L, theta0, Ed0, eps_nn, eps_nnn, beta, k0, targetCoverage)
% MONTECARLO_HEX_JUL25 
% Monte Carlo simulation of first-order TPD on a hexagonal lattice with nearest-neighbor 
% and next-nearest-neighbor interactions (Meng & Weinberg, J. Chem. Phys. 100, 1994).
%
% Inputs:
%   L             - Lattice size (L x L sites, hexagonal coordination via periodic BCs)
%   theta0        - Initial fractional coverage (0 < theta0 <= 1)
%   Ed0           - Desorption energy for an isolated molecule (e.g., in eV)
%   eps_nn        - Interaction energy per nearest neighbor (repulsive if positive, eV)
%   eps_nnn       - Interaction energy per next-nearest neighbor (eV)
%   beta          - Heating rate (K/s) 
%   k0            - Pre-exponential factor for desorption (s^-1)
%   targetCoverage (optional) - Coverage at which to report average Ed (default 0.5)
%
% Outputs:
%   TPD_data    - Matrix with two columns [T, theta] giving surface coverage theta vs. temperature T during desorption
%   energyStats - (Optional) 1x2 vector [Ed_avg, Ed_std] = average desorption energy and standard deviation for desorption events around targetCoverage.
%
% The code preserves the original output (TPD_data) and adds a second optional output (energyStats) for the energy statistics.
%
% Reference: B. Meng and W.H. Weinberg, J. Chem. Phys. 100, 5280 (1994).

    if nargin < 8
        targetCoverage = 0.5;  % default target coverage for energy analysis
    end
    % Set default values for optional inputs if not provided
    if nargin < 7 || isempty(k0)
        k0 = 1e13;  % s^-1, typical first-order desorption pre-factor
    end
    if nargin < 6 || isempty(beta)
        beta = 5;   % K/s, default linear heating rate
    end

    % Constants and initialization
    kB = 8.617333e-5;  % Boltzmann constant in eV/K (for energies in eV and T in K)
    N_sites = L * L;                      % total number of adsorption sites
    N0 = round(theta0 * N_sites);         % number of molecules initially adsorbed
    % Ensure N0 does not exceed lattice
    N0 = min(N0, N_sites);
    
    % Initialize lattice occupancy (hexagonal lattice with periodic boundaries)
    occ = false(L);           % occupancy matrix (false = empty, true = occupied)
    % Randomly populate N0 sites with adsorbates
    randIndices = randperm(N_sites, N0);
    occ(randIndices) = true;
    currentN = N0;            % current number of adsorbed molecules
    
    % If needed, relax initial configuration via surface diffusion to equilibrium
    % (Assuming rapid surface diffusion as in Meng & Weinberg's algorithm:contentReference[oaicite:5]{index=5}:contentReference[oaicite:6]{index=6})
    % *** Existing code for diffusion relaxation would go here, if implemented ***
    
    % Preallocate history arrays for performance
    maxEvents = N0; 
    T_history    = zeros(maxEvents+1, 1);  % Temperature at each step (including initial)
    cov_history  = zeros(maxEvents+1, 1);  % Coverage at each step (including initial)
    Ed_events    = zeros(maxEvents, 1);    % Activation energy of each desorption event
    cov_events   = zeros(maxEvents, 1);    % Coverage just before each desorption event
    
    % Initial state
    time = 0;                % simulation time (s)
    T_history(1)   = 0;      % starting temperature (assume ~0 K for simulation start)
    cov_history(1) = currentN / N_sites;  % initial fractional coverage (theta0)
    T_current = 0;           % initial temperature

    % Initialization debug
    fprintf('Starting KMC TPD simulation on %dx%d hexagonal lattice...\n', L, L);
    fprintf('Initial coverage: %.2f, Initial adsorbates: %d\n', theta0, N0);
    fprintf('Parameters: Ed0=%.3f eV, eps_nn=%.4f eV, eps_nnn=%.4f eV, beta=%.2f K/s, k0=%.2e s^-1\n', Ed0, eps_nn, eps_nnn, beta, k0);
    
    % Calculate initial desorption rates for all molecules
    [r_max, ~] = computeDesorptionRates(T_current);   % pass current temperature
    
    % Main KMC loop for desorption events
    eventCount = 0;
    fprintf('Beginning desorption loop...\n');
    progressCheckpoints = round(linspace(1, N0, 11)); % 0%,10%,...,100%
    while currentN > 0
        % Choose a random occupied site (candidate for desorption)
        idx = chooseRandomOccupiedSite();  % (function returns linear index of a random occupied site)
        [nn_count, nnn_count] = countNeighbors(idx);           % get local neighbor counts for this site
        Ed = Ed0 - eps_nn*nn_count - eps_nnn*nnn_count;        % effective desorption barrier for this site
        r_idx = k0 * exp(-Ed / (kB * T_current)); % desorption rate for this site at current T
        P = r_idx / r_max;  % acceptance probability for this site
        
        if rand() < P
            % Desorption event occurs
            eventCount = eventCount + 1;
            % Log the coverage and Ed for this desorption event
            cov_events(eventCount) = currentN / N_sites;  % coverage *before* desorption
            Ed_events(eventCount)  = Ed;
            % Remove the molecule from the lattice
            occ(idx) = false;
            currentN = currentN - 1;
            % Update simulation time by time increment Δt = 1 / (sum of all rates)
            totalRate = computeTotalRate(T_current);    % pass current temperature
            if totalRate > 0
                dt = 1 / totalRate;
            else
                dt = 1e6;  % (arbitrarily large time step)
            end
            time = time + dt;
            % Update temperature based on heating rate: T = T0 + β * time
            T_current = beta * time; 
            % Record the new state after this event
            T_history(eventCount+1)   = T_current;
            cov_history(eventCount+1) = currentN / N_sites;
            % Recalculate r_max after removal (neighbors' environments changed)
            [r_max, ~] = computeDesorptionRates(T_current);  % pass current temperature
            % Optionally, relax the surface via diffusion to equilibrium after desorption
            % *** Existing diffusion equilibration code would be called here, if implemented ***
            % Progress debug
            if any(eventCount == progressCheckpoints)
                fprintf('  Desorption progress: %3.0f%% (Event %d/%d, Coverage: %.3f, T: %.1f K)\n', ...
                    100*eventCount/N0, eventCount, N0, cov_events(eventCount), T_current);
            end
        end
        % If the event did not occur (rand >= P), loop and try another random site without time advancement
    end
    fprintf('Desorption complete. Total events: %d\n', eventCount);
    % Trim history arrays to actual number of events + 1
    T_history   = T_history(1:eventCount+1);
    cov_history = cov_history(1:eventCount+1);
    Ed_events   = Ed_events(1:eventCount);
    cov_events  = cov_events(1:eventCount);
    
    % Prepare primary output: coverage vs temperature data
    TPD_data = [T_history, cov_history];
    
    % --- Begin energy statistics aggregation (new) ---
    % Compute average and std of Ed for events occurring at targetCoverage (±0.005)
    target = targetCoverage;
    width  = 0.01;  % coverage window width for averaging (e.g., ±0.005 around target)
    cov_min = target - width/2;
    cov_max = target + width/2;
    % Find all desorption events where coverage before desorption was in [cov_min, cov_max]
    idx_events = find(cov_events >= cov_min & cov_events <= cov_max);
    if isempty(idx_events)
        % If no events in this range (e.g., coverage changes in bigger steps), include nearest event
        [~, nearestIdx] = min(abs(cov_events - target));
        idx_events = nearestIdx;
    end
    selected_Eds = Ed_events(idx_events);
    % Calculate statistics
    Ed_avg = mean(selected_Eds);
    Ed_std = std(selected_Eds, 0);  % standard deviation (population std since sample is the actual events)
    % Set optional output
    energyStats = [Ed_avg, Ed_std];
    % --- End energy statistics aggregation ---
    % Add outputs for event data
    % cov_events and Ed_events are already trimmed to eventCount
    % (no further action needed)

    % Helper function to compute desorption rates and maximum rate
    function [r_max, rates] = computeDesorptionRates(T)
        rates = [];  %#ok<NASGU>  % (unused output in this simplified implementation)
        r_max = 0;
        % Loop over all occupied sites
        occ_linear_indices = find(occ);
        for ii = 1:length(occ_linear_indices)
            site = occ_linear_indices(ii);
            [nn, nnn] = countNeighbors(site);
            E_site = Ed0 - eps_nn*nn - eps_nnn*nnn;
            r_site = k0 * exp(-E_site / (kB * T));
            if r_site > r_max
                r_max = r_site;
            end
        end
    end

    % Helper function to compute total desorption rate (sum over all occupied sites)
    function totalR = computeTotalRate(T)
        totalR = 0;
        occ_linear_indices = find(occ);
        for ii = 1:length(occ_linear_indices)
            site = occ_linear_indices(ii);
            [nn, nnn] = countNeighbors(site);
            E_site = Ed0 - eps_nn*nn - eps_nnn*nnn;
            r_site = k0 * exp(-E_site / (kB * T));
            totalR = totalR + r_site;
        end
    end

    % Helper function to choose a random occupied site (returns linear index)
    function idx = chooseRandomOccupiedSite()
        occ_indices = find(occ);
        randIdx = randi(length(occ_indices));
        idx = occ_indices(randIdx);
    end

    % Helper function to count nearest and next-nearest neighbors for a given site index
    function [nn_count, nnn_count] = countNeighbors(siteIndex)
        % Convert linear index to 2D coordinates (row, col)
        [row, col] = ind2sub([L, L], siteIndex);
        % Define neighbor displacements for hexagonal lattice (axial coordinates representation)
        % Nearest neighbors (6 offsets):
        nn_offsets = [ 0  1;  0 -1;  1  0; -1  0;  1 -1; -1  1 ];
        % Next-nearest neighbors (6 offsets, second ring):
        nnn_offsets = [ 1  1; -1 -1; 2  0; -2  0;  1 -2; -1  2 ];
        % Count occupied neighbors applying periodic boundary conditions
        nn_count = 0;
        nnn_count = 0;
        for k = 1:size(nn_offsets,1)
            rr = mod(row + nn_offsets(k,1) - 1, L) + 1;
            cc = mod(col + nn_offsets(k,2) - 1, L) + 1;
            if occ(rr, cc)
                nn_count = nn_count + 1;
            end
        end
        for k = 1:size(nnn_offsets,1)
            rr = mod(row + nnn_offsets(k,1) - 1, L) + 1;
            cc = mod(col + nnn_offsets(k,2) - 1, L) + 1;
            if occ(rr, cc)
                nnn_count = nnn_count + 1;
            end
        end
    end
end