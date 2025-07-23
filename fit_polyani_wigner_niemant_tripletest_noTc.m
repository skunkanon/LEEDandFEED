function [best_params, fit_error, hessian, cov_matrix, stderr] = fit_polyani_wigner_niemant_tripletest_noTc(exp_data, init_params, Tc_fixed, beta, init_tmp, max_tmp)
% FIT_POLYANI_WIGNER_NIEMANT_TRIPLETEST_NOTC Fits multiple TDS spectra using the same kinetic parameters with weighted least squares and parameter scaling, with T_c fixed
% Inputs:
%   exp_data: Cell array of experimental data, each cell contains {time, rate, N_0, uncertainties}
%   init_params: Initial guess for [Ea, w, preexponent]
%   Tc_fixed: Fixed value for T_c (not optimized)
%   beta: Heating rate (K/s)
%   init_tmp: Initial temperature (K)
%   max_tmp: Maximum temperature (K)
% Outputs:
%   best_params: Best fit parameters [Ea, w, preexponent] (T_c fixed)
%   fit_error: Final fitting error (sum across all spectra)
%   hessian: Hessian matrix at the solution (3x3)
%   cov_matrix: Covariance matrix for the three parameters
%   stderr: Standard errors for the three parameters

% Define scaling factors to normalize parameters to similar magnitudes
scale_factors = [1e5, 1e4, 1e7]; % Scale to roughly 1-10 range
scaled_init_params = init_params ./ scale_factors;

% Define the objective function (sum of weighted squared errors across all spectra)
objective = @(scaled_params) calculate_error_multiple_weighted_scaled_noTc(scaled_params, exp_data, Tc_fixed, beta, init_tmp, max_tmp, scale_factors);

% Set up optimization options for fmincon
options = optimoptions('fmincon', 'Display', 'iter', ...
                      'MaxIterations', 400, ...
                      'MaxFunctionEvaluations', 800, ...
                      'OptimalityTolerance', 1e-12, ...
                      'StepTolerance', 1e-12, ...
                      'FunctionTolerance',1e-12, ... 
                      'Algorithm', 'interior-point');

% Set bounds: [Ea, w, preexponent] (scaled)
lb = [140 * 1000, 0, 0] ./ scale_factors;    % Lower bounds
ub = [280 * 1000 , 100*1000, inf] ./ scale_factors;     % Upper bounds

[scaled_best_params, fit_error, ~, ~, ~, ~, scaled_hessian] = fmincon(objective, scaled_init_params, [], [], [], [], lb, ub, [], options);

% Rescale parameters back to original units
best_params = scaled_best_params .* scale_factors;

% Calculate number of data points (for all spectra)
n_data = 0;
for i = 1:length(exp_data)
    n_data = n_data + length(exp_data{i}{2});
end
n_params = 3; % Only three parameters
residual_variance = fit_error / (n_data - n_params);

% Rescale Hessian back to original parameter space
if ~isempty(scaled_hessian) && all(size(scaled_hessian) == [3 3]) && all(~isnan(scaled_hessian(:)))
    % Create scaling matrix
    scale_matrix = diag(scale_factors);
    hessian = scale_matrix * scaled_hessian * scale_matrix;
    
    % Calculate covariance and standard errors
    cov_matrix = residual_variance * inv(hessian);
    stderr = sqrt(diag(cov_matrix));
else
    hessian = nan(3);
    cov_matrix = nan(3);
    stderr = nan(3,1);
end

% Plot the results for all spectra if fit succeeded
if all(~isnan(best_params))
    plot_results_multiple_noTc(exp_data, best_params, Tc_fixed, beta, init_tmp, max_tmp);
end

% Print standard errors
fprintf('\nStandard errors (Ea, w, preexponent):\n');
disp(stderr');

% Print condition number of Hessian
if ~isempty(hessian) && all(size(hessian) == [3 3]) && all(~isnan(hessian(:)))
    eigenvals = eig(hessian);
    condition_number = max(eigenvals) / min(eigenvals);
    fprintf('\nHessian condition number: %.2e\n', condition_number);
    fprintf('Eigenvalues: [%.2e, %.2e, %.2e]\n', eigenvals(1), eigenvals(2), eigenvals(3));
end
end

function error = calculate_error_multiple_weighted_scaled_noTc(scaled_params, exp_data, Tc_fixed, beta, init_tmp, max_tmp, scale_factors)
    % Rescale parameters back to original units for calculation
    params = scaled_params .* scale_factors;
    
    % Extract parameters (3 parameters + fixed T_c)
    Ea = params(1);
    w = params(2);
    preexponent = params(3);
    T_c = Tc_fixed;
    
    total_error = 0;
    
    % Loop through all experimental spectra
    for i = 1:length(exp_data)
        exp_time = exp_data{i}{1};
        exp_rate = exp_data{i}{2};
        exp_uncertainty = exp_data{i}{4}; % Uncertainties from fourth element
        N_0 = exp_data{i}{3};
        
        % Get model prediction for this coverage
        [model_time, ~, model_rate, ~] = polyani_wigner_niemant(beta, init_tmp, Ea, w, preexponent, N_0, max_tmp, T_c);
        
        % Ensure all inputs are column vectors
        exp_time = exp_time(:);
        exp_rate = exp_rate(:);
        exp_uncertainty = exp_uncertainty(:);
        model_time = model_time(:);
        model_rate = model_rate(:);
        
        % Interpolate model results to match experimental time points
        model_rate_interp = interp1(model_time, model_rate, exp_time, 'linear', 'extrap');
        
        % Remove any NaN values that might occur during interpolation
        valid_idx = ~isnan(model_rate_interp);
        model_rate_interp = model_rate_interp(valid_idx);
        exp_rate_valid = exp_rate(valid_idx);
        exp_uncertainty_valid = exp_uncertainty(valid_idx);
        
        % Weighted least squares
        weights = 1 ./ (exp_uncertainty_valid.^2); % Inverse variance weighting
        residuals = (model_rate_interp - exp_rate_valid) .* sqrt(weights);
        
        % Calculate weighted error for this spectrum
        spectrum_error = sum(residuals.^2);
        total_error = total_error + spectrum_error;
    end
    
    error = total_error;
end

function plot_results_multiple_noTc(exp_data, params, Tc_fixed, beta, init_tmp, max_tmp)
    % Extract parameters (3 parameters + fixed T_c)
    Ea = params(1);
    w = params(2);
    preexponent = params(3);
    T_c = Tc_fixed;
    
    % Create figure
    figure; clf;
    cmap = lines(length(exp_data)); % Different markers for each spectrum
    
    % Plot each spectrum
    for i = 1:length(exp_data)
        exp_time = exp_data{i}{1};
        exp_rate = exp_data{i}{2};
        N_0 = exp_data{i}{3};
        
        % Get model prediction
        [model_time, ~, model_rate, ~] = polyani_wigner_niemant(beta, init_tmp, Ea, w, preexponent, N_0, max_tmp, T_c);
        
        % Plot experimental data
        plot(exp_time, exp_rate, 'Color', cmap(i,:), 'DisplayName', sprintf('Experimental N_0=%.3f', N_0));
        hold on;
        
        % Plot fitted data
        plot(model_time, model_rate, 'Color', cmap(i,:), 'DisplayName', sprintf('Fitted N_0=%.3f', N_0));
    end
    
    xlabel('Time (s)');
    ylabel('Desorption Rate');
    title('Niemantsverdriet Multi-Spectrum Fit (Weighted, T_c fixed)');
    legend('show');
    
    % Display parameters
    fprintf('\nFitted Parameters:\n');
    fprintf('Ea = %.14f kJ/mol\n', Ea/1000);
    fprintf('w = %.14f J/mol\n', w);
    fprintf('Preexponential = %.14e s^-1\n', preexponent);
    fprintf('T_c = %.14f K (fixed)\n', T_c);
    
    % Display coverages used
    fprintf('\nCoverages used:\n');
    for i = 1:length(exp_data)
        fprintf('Spectrum %d: N_0 = %.3f\n', i, exp_data{i}{3});
    end

end
