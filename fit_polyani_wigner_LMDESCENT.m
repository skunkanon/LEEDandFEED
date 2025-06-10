function [best_params, fit_error] = fit_polyani_wigner(exp_time, exp_rate, init_params, beta, init_tmp, max_tmp, N_0)
% FIT_POLYANI_WIGNER Fits TDS data using the Polanyi-Wigner equation
% Inputs:
%   exp_time: Experimental time data points
%   exp_rate: Experimental desorption rate data points
%   init_params: Initial guess for [Ea, gam, preexponent]
%   beta: Heating rate (K/s)
%   init_tmp: Initial temperature (K)
%   max_tmp: Maximum temperature (K)
%   N_0: Initial coverage
% Outputs:
%   best_params: Best fit parameters [Ea, gam, preexponent]
%   fit_error: Final fitting error

% Set up optimization options for Levenberg-Marquardt
options = optimoptions('lsqnonlin', ...
                      'Display', 'iter', ...
                      'MaxIterations', 100, ...
                      'MaxFunctionEvaluations', 200, ...
                      'FunctionTolerance', 1e-6, ...
                      'StepTolerance', 1e-3);

% Define residuals function for Levenberg-Marquardt
residuals = @(params) calculate_residuals(params, exp_time, exp_rate, beta, init_tmp, max_tmp, N_0);

% Run the optimization
[best_params, fit_error] = lsqnonlin(residuals, init_params, [], [], options);
fit_error = fit_error^2;  % Convert to sum of squares

% Plot the results
plot_results(exp_time, exp_rate, best_params, beta, init_tmp, max_tmp, N_0);
end

function residuals = calculate_residuals(params, exp_time, exp_rate, beta, init_tmp, max_tmp, N_0)
    try
        % Extract parameters
        Ea = params(1);
        gam = params(2);
        preexponent = params(3);
        
        % Get model prediction
        [model_time, ~, model_rate, ~] = polyani_wigner(beta, init_tmp, Ea, gam, preexponent, N_0, max_tmp);
        
        % Ensure all inputs are column vectors
        exp_time = exp_time(:);
        exp_rate = exp_rate(:);
        model_time = model_time(:);
        model_rate = model_rate(:);
        
        % Interpolate model results to match experimental time points
        model_rate_interp = interp1(model_time, model_rate, exp_time, 'linear', 'extrap');
        
        % Remove any NaN values that might occur during interpolation
        valid_idx = ~isnan(model_rate_interp);
        model_rate_interp = model_rate_interp(valid_idx);
        exp_rate_valid = exp_rate(valid_idx);
        
        % Calculate residuals (differences)
        residuals = model_rate_interp - exp_rate_valid;
        
    catch ME
        % If any error occurs, return large residuals
        fprintf('Error in residuals function: %s\n', ME.message);
        residuals = 1e5 * ones(size(exp_rate));
    end
end

function plot_results(exp_time, exp_rate, params, beta, init_tmp, max_tmp, N_0)
    % Extract parameters
    Ea = params(1);
    gam = params(2);
    preexponent = params(3);
    
    % Get model prediction
    [model_time, ~, model_rate, ~] = polyani_wigner(beta, init_tmp, Ea, gam, preexponent, N_0, max_tmp);
    
    % Create figure
    figure;
    plot(exp_time, exp_rate, 'ko', 'DisplayName', 'Experimental');
    hold on;
    plot(model_time, model_rate, 'r-', 'DisplayName', 'Fitted');
    xlabel('Time (s)');
    ylabel('Desorption Rate');
    title('TDS Fit Results (Levenberg-Marquardt)');
    legend('show');
    
    % Display parameters
    fprintf('\nFitted Parameters:\n');
    fprintf('Ea = %.2f kJ/mol\n', Ea/1000);
    fprintf('gam = %.2f kJ/mol\n', gam/1000);
    fprintf('Preexponential = %.2e s^-1\n', preexponent);
end