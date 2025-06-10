function [best_params, fit_error] = fit_polyani_wigner(exp_time, exp_rate, init_params, beta, init_tmp, max_tmp, N_0)
% FIT_POLYANI_WIGNER Fits TDS data using the Polanyi-Wigner equation with Nelder-Mead optimization
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

% Define the objective function (sum of squared errors)
objective = @(params) calculate_error(params, exp_time, exp_rate, beta, init_tmp, max_tmp, N_0);

% Set up optimization options
options = optimset('Display', 'iter', ...
                  'MaxIter', 1000, ...
                  'MaxFunEvals', 2000, ...
                  'TolFun', 1e-6, ...
                  'TolX', 1e-3);

% Run the optimization using Nelder-Mead method
[best_params, fit_error] = fminsearch(objective, init_params, options);

% Plot the results
plot_results(exp_time, exp_rate, best_params, beta, init_tmp, max_tmp, N_0);
end

function error = calculate_error(params, exp_time, exp_rate, beta, init_tmp, max_tmp, N_0)
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
        
        % Calculate error (sum of squared differences)
        error = sum((model_rate_interp - exp_rate_valid).^2);
        
        % Ensure error is a scalar
        if ~isscalar(error)
            error = sum(error);
        end
        
    catch ME
        % If any error occurs, return a large error value
        fprintf('Error in objective function: %s\n', ME.message);
        error = 1e10;
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
    title('TDS Fit Results (Nelder-Mead)');
    legend('show');
    
    % Display parameters
    fprintf('\nFitted Parameters:\n');
    fprintf('Ea = %.14f kcal/mol\n', Ea/4184);
    fprintf('gam = %.14f kcal/mol\n', gam/4184);
    fprintf('Preexponential = %.14e s^-1\n', preexponent);
end