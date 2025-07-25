function [best_params, fit_error] = fit_polyani_wigner_comp1(exp_time, exp_rate, init_params, beta, init_tmp, max_tmp, N_0)
% FIT_POLYANI_WIGNER_COMP1 Fits TDS data using the Polanyi-Wigner equation with compensation effect
% Inputs:
%   exp_time: Experimental time data points
%   exp_rate: Experimental desorption rate data points
%   init_params: Initial guess for [Ea, gam, w, preexponent]
%   beta: Heating rate (K/s)
%   init_tmp: Initial temperature (K)
%   max_tmp: Maximum temperature (K)
%   N_0: Initial coverage
% Outputs:
%   best_params: Best fit parameters [Ea, gam, w, preexponent]
%   fit_error: Final fitting error

% Define the objective function (sum of squared errors)
objective = @(params) calculate_error(params, exp_time, exp_rate, beta, init_tmp, max_tmp, N_0);

% Set up optimization options
options = optimset('Display', 'iter', ...
                  'MaxIter', 10000, ...
                  'MaxFunEvals', 20000, ...
                  'TolFun', 1e-10, ...
                  'TolX', 1e-8);

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
        w = params(3);
        preexponent = params(4);
        
        % Get model prediction
        [model_time, ~, model_rate, ~] = polyani_wigner_comp1(beta, init_tmp, Ea, gam, w, preexponent, N_0, max_tmp);
        
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
    w = params(3);
    preexponent = params(4);
    
    % Get model prediction
    [model_time, ~, model_rate, ~] = polyani_wigner_comp1(beta, init_tmp, Ea, gam, w, preexponent, N_0, max_tmp);
    
    % Create figure
    figure; clf;
    plot(exp_time, exp_rate, 'ko', 'DisplayName', 'Experimental');
    hold on;
    plot(model_time, model_rate, 'r-', 'DisplayName', 'Fitted');
    xlabel('Time (s)');
    ylabel('Desorption Rate');
    title('TDS Fit, log10(v) and Ea linear with coverage');
    legend('show');
    
    % Display parameters
    fprintf('\nFitted Parameters:\n');
    fprintf('Ea = %.14f kcal/mol\n', Ea/4184);
    fprintf('gam = %.14f J/coverage \n', gam);
    fprintf('w = %.14f \n J/mol', w);
    fprintf('\n Preexponential = %.14e s^-1\n', preexponent);
end 