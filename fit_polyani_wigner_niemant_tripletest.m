function [best_params, fit_error] = fit_polyani_wigner_niemant_tripletest(exp_data, init_params, beta, init_tmp, max_tmp)
% FIT_POLYANI_WIGNER_NIEMANT_TRIPLETEST Fits multiple TDS spectra using the same kinetic parameters
% Inputs:
%   exp_data: Cell array of experimental data, each cell contains {time, rate, N_0}
%             exp_data{1} = {time1, rate1, N_0_1}
%             exp_data{2} = {time2, rate2, N_0_2}
%             exp_data{3} = {time3, rate3, N_0_3}
%   init_params: Initial guess for [Ea, w, preexponent, T_c]
%   beta: Heating rate (K/s)
%   init_tmp: Initial temperature (K)
%   max_tmp: Maximum temperature (K)
% Outputs:
%   best_params: Best fit parameters [Ea, w, preexponent, T_c]
%   fit_error: Final fitting error (sum across all spectra)

% Define the objective function (sum of squared errors across all spectra)
objective = @(params) calculate_error_multiple(params, exp_data, beta, init_tmp, max_tmp);

% Set up optimization options for fmincon
options = optimoptions('fmincon', 'Display', 'iter', ...
                      'MaxIterations', 300, ...
                      'MaxFunctionEvaluations', 600, ...
                      'OptimalityTolerance', 1e-16, ...
                      'StepTolerance', 1e-16, ...
                      'FunctionTolerance',1e-16, ... 
                      'Algorithm', 'interior-point');

% Set bounds: [Ea, w, preexponent, T_c]

%lb = [init_params(1) - 20*1000, -50 * 1000, 10^6.8, 0];    % Lower bounds
%ub = [init_params(1) + 20*1000, 0, 10^7.8, inf];     % Upper bounds

lb = [140 * 1000, 0, 0, 0];    % Lower bounds
ub = [inf , 100*1000, inf, inf];     % Upper bounds

% Run the optimization using fmincon with bounds
[best_params, fit_error] = fmincon(objective, init_params, [], [], [], [], lb, ub, [], options);

% Plot the results for all spectra
plot_results_multiple(exp_data, best_params, beta, init_tmp, max_tmp);
end

function error = calculate_error_multiple(params, exp_data, beta, init_tmp, max_tmp)
    try
        % Extract parameters
        Ea = params(1);
        w = params(2);
        preexponent = params(3);
        T_c = params(4);
        
        total_error = 0;
        
        % Loop through all experimental spectra
        for i = 1:length(exp_data)
            exp_time = exp_data{i}{1};
            exp_rate = exp_data{i}{2};
            N_0 = exp_data{i}{3};
            
            % Get model prediction for this coverage
            [model_time, ~, model_rate, ~] = polyani_wigner_niemant(beta, init_tmp, Ea, w, preexponent, N_0, max_tmp, T_c);
            
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
            
            % Calculate error for this spectrum
            spectrum_error = sum((model_rate_interp - exp_rate_valid).^2);
            total_error = total_error + spectrum_error;
        end
        
        error = total_error;
        
    catch ME
        % If any error occurs, return a large error value
        fprintf('Error in objective function: %s\n', ME.message);
        error = 1e10;
    end
end

function plot_results_multiple(exp_data, params, beta, init_tmp, max_tmp)
    % Extract parameters
    Ea = params(1);
    w = params(2);
    preexponent = params(3);
    T_c = params(4);
    
    % Create figure
    figure; clf;
    colors = {'ko', 'bs', 'rd'}; % Different markers for each spectrum
    
    % Plot each spectrum
    for i = 1:length(exp_data)
        exp_time = exp_data{i}{1};
        exp_rate = exp_data{i}{2};
        N_0 = exp_data{i}{3};
        
        % Get model prediction
        [model_time, ~, model_rate, ~] = polyani_wigner_niemant(beta, init_tmp, Ea, w, preexponent, N_0, max_tmp, T_c);
        
        % Plot experimental data
        plot(exp_time, exp_rate, colors{i}, 'DisplayName', sprintf('Experimental N_0=%.3f', N_0));
        hold on;
        
        % Plot fitted data
        plot(model_time, model_rate, colors{i}(1), 'DisplayName', sprintf('Fitted N_0=%.3f', N_0));
    end
    
    xlabel('Time (s)');
    ylabel('Desorption Rate');
    title('Niemantsverdriet Multi-Spectrum Fit');
    legend('show');
    
    % Display parameters
    fprintf('\nFitted Parameters:\n');
    fprintf('Ea = %.14f kJ/mol\n', Ea/1000);
    fprintf('w = %.14f J/mol\n', w);
    fprintf('Preexponential = %.14e s^-1\n', preexponent);
    fprintf('T_c = %.14f K \n', T_c);
    
    % Display coverages used
    fprintf('\nCoverages used:\n');
    for i = 1:length(exp_data)
        fprintf('Spectrum %d: N_0 = %.3f\n', i, exp_data{i}{3});
    end

    
end