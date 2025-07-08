function [best_params, fit_error] = fit_polyani_wigner_2_tripletest(exp_data, init_params, beta, init_tmp, max_tmp)
% FIT_POLYANI_WIGNER_2_TRIPLETEST Fits multiple TDS spectra using the same kinetic parameters for polyani_wigner_2
% Inputs:
%   exp_data: Cell array of experimental data, each cell contains {time, rate, N_0}
%   init_params: Initial guess for [Ea, w_1, w_2, preexponent, T_c]
%   beta: Heating rate (K/s)
%   init_tmp: Initial temperature (K)
%   max_tmp: Maximum temperature (K)
% Outputs:
%   best_params: Best fit parameters [Ea, w_1, w_2, preexponent, T_c]
%   fit_error: Final fitting error (sum across all spectra)

objective = @(params) calculate_error_multiple_2(params, exp_data, beta, init_tmp, max_tmp);

options = optimoptions('fmincon', 'Display', 'iter', ...
                      'MaxIterations', 300, ...
                      'MaxFunctionEvaluations', 600, ...
                      'OptimalityTolerance', 1e-8, ...
                      'StepTolerance', 1e-6, ...
                      'FunctionTolerance', 1e-8, ... 
                      'Algorithm', 'interior-point');

lb = [140 * 1000, 0 * 1000, -1e6, 0, 0];    % Lower bounds for [Ea, w_1, w_2, preexponent, T_c]
ub = [180* 1000 , 60 * 1000 , 1e6, inf, inf];     % Upper bounds

[best_params, fit_error] = fmincon(objective, init_params, [], [], [], [], lb, ub, [], options);

plot_results_multiple_2(exp_data, best_params, beta, init_tmp, max_tmp);
end

function error = calculate_error_multiple_2(params, exp_data, beta, init_tmp, max_tmp)
    try
        Ea = params(1);
        w_1 = params(2);
        w_2 = params(3);
        preexponent = params(4);
        T_c = params(5);
        total_error = 0;
        for i = 1:length(exp_data)
            exp_time = exp_data{i}{1};
            exp_rate = exp_data{i}{2};
            N_0 = exp_data{i}{3};
            [model_time, ~, model_rate, ~] = polyani_wigner_niemant_2(beta, init_tmp, Ea, w_1, w_2, preexponent, N_0, max_tmp, T_c);
            exp_time = exp_time(:);
            exp_rate = exp_rate(:);
            model_time = model_time(:);
            model_rate = model_rate(:);
            model_rate_interp = interp1(model_time, model_rate, exp_time, 'linear', 'extrap');
            valid_idx = ~isnan(model_rate_interp);
            model_rate_interp = model_rate_interp(valid_idx);
            exp_rate_valid = exp_rate(valid_idx);
            spectrum_error = sum((model_rate_interp - exp_rate_valid).^2);
            total_error = total_error + spectrum_error;
        end
        error = total_error;
    catch ME
        fprintf('Error in objective function: %s\n', ME.message);
        error = 1e10;
    end
end

function plot_results_multiple_2(exp_data, params, beta, init_tmp, max_tmp)
    Ea = params(1);
    w_1 = params(2);
    w_2 = params(3);
    preexponent = params(4);
    T_c = params(5);
    figure; clf;
    colors = {'ko', 'bs', 'rd'};
    for i = 1:length(exp_data)
        exp_time = exp_data{i}{1};
        exp_rate = exp_data{i}{2};
        N_0 = exp_data{i}{3};
        [model_time, ~, model_rate, ~] = polyani_wigner_niemant_2(beta, init_tmp, Ea, w_1, w_2, preexponent, N_0, max_tmp, T_c);
        plot(exp_time, exp_rate, colors{i}, 'DisplayName', sprintf('Experimental N_0=%.3f', N_0));
        hold on;
        plot(model_time, model_rate, colors{i}(1), 'DisplayName', sprintf('Fitted N_0=%.3f', N_0));
    end
    xlabel('Time (s)');
    ylabel('Desorption Rate');
    title('Polyani Wigner Niemantsverdriet (2nd Order) Multi-Spectrum Fit');
    legend('show');
    fprintf('\nFitted Parameters:\n');
    fprintf('Ea = %.14f kJ/mol\n', Ea/1000);
    fprintf('w_1 = %.14f J/mol\n', w_1);
    fprintf('w_2 = %.14f J/mol\n', w_2);
    fprintf('Preexponential = %.14e s^-1\n', preexponent);
    fprintf('T_c = %.14f K \n', T_c);
    fprintf('\nCoverages used:\n');
    for i = 1:length(exp_data)
        fprintf('Spectrum %d: N_0 = %.3f\n', i, exp_data{i}{3});
    end
end
