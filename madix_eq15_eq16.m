%% EQ 15 and 16 IMPLEMENTATION - FIXED VERSION
% Generate Desorption Spectra 
fig5_beta = linspace(0.2, 35, 18); 
[e16_x_5, e16_y_5, e16_t_5] = deal(cell(1, length(fig5_beta)));
[e16_Tp_5, tempindex, int_list] = deal(zeros(1, length(fig5_beta)));

% Data generation remains the same
fprintf("\nFIGURE 5\n");
for i = 1:length(fig5_beta)
    [e16_t_5{i}, e16_x_5{i}, e16_y_5{i}] = madix_sim_eq16(fig5_beta(i));
    [~, tempindex(i)] = max(e16_y_5{i});
    e16_Tp_5(i) = e16_x_5{i}(tempindex(i));
    fprintf("Peak Temp for B = %.4f: %.2f°C\n", fig5_beta(i), e16_Tp_5(i)-273);
end

%% Show Desorption Plots (Fixed)
figure(1); clf;  % Explicit figure number + clear existing
%set(gcf, 'Position', [100 100 800 600]);  % Set figure size
hold on;

% Use line instead of plot for better performance
colors = parula(length(fig5_beta));  % Better colormap
for i = 1:length(fig5_beta)
    line(e16_x_5{i}, e16_y_5{i}, 'Color', colors(i,:), 'LineWidth', 1.5,...
        'DisplayName', sprintf('β = %.1f K/s', fig5_beta(i)));
end

hold off;
xlabel('Temperature (K)');
ylabel('Desorption Rate (-dθ/dt)');
title('TPD Spectra: Variation with Heating Rate (β)');
legend('show', 'Location', 'eastoutside', 'NumColumns', 2);  % Better legend placement
grid on;
drawnow;  % Force immediate update

%% Equation 16 Plot (Fixed)
figure(2); clf;
%set(gcf, 'Position', [200 200 700 500]);


% Enhanced plot
scatter(x_e16, y_e16, 80, 'filled', 'MarkerEdgeColor', 'k',...
    'MarkerFaceColor', [0.2 0.6 0.9], 'DisplayName', 'Experimental');
hold on;
plot(x_e16, y_e16_fit, 'r--', 'LineWidth', 2, 'DisplayName', 'Linear Fit');

% Add equation text
eq_text = sprintf('y = %.2fx + %.2f\nE_a = %.1f kcal/mol',...
    e16_m, e16_b, (-e16_m*R)/4184);
text(0.6*max(x_e16), 0.8*max(y_e16), eq_text,...
    'FontSize', 12, 'BackgroundColor', 'w');

hold off;
xlabel('1/T_p (K^{-1})');
ylabel('ln(-dθ/dt|_{peak})');
title('Arrhenius Analysis (Equation 16)');
legend('show', 'Location', 'northwest');
grid on;
drawnow;

%% Equation 15 Plot (Fixed)
figure(3); clf;
%set(gcf, 'Position', [300 300 700 500]);

% ... [keep existing calculations] ...

% Enhanced visualization
scatter(tempx, tempy, 80, 'filled', 'MarkerEdgeColor', 'k',...
    'MarkerFaceColor', [0.9 0.4 0.3], 'DisplayName', 'Data');
hold on;
plot(tempx, fitted_y, 'b-', 'LineWidth', 2, 'DisplayName', 'Linear Fit');

% Add equation text
eq_text = sprintf('y = %.2fx + %.2f\nE_a = %.1f kcal/mol', m, b, Ea_calc);
text(0.6*max(tempx), 0.8*max(tempy), eq_text,...
    'FontSize', 12, 'BackgroundColor', 'w');

hold off;
xlabel('1/T_p (K^{-1})');
ylabel('ln(β/T_p^2)');
title('Pre-exponential Analysis (Equation 15)');
legend('show', 'Location', 'southeast');
grid on;
drawnow;