function [] = polyani_plot(x_spans, y_spans, name_span, plot_title)

figure;
hold on;
for i = 1:length(x_spans) % make curves
    plot(x_spans{i},y_spans{i}, 'DisplayName', string(name_span(i)));  

end


for i = 1:length(y_spans)
    [M, I] = max(y_spans{i});
    plot(x_spans{i}(I), M, 'o', 'MarkerSize',5 ,'DisplayName', name_span(i));


end

hold off;





xlabel('Temperature (K)');
ylabel('Rate = -dN/dt');
title(plot_title);
legend show;
grid on;

end

