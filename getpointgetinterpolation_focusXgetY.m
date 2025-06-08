function [x_actual, y_actual] = getpointgetinterpolation_focusXgetY(x_span, y_span, x_point, tolerance)
%'focusYgetX' takes in a dependent variable for the center of
%interpolation. This takes an independent variable instead. 

%First on 'erley80_test'. 6/7
desired_x = x_point; %Code definition redunancy to remember original code
[~,ind_center] = min(abs(x_span - desired_x));
%tolerance = 0.5;
avg_spacing = abs(x_span(1) - x_span(end)) / size(x_span,2);

if avg_spacing > tolerance
    tolerance = avg_spacing * 2;
end

%fprintf('%.4f\n', avg_spacing);
ind_range = round((tolerance/avg_spacing) / 2);
fprintf('%.4f\n', ind_range);
x_span_focus = x_span(ind_center - ind_range : ind_center + ind_range);
y_span_focus =  y_span(ind_center - ind_range :  ind_center + ind_range);


figure(7); clf;
hold on;
plot(x_span(ind_center - ind_range : ind_center + ind_range) , y_span(ind_center - ind_range :  ind_center + ind_range), 'b--','LineWidth',3);
hold off;


x_span_focus_interp = linspace(min(x_span_focus), max(x_span_focus),2000); 
%Raise the number at the end for more accuracy. Not sure what's too high to become unreasonably slow. 
%Looks like returns diminish hard past ~5000 points. ~7-fold improvement
%from 1000 to 2000, but only ~2/3 fold improvement from 2000 to 5000. Kept
%at 2000 as of 6/7. 

y_span_focus_interp = interp1(x_span_focus,y_span_focus, x_span_focus_interp, 'pchip');

figure(7);
hold on;
plot(x_span_focus_interp, y_span_focus_interp, 'r','LineWidth',2);
plot(x_span, y_span, 'k--')
hold off;

[~,ind_actual] = min(abs(x_span_focus_interp - desired_x));
x_actual = x_span_focus_interp(ind_actual);
y_actual = y_span_focus_interp(ind_actual);

fprintf('%.4e\n', x_actual - x_point);
end