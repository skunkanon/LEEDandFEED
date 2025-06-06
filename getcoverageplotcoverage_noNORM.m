function [coverage_span_actual] = getcoverageplotcoverage_noNORM(time_span, signal_span, N0, color)
%6/5 - Commenting out normalization line.

coverage_span_raw = zeros(1,length(time_span));

for i = 1:length(time_span) - 1 
    coverage_span_raw(i) = trapz(time_span(i:end) , signal_span(i:end)); %gets the area underneath the curve past the signal (new_y) at index 'i'
end


fprintf('NEW INSTance \n');

[max_coverage_span_raw, ~] = max(coverage_span_raw); %gets maximum of areas for below normalization 
coverage_span_actual = coverage_span_raw;
%coverage_span_actual = coverage_span_raw * (N0 / max_coverage_span_raw); %trapz assumes a spacing of 1, so if the actual spacing is lower the reported area's going to be higher
 %more on that, it makes the maximum for new_N_ (the first value)be equal to the initial coverage, N0_
rate_span = signal_span;

 %rate_span = signal_span * (N0/max_coverage_span_raw);
 %Works. 


%plot(time_span, signal_span, 'k--','LineWidth',5); %Raw rate
%plot(time_span, c_Pt_01, 'k--', 'LineWidth',5); %Raw rate coverage 
%plot(time_span, signal_span,'g--',time_span, coverage_span_raw, 'g--'); %Transformed raw rate and coverage
plot(time_span, rate_span,color) %Normalized transformed
plot(time_span, coverage_span_actual,color, 'LineWidth',2);


%Works.6/4 

end