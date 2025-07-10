function [temperature_span,signal_span, BG_std] = reduceBGgetspectra(x_raw, y_raw,lowTtarget, highTtarget, color)
%Variable definitions say temperature, but can take time bounds as well as temperature bounds. 

templowBGtarget = lowTtarget; 
temphighBGtarget = highTtarget; 
new_x_raw = x_raw; 
new_y_raw = y_raw; 
%Redundancy in variable definition just to show where code's coming from, CLAA_BG_reduce. 


[~, templowBG_ind] = min(abs(new_x_raw-templowBGtarget)); 
new_templowBG = new_x_raw(templowBG_ind);
new_BG_pre_raw = new_y_raw(1:templowBG_ind);
new_BG_pre = mean(new_BG_pre_raw);

[~, temphighBG] = min(abs(new_x_raw-temphighBGtarget));
new_temphighBG = new_x_raw(temphighBG);
new_BG_post_raw = new_y_raw(temphighBG: length(new_y_raw));
new_BG_post = mean(new_BG_post_raw);

new_BG = @(x) new_BG_pre + (x - new_templowBG) * ...
    ((new_BG_post - new_BG_pre) ...
    /(new_temphighBG - new_templowBG));

new_BG_span = linspace(new_templowBG , new_temphighBG ...
    ,(temphighBG - templowBG_ind + 1));




temperature_span = new_x_raw(templowBG_ind:temphighBG);
signal_span = new_y_raw(templowBG_ind : temphighBG) - new_BG(new_BG_span);
fprintf('%.4f\n' , trapz(signal_span, temperature_span));
plot(x_raw, new_BG_pre * ones(size(x_raw)), 'k--');
plot(x_raw, new_BG_post * ones(size(x_raw)), 'k--');
plot(new_BG_span, new_BG(new_BG_span), color);
plot(x_raw, y_raw, color);
plot(temperature_span, signal_span, color, 'LineWidth',2);

%FOR WEIGHTED FITTING, 7/9

BG_std = std(new_BG_pre_raw) * sqrt(1 + abs(signal_span) / new_BG_pre);


end