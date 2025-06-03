%Time trial - with insights from 'erley80test', trying to do Arrhenius
%analysis with actual time integrals 
%Run CLAABGreduce first. 



%Plotting the background-reduced signals 
n_trace = 8;
cmap = parula(n_trace);

temp_init = 300; %K 
beta = 50; %K/s
time = @(x) (x-temp_init)/beta;

figure(1); clf;
hold on;
plot(time(new_x_0p4), new_y_0p4, 'o-');
plot(time(new_x_0p8), new_y_0p8, 'o-');
plot(time(new_x_1p2), new_y_1p2, 'o-');
plot(time(new_x_1p6), new_y_1p6, 'o-');
plot(time(new_x_2p0), new_y_2p0, 'o-');
plot(time(new_x_2p8), new_y_2p8, 'o-');
plot(time(new_x_4p0), new_y_4p0, 'o-');
plot(time(new_x_8p0), new_y_8p0, 'o-');

set(gca, 'ColorOrder', cmap, 'NextPlot', 'replacechildren');
legend('0.4', '0.8', '1.2', '1.6', '2.0', '2.8', '4.0', '8.0');
xlabel('Time (s)');
ylabel('Spectrometer Signal');
title('Desorption Traces, Signal vs Time')
hold off;

%% Getting outputs of time() of lower coverage desorption spectra

t_0p4 = time(new_x_0p4);
t_0p8 = time(new_x_0p8);
t_1p2 = time(new_x_1p2);
t_1p6 = time(new_x_1p6);

%%
fprintf('NEW INSTANCE \n');

