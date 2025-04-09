


%% 3/18 - Cox Lambert '81 Test 

[t,x,y] = polyani_wigner(50,50,247*10^3, 0,10^13,1.3,1500); %third variable is Ea in j/mol

plot(x,y);

%3/18, It works

%% 3/18 - Fig 5, Cox Lambert '81
fprintf('NEW INSTANCE \n');

Rh111_rho = 1.6 * 10 ^19; %Looked up 111 FCC surface area and Rh lattice constant 

Cl_dose = [0.4,0.8,1.2,2.0,2.8,4.0,8.0] .* 10^19; 

N0_span = Cl_dose ./ Rh111_rho;
%%  3/18 - Makes all coverages above saturation become 1. 


for i = 1:length(N0_span)
    if N0_span(i) > 1
        N0_span(i) = 1;
    end
end

%% 3/18 - Make Plot 
%Authors report Rh desorption energy drops from 249 to 237 kj/mol from ~0
%to 1/3 coverage. That's a 12 kJ difference from 0 to 1/3, so 36 kJ from 0
%to 1. Using 36 kJ or 36 * 10^3 J (unit consistency) for gamma (y_E in the
%function). 


[ti_span, x_span,y_span] = deal(cell(1,length(N0_span)));

figure;
hold on;
for i=1:length(N0_span)
    [ti_span{i},x_span{i},y_span{i}] = polyani_wigner(50, 0, 247*10^3, 36 * 10^3, 10^13, N0_span(i), 1500);
    plot(x_span{i},y_span{i},'DisplayName', sprintf('N = %.3f', N0_span(i)));
end

hold off;

%3/18 - Coverage - energy relation isn't checking out, and kind of clueless as to what happens
%past saturation. Moving on to another publication for now. 
%Also, don't know how temperature at exposure is going to be taken into
%account in the math. 


%% 3/28 - Digitizing scanned plots. Works. Below is Fig 5 for 300 K. 

figure;
hold on;

%plot(fig5_4(:, 1), fig5_4(:, 2));




%plot(fig5_8(:, 1), fig5_8(:, 2));

%plot(fig5_2p8(:, 1), fig5_2p8(:, 2));

%plot(fig5_2p0(:, 1), fig5_2p0(:, 2));
%plot(fig5_1p6(:, 1), fig5_1p6(:, 2));
plot(fig5_1p2(:, 1), fig5_1p2(:, 2), 'DisplayName', sprintf('1.2 x 10^{19} Dose / 0.28 Coverage'), 'LineWidth', 2);

plot(fig5_0p8(:, 1), fig5_0p8(:, 2), 'DisplayName',sprintf('0.8 x 10^{19} Dose / 0.19 Coverage'), 'LineWidth',2);
plot(fig5_0p4(:, 1), fig5_0p4(:, 2), 'DisplayName', sprintf('0.4 x 10^{19} Dose / 0.095 Coverage '), 'LineWidth',2);
xlabel('Temperature (K)')
ylabel('Normalized 35 amu Mass Spec Signal Intensity')
title('Digitized Fig 5a, Coverage < 1/3 ')
hold off;
legend('show');


%% - 4/8, FUCK the cell array lol 

%doing 0.8 first 
x_0p8 = fig5_0p8(:,1);
y_0p8 = fig5_0p8(:,2);
[min_0p8, ~] = min(y_0p8);
[M_0p8, ~] = max(y_0p8);
cutoff_0p8 = (M_0p8-min_0p8)/10;
%{

for i=1:length(y_0p8)
    y_0p8(i) = y_0p8(i) - min_0p8;
    if y_0p8(i) < cutoff_0p8
        y_0p8(i) = 0;
    end
end

%}
%doing 1.2 
x_1p2 = fig5_1p2(:,1);
y_1p2 = fig5_1p2(:,2);
[min_1p2, ~] = min(y_1p2);
[M_1p2, ~] = max(y_1p2);
cutoff_1p2 = (M_1p2-min_1p2)/10;
%{

for i=1:length(y_1p2)
    y_1p2(i) = y_1p2(i) - min_1p2;
    if y_1p2(i) < cutoff_1p2
        y_1p2(i) = 0;
   
    end
end
%}

%doing 0.4
x_0p4 = fig5_0p4(:,1);
y_0p4 = fig5_0p4(:,2);
[min_0p4, ~] = min(y_0p4);
[M_0p4, ~] = max(y_0p4);
cutoff_0p4 = (M_0p4-min_0p4)/10;

%{
for i=1:length(y_0p4)
    y_0p4(i) = y_0p4(i) - min_0p4;
    if y_0p4(i) < cutoff_0p4
        y_0p4(i) = 0;
   
    end
end
%}
figure;
hold on;
plot(x_1p2, y_1p2, 'DisplayName', sprintf('1.2 x 10^{19} Dose / 0.28 Coverage'), 'LineWidth', 2);

plot(x_0p8, y_0p8, 'DisplayName',sprintf('0.8 x 10^{19} Dose / 0.19 Coverage'), 'LineWidth',2);
plot(x_0p4, y_0p4, 'DisplayName', sprintf('0.4 x 10^{19} Dose / 0.095 Coverage '), 'LineWidth',2);
xlabel('Temperature (K)')
ylabel('Normalized 35 amu Mass Spec Signal Intensity')
title('Digitized Fig 5a, Coverage < 1/3 ')
hold off;
legend('show');

%% 4/8, finishing background subtraction /yeah no this sucks
%{
x_0p4 = x_0p4(208:338);
y_0p4 = y_0p4(208:338);

x_1p2 = x_1p2(163:273);
y_1p2 = y_1p2(163:273);

x_0p8 = x_0p8(203:339);
y_0p8 = y_0p8(203:339);


[min_0p8,~] = min(y_0p8);

for i=1:length(y_0p8)
    y_0p8(i) = y_0p8(i) - min_0p8;
end

x_1p2 = x_1p2(1:105);
y_1p2 = y_1p2(1:105);


[min_1p2,~] = min(y_1p2);



for i=1:length(y_1p2)
    y_1p2(i) = y_1p2(i) - min_1p2;
end


figure;
hold on;
scatter(x_1p2,y_1p2, 'DisplayName', '');
scatter(x_0p8,y_0p8);
scatter(x_0p4,y_0p4);
hold off;

%}

figure;
hold on;
scatter(x_1p2, y_1p2, 'DisplayName', sprintf('1.2 x 10^{19} Dose / 0.28 Coverage'));

scatter(x_0p8, y_0p8, 'DisplayName',sprintf('0.8 x 10^{19} Dose / 0.19 Coverage'));
scatter(x_0p4, y_0p4, 'DisplayName', sprintf('0.4 x 10^{19} Dose / 0.095 Coverage '));
xlabel('Temperature (K)')
ylabel('Normalized 35 amu Mass Spec Signal Intensity')
title('Background-Subtracted Fig 5a, Coverage < 1/3 ')
hold off;
legend('show');
%% 4-8 THE MATH PART WOOOO

f_0p8 = y_0p8';
S_0p8 = zeros(1,length(f_0p8));

f_0p4 = y_0p4';
S_0p4 = zeros(1,length(f_0p4));

f_1p2 = y_1p2';
S_1p2 = zeros(1,length(f_1p2));

for i=1:length(f_0p8)
    S_0p8(i) = trapz(f_0p8(1:(length(f_0p8) - i)));
end


for i=1:length(f_0p4)
    S_0p4(i) = trapz(f_0p4(1:(length(f_0p4) - i)));
end

for i=1:length(f_1p2)
    S_1p2(i) = trapz(f_1p2(1:(length(f_1p2) - i)));
end


%%
fprintf('NEW INSTANCE \n');
x1 = 1 ./(x_0p4');
y1 = log( 50 .* f_0p4 ./ S_0p4);

x2 = 1 ./(x_0p8');
y2 = log( 50 .* f_0p8 ./ S_0p8);

x3 = 1 ./(x_1p2');
y3 = log( 50 .* f_1p2 ./ S_1p2);
%%
figure;
hold on;
scatter(x1,y1, 'DisplayName',sprintf('1.2 x 10^{19} Dose / 0.28 Coverage'));
scatter(x2,y2, 'DisplayName',sprintf('0.8 x 10^{19} Dose,/ 0.19 Coverage'));
scatter(x3,y3, 'DisplayName', sprintf('0.4 x 10^{19} Dose / 0.095 Coverage'));
set(gca, 'YDir','reverse');
xlabel('1/Temperature (K)');
ylabel('ln(Beta*f/S)');
legend('show')
title('Raw Analysis')
hold off;


%% ok actual BG subtraction now 
fprintf('NEW INSTANCE \n');

%doing 0.8 first 
x_0p8_BG = fig5_0p8(:,1);
y_0p8_BG = fig5_0p8(:,2);
[min_0p8_BG, ~] = min(y_0p8_BG);
[M_0p8_BG, ~] = max(y_0p8_BG);
cutoff_0p8 = (M_0p8_BG-min_0p8_BG)/10;


for i=1:length(y_0p8_BG)
    y_0p8_BG(i) = y_0p8_BG(i) - min_0p8;
    if y_0p8_BG(i) < cutoff_0p8
        y_0p8_BG(i) = 0;
    end
end


%0.12
x_1p2_BG = fig5_1p2(:,1);
y_1p2_BG = fig5_1p2(:,2);
[min_1p2_BG, ~] = min(y_1p2_BG);
[M_1p2_BG, ~] = max(y_1p2_BG);
cutoff_1p2 = (M_1p2_BG-min_1p2_BG)/10;


for i=1:length(y_1p2_BG)
    y_1p2_BG(i) = y_1p2_BG(i) - min_1p2;
    if y_1p2_BG(i) < cutoff_1p2
        y_1p2_BG(i) = 0;
    end
end

%}
%0.4 
x_0p4_BG = fig5_0p4(:,1);
y_0p4_BG = fig5_0p4(:,2);
[min_0p4_BG, ~] = min(y_0p4_BG);
[M_0p4_BG, ~] = max(y_0p4_BG);
cutoff_0p4 = (M_0p4_BG-min_0p4_BG)/10;


for i=1:length(y_0p4_BG)
    y_0p4_BG(i) = y_0p4_BG(i) - min_0p4;
    if y_0p4_BG(i) < cutoff_0p4
        y_0p4_BG(i) = 0;
    end
end


%}
%%


%doing 0.8 first 
x_0p8_BG = fig5_0p8(:,1);
y_0p8_BG = fig5_0p8(:,2);
[min_0p8_BG, ~] = min(y_0p8_BG);
[M_0p8_BG, ~] = max(y_0p8_BG);
cutoff_0p8 = (M_0p8_BG-min_0p8_BG)/10;


for i=1:length(y_0p8_BG)
    y_0p8_BG(i) = y_0p8_BG(i) - min_0p8;
    if y_0p8_BG(i) < cutoff_0p8
        y_0p8_BG(i) = 0;
    end
end


%%
figure;
hold on;
scatter(x_1p2_BG, y_1p2_BG,'DisplayName',sprintf('1.2 x 10^{19} Dose / 0.28 Coverage'));
scatter(x_0p8_BG, y_0p8_BG,'DisplayName',sprintf('0.8 x 10^{19} Dose / 0.19 Coverage'));

scatter(x_0p4_BG, y_0p4_BG,'DisplayName',sprintf('0.4 x 10^{19} Dose / 0.095 Coverage'));

hold off;
%% WHY DOES IT LOOK LIKE THAT FUCK ME  


fprintf('NEW INSTANCE \n');
x_0p8_BG = x_0p8_BG(51:81);
y_0p8_BG = y_0p8_BG(51:81);

%% MADE THE PLOT WITH 1.2 LOWER THAN THE REST OF THEM SOMEHOW WHAT 

fprintf('NEW INSTANCE \n');
f_0p8_BG = y_0p8_BG';
S_0p8_BG = zeros(1,length(f_0p8_BG));


for i=1:length(f_0p8_BG)
    S_0p8_BG(i) = trapz(f_0p8_BG(1:(length(f_0p8_BG) - i)));
end


x2_BG = 1 ./(x_0p8_BG');
y2_BG = log( 50 .* f_0p8_BG ./ S_0p8_BG);
%{
% 0.4
f_0p4_BG = y_0p4_BG';
S_0p4_BG = zeros(1,length(f_0p4_BG));


for i=1:length(f_0p4_BG)
    S_0p4_BG(i) = trapz(f_0p4_BG(1:(length(f_0p4_BG) - i)));
end


x1_BG = 1 ./(x_0p4_BG');
y1_BG = log( 50 .* f_0p4_BG ./ S_0p4_BG);

% 1.2
f_1p2_BG = y_1p2_BG';
S_1p2_BG = zeros(1,length(f_1p2_BG));


for i=1:length(f_1p2_BG)
    S_1p2_BG(i) = trapz(f_1p2_BG(1:(length(f_1p2_BG) - i)));
end


x3_BG = 1 ./(x_1p2_BG');
y3_BG = log( 50 .* f_1p2_BG ./ S_1p2_BG);
 
%}
x2_BG = x2_BG(1:29);
y2_BG = y2_BG(1:29);

p_0p8 = polyfit(x2_BG, y2_BG,1);
m_0p8 = p_0p8(1);
b_0p8 = p_0p8(2);
y2_BG_fit = polyval(p_0p8,x2_BG);


figure;
hold on;
%scatter(x3_BG, y3_BG,'DisplayName',sprintf('1.2 x 10^{19} Dose / 0.28 Coverage'));
scatter(x2_BG, y2_BG,'DisplayName',sprintf('0.8 x 10^{19} Dose / 0.19 Coverage')); %huh?????
%scatter(x1_BG, y1_BG,'DisplayName',sprintf('0.4 x 10^{19} Dose / 0.095 Coverage'));
plot(x2_BG, y2_BG_fit, 'DisplayName', sprintf('Regression: y = %.4fx + %.4f',m_0p8,b_0p8));

set(gca, 'YDir','reverse');
title('Analysis')
legend('show')
hold off;
%%

x2_BG = x2_BG(1:29);
y2_BG = y2_BG(1:29);

% Linear regression with stats
X = [ones(length(x2_BG), 1), x2_BG];  % Design matrix
[b, b_int, ~, ~, stats] = regress(y2_BG, X);  % b = [intercept; slope]

% Extract results
m_0p8 = b(2);
b_0p8 = b(1);
R2 = stats(1);
slope_CI = b_int(2, :);
intercept_CI = b_int(1, :);

% Predicted fit
y2_BG_fit = X * b;

% Plot
figure;
hold on;
scatter(x2_BG, y2_BG, 'DisplayName', sprintf('0.8 x 10^{19} Dose / 0.19 Coverage'));
plot(x2_BG, y2_BG_fit, 'DisplayName', ...
    sprintf('Regression: y = %.4fx + %.4f', m_0p8, b_0p8));
set(gca, 'YDir','reverse');
title('Analysis');
legend('show');
hold off;

% Print stats
fprintf('Slope: %.4f\n', m_0p8);
fprintf('95%% CI for slope: [%.4f, %.4f]\n', slope_CI(1), slope_CI(2));
fprintf('Intercept: %.4f\n', b_0p8);
fprintf('95%% CI for intercept: [%.4f, %.4f]\n', intercept_CI(1), intercept_CI(2));
fprintf('R² = %.4f\n', R2);



%% Making cell array from digitized points 


vars = who('fig5_8', 'fig5_4', 'fig5_2p8', 'fig5_2p4' , 'fig5_2p0', 'fig5_1p6', 'fig5_1p2', 'fig5_0p8', 'fig5_0p4');

% Initialize cell array
data_cell = cell(1, numel(vars));

% Populate cell array
for i = 1:numel(vars)
    data_cell{i} = eval(vars{i});
end

%%
testtest = data_cell{1};
testtesttest = testtest(:, 1)';




%% Plotting cell array 

fprintf('NEW INSTANCE \n');

figure;
hold on;  % Enable overlaying plots
colors = parula(numel(data_cell));  % Assign distinct colors

for i = 1:numel(data_cell)
    x = data_cell{i}(:, 1);  % First column (x-values)
    y = data_cell{i}(:, 2);  % Second column (y-values)
    plot(x, y, 'LineWidth', 1.5, 'Color', colors(i, :), ...
         'DisplayName', sprintf('Dataset %d', i));
end

hold off;
grid on;
xlabel('X-axis');
ylabel('Y-axis');
title('Multiple Datasets');
legend('show');  % Show legend with dataset labels


%% 4/2 - Deriving activation energy, pre-exponential 
fprintf('\n NEW INSTANCE \n');
%Getting peak temp

[Tp, Tp_index] = deal(zeros(1,numel(data_cell)));


%Getting peak temps 
for i = 1:numel(data_cell) 
    x = data_cell{i}(:, 1);  % First column (x-values)
    y = data_cell{i}(:, 2);  % Second column (y-values)
    [~, Tp_index(i)] = max(y);
    Tp(i) = x(Tp_index(i));
end

%Calculating plot values

[log1_x, log1_y] = deal(zeros(1,numel(data_cell)));
beta = 50; %K/s
for i = 1:numel(data_cell)
    log1_x(i) = 1/Tp(i);
    log1_y(i) = log(beta/Tp(i)^2);
end

%Plotting 


p = polyfit(log1_x, log1_y,1);

m = p(1);
b = p(2);

fit_log_y = polyval(p,log1_x);

figure;
hold on;
plot(log1_x, log1_y, 'bo', 'DisplayName', 'Raw Data');
plot(log1_x, fit_log_y, 'r', 'DisplayName', sprintf('Fit Line: y = %.2fx + %.2f', m , b));
xlabel('1/Tp');
ylabel('ln(Beta/Tp^2');
legend('show');
hold off;


%Calculate Ea

R = 8.314; % J/molK
%kcal_to_j = 4184;
CoxL_Ea = (-m*R); %okay this sucks it's way too low and NEGATIVE what 


 
%% 4-3 - Deriving coverage from LEED 

%