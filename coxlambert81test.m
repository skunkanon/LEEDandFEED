


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

plot(fig5_4(:, 1), fig5_4(:, 2));




plot(fig5_8(:, 1), fig5_8(:, 2));

plot(fig5_2p8(:, 1), fig5_2p8(:, 2));

plot(fig5_2p0(:, 1), fig5_2p0(:, 2));
plot(fig5_1p6(:, 1), fig5_1p6(:, 2));
plot(fig5_1p2(:, 1), fig5_1p2(:, 2));

plot(fig5_0p8(:, 1), fig5_0p8(:, 2));
plot(fig5_0p4(:, 1), fig5_0p4(:, 2));


hold off;

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