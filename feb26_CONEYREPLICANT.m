%{ 

%2/27 - ALL COMMENTED: JUST FOR SCANNING FROM DEANGELIS 92 

CoNi_SCAN = sortrows(CoNi_SCAN,1);
ratespan_CONI = CoNi_SCAN(:,2);
timespan_CONI = CoNi_SCAN(:,1);
scatter(timespan_CONI, ratespan_CONI);

%%

% x,y are column vectors (same length), ordered along the curve
t = [0; cumsum(hypot(diff(timespan_CONI), diff(ratespan_CONI)))];   % arc-length parameter

p = 0.95;  % smoothing parameter: 1 = interpolate, smaller = smoother
spx = csaps(t, timespan_CONI, p);   % requires Spline Toolbox
spy = csaps(t, ratespan_CONI, p);

tq = linspace(t(1), t(end), 500);
timespan_CONI = fnval(spx, tq);
ratespan_CONI = fnval(spy, tq);

%plot(timespan_CONI, ratespan_CONI, '.', timespan_CONI, ratespan_CONI, '-', 'LineWidth', 1.5);

%% final bit of scanning thing 
ratespan_CONI(ratespan_CONI < 0) = 0;
figure(1); clf;
hold on;
scatter(timespan_CONI,ratespan_CONI);
plot(timespan_CONI, ratespan_CONI);
hold off;

%%

ratespan_CONI = ratespan_CONI';
timespan_CONI = timespan_CONI';


%%

figure(1); clf;
hold on;
scatter(timespan_CONI,ratespan_CONI);
plot(timespan_CONI, ratespan_CONI);
hold off;

%%

ratespan_CONI = ratespan_CONI(:);
timespan_CONI = timespan_CONI(:);

data = [timespan_CONI ratespan_CONI];
writematrix(data, 'CONEYISLAND.csv', 'Delimiter',',');
%%
raw = clipboard('paste');              % grabs your tab/space separated text

v = sscanf(raw, '%f%f');               % reads pairs of numbers (whitespace-separated)
data = reshape(v, 2, []).';            % Nx2 matrix: [x y]

% Make MATLAB matrix literal text
lines = compose('  %.15g, %.15g;', data(:,1), data(:,2));
txt = "data = [" + newline + strjoin(lines, newline) + newline + "];";

clipboard('copy', txt);                % now paste into your script/editor

%}

CoNi_SCAN = sortrows(CoNi_SCAN,1);
ratespan_CONI = CoNi_SCAN(:,2);
timespan_CONI = CoNi_SCAN(:,1);
scatter(timespan_CONI, ratespan_CONI);

ratespan_CONI = ratespan_CONI / trapz(timespan_CONI, ratespan_CONI); 
% integrates to 1 for initial coverage at 1 
%% coverage vs time 3/2/26


% cumulative desorbed coverage:  ∫_t0^t r(t') dt'
theta_desorbed  = cumtrapz(timespan_CONI, ratespan_CONI);

% remaining coverage: theta(t) = 1 - ∫ r dt
theta_remaining = 1 - theta_desorbed;

figure(1); clf;
plot(timespan_CONI, theta_remaining, 'LineWidth', 1.5);
xlabel('time');
ylabel('remaining coverage \theta(t)');
ylim([0 1]);
grid on;

%%  Rate vs remaining coverage (rate is dependent variable), fig 4 

figure(2); clf;
scatter(theta_remaining, ratespan_CONI, 'LineWidth', 1.5);
xlabel('remaining coverage \theta');
ylabel('rate r(\theta)');
xlim([0 1]);
grid on;

% TO DO - polyani wigner fitting algorithm for figure 3 (get prefactor + Ea
% right there); 




%%
hodl = struct();
hodl.t1 = [35 95]; hodl.T1 = 200; 
hodl.t2 = [112 190]; hodl.T2 = 300;



exp_data = cell(1,1);
exp_data{1} = {timespan_CONI, ratespan_CONI, 1, ones(size(ratespan_CONI)), hodl};

init1 = [30e3, -20e3, 1e15, 1000];
init2 = [130e3, -20e3, 1e15, 1000];

des_order = 1; 

[best1, best2, fit_error, hess1, hess2, cov1, cov2, se1, se2] = ...
    fit_pw_niemant_isothermal_twostep(exp_data, init1, init2, des_order);

