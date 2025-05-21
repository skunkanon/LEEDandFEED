%% 4/11: Cleaning up first attempt. First, making a plot of coverage vs temperature. 
%Run 'coxlambertscan.m' first. 
x_0p8_raw = fig5_0p8(:,1)';
y_0p8_raw = fig5_0p8(:,2)';

scatter(x_0p8_raw, y_0p8_raw);
offset_0p8 = 5;

bg_pre_raw = y_0p8_raw(1:51-offset_0p8); %51 = 1000 K
bg_post_raw = y_0p8_raw(85:100); %85 = 1340 K 

%bg_post_raw = y_0p8_raw(85+offset_0p8:100); %85 = 1340 K 


bg_pre = mean(bg_pre_raw);
bg_post = mean(bg_post_raw);

% Linear background reduction 
templowBG_0p8 = x_0p8_raw(51-offset_0p8);
temphighBG_0p8 = x_0p8_raw(85);

%temphighBG_0p8 = x_0p8_raw(85+offset_0p8);

%test test 

bg = @(x) bg_pre + (x-templowBG_0p8)*((bg_post-bg_pre)/(temphighBG_0p8-templowBG_0p8)); %Temps that index 51 and 85 correspond to 

bg_span = linspace(templowBG_0p8,temphighBG_0p8,85-51+1+offset_0p8); %Number of indices from 51 to 85 
%bg_span = linspace(templowBG_0p8,temphighBG_0p8,85-51+1+offset_0p8*2); %Number of indices from 51 to 85 

x_0p8 = x_0p8_raw(51-offset_0p8:85); %restricting temperature range 
%x_0p8 = x_0p8_raw(51-offset_0p8:85+offset_0p8); %restricting temperature range 
figure(1);
hold on;
plot(x_0p8_raw, bg_pre * ones(size(x_0p8_raw)), '--k', 'LineWidth', 1.5);  % dashed black line
plot(x_0p8_raw, bg_post * ones(size(x_0p8_raw)), '--k', 'LineWidth', 1.5);
plot(bg_span,bg(bg_span),'b');
plot(x_0p8_raw,y_0p8_raw,'b');
hold off;

figure(1);
hold on;
y_0p8 = y_0p8_raw(51-offset_0p8:85) - bg(bg_span);
scatter(x_0p8,y_0p8,'b');
hold off;
%%
% Doing 0.4 




x_0p4_raw = fig5_0p4(:,1)';
y_0p4_raw = fig5_0p4(:,2)';
figure(1);
hold on;
scatter(x_0p4_raw, y_0p4_raw,'r');
hold off;
% 5/20, no offset, lowest coverage 
bg_pre_raw_04 = y_0p4_raw(1:49); %49 = 1000 K
bg_post_raw_04 = y_0p4_raw(83:98); %83 = 1340 K 


bg_pre_04 = mean(bg_pre_raw_04);
bg_post_04 = mean(bg_post_raw_04);

%Linear background reduction 

bg_04 = @(x) bg_pre_04 + (x-1000)*((bg_post_04-bg_pre_04)/(1340-1000)); %Temps that index 51 and 85 correspond to 

bg_span_04 = linspace(1000,1340,83-49+1); %Number of indices from 51 to 85 

x_0p4 = x_0p4_raw(49:83); %restricting temperature range 
figure(1);
hold on;
plot(x_0p4_raw, bg_pre_04 * ones(size(x_0p4_raw)), '--k', 'LineWidth', 1.5);  % dashed black line
plot(x_0p4_raw, bg_post_04 * ones(size(x_0p4_raw)), '--k', 'LineWidth', 1.5);
plot(bg_span_04,bg_04(bg_span_04),'r');
plot(x_0p4_raw,y_0p4_raw,'r');
hold off;

figure(1);
hold on;
y_0p4 = y_0p4_raw(49:83) - bg_04(bg_span_04);
scatter(x_0p4,y_0p4,'r');
hold off;
%%
% 4/11 1.2


x_1p2_raw = fig5_1p2(:,1)';
y_1p2_raw = fig5_1p2(:,2)';
figure(1);
hold on;
scatter(x_1p2_raw, y_1p2_raw,'m');
hold off;
offset_1p2 = 9;
bg_pre_raw_12 = y_1p2_raw(1:51-offset_1p2); %51 = 1000 K
bg_post_raw_12 = y_1p2_raw(84+offset_1p2:95); %84 = 1350 K 
templowBG_1p2 = x_1p2_raw(51-offset_1p2);
temphighBG_1p2 = x_1p2_raw(84+offset_1p2);

bg_pre_12 = mean(bg_pre_raw_12);
bg_post_12 = mean(bg_post_raw_12);

% Linear background reduction 

bg_12 = @(x) bg_pre_12 + (x-templowBG_1p2)*((bg_post_12-bg_pre_12)/(temphighBG_1p2-templowBG_1p2)); %Temps that index 51 and 85 correspond to 

bg_span_12 = linspace(templowBG_1p2,temphighBG_1p2,84-51+1+offset_1p2*2); %Number of indices from 51 to 85 

x_1p2 = x_1p2_raw(51-offset_1p2:84+offset_1p2); %restricting temperature range 
figure(1);
hold on;
plot(x_1p2_raw, bg_pre_12 * ones(size(x_1p2_raw)), '--k', 'LineWidth', 1.5);  % dashed black line
plot(x_1p2_raw, bg_post_12 * ones(size(x_1p2_raw)), '--k', 'LineWidth', 1.5);
plot(bg_span_12,bg_12(bg_span_12),'m');
plot(x_1p2_raw,y_1p2_raw,'m');
hold off;

figure(1);
hold on;
y_1p2 = y_1p2_raw(51:84) - bg_12(bg_span_12);
scatter(x_1p2,y_1p2,'m');
hold off;

title('Linear Background Reduction');
ylim([0 1]);

%% 4/15 - Determining coverage from dosage, with relation from LEED 
% Copied from 'appendthistoCL81...m'. Run 'coxlambertscan_fig2cov.m first'.
%Uses Figure 2's Auger ratio. Same as TDS up until sqrt(3) x sqrt(3)
%pattern. 


%scatter(normcoverage(:, 1), normcoverage(:, 2)); %regular plot function loops back in on itself 

% very crude method for calibrating the axes lol. bring the last number in x_interp at line 15 to 60,000 and find where it's 1.4 

% Step 1: Remove duplicate x-values before interpolation
[unique_x, idx] = unique(normcoverage(:,1)); 
unique_y = normcoverage(idx,2);

% Step 2: Sort the unique values
SNC = sortrows([unique_x, unique_y], 1);

% Step 3: Create finer interpolation grid 
x_interp = linspace(min(SNC(:,1)), max(SNC(:,1)), 2000)';

% Step 4: Interpolate (with extrapolation disabled for safety)
y_interp = interp1(SNC(:,1), SNC(:,2), x_interp, 'pchip');

% Step 5: Combine interpolated data
NC_interp = [x_interp, y_interp];


%Adjusting scale based off of LEED
scaleFactor = (1/3)/(0.589311184524498);



% Step 6: Plotting
figure(2);
hold on;
scatter(normcoverage(:,1), normcoverage(:,2) * scaleFactor, 'b'); % Original data
%plot(SNC(:,1), SNC(:,2), 'r-', 'LineWidth', 2); % Sorted unique data
plot(x_interp, y_interp * scaleFactor, 'r', 'LineWidth', 1.5); % Interpolated curve
hold off;

xlabel('Dosage (molecules/cm^2  * 10^{19})');
ylabel('Normalized Coverage');
title('Dosage vs Coverage');
legend('Original Data', 'Interpolated Original Data', 'Location', 'best');
grid on;

%% 4/15 Now have coverage from dosage. Zoom in to get equal coverage. 
% Kind of redundant with following section because of proportional
% relationship between desorption rate and mass spec signal.

N0_1p2 = 0.2877;
N0_0p4 = 0.0908;
N0_0p8 = 0.2004;

figure(3);
hold on;

N_0p8 = zeros(1,length(x_0p8));
for i = 1:length(x_0p8)
    N_0p8(i) = trapz(y_0p8(1:length(x_0p8) +1 - i));
end
[max_0p8, ~] = max(N_0p8);
N_0p8 = N_0p8 * (N0_0p8 / max_0p8);
plot(x_0p8, N_0p8, 'b');
scatter(x_0p8, N_0p8 , 'b');



N_0p4 = zeros(1,length(x_0p4));
for i = 1:length(x_0p4)
    N_0p4(i) = trapz(y_0p4(1:length(x_0p4) + 1 - i));
end
[max_0p4, ~] = max(N_0p4);
N_0p4 = N_0p4 * (N0_0p4 / max_0p4);
scatter(x_0p4, N_0p4, 'r');
plot(x_0p4 , N_0p4, 'r');


N_1p2 = zeros(1,length(x_1p2));
for i = 1:length(x_1p2)
    N_1p2(i) = trapz(y_1p2(1:length(x_1p2) + 1 - i));
end
[max_1p2, ~] = max(N_1p2);
N_1p2 = N_1p2 * (N0_1p2 / max_1p2);
scatter(x_1p2, N_1p2, 'm');
plot(x_1p2 , N_1p2 , 'm');

title('Coverage vs temperature');
ylim([0 1/3]);


%% 4/15  Taking background reduction from figure 1 and adding integrals of each isostere. 
fprintf('########### NEW INSTANCE \n ')

figure(4);
hold on;

S_0p4 = zeros(1,length(x_0p4));
for i = 1:length(x_0p4)
    S_0p4(i) = trapz(y_0p4(1:length(x_0p4) +1 - i))/length(S_0p4);
end
int_0p4 = trapz(x_0p4, y_0p4);
S_0p4 = S_0p4 * int_0p4/max(S_0p4);
scatter(x_0p4, S_0p4, 'r');
plot(x_0p4, S_0p4, 'r');
plot(x_0p4, y_0p4,'r');

S_0p8 = zeros(1,length(x_0p8));
for i = 1:length(x_0p8)
    S_0p8(i) = trapz(y_0p8(1:length(x_0p8) +1 - i))/length(S_0p8);
end
int_0p8 = trapz(x_0p8, y_0p8);
S_0p8 = S_0p8 * int_0p8/max(S_0p8);
scatter(x_0p8, S_0p8,'b');
plot(x_0p8, S_0p8,'b');
plot(x_0p8, y_0p8, 'b');


S_1p2 = zeros(1,length(x_1p2));
for i = 1:length(x_1p2)
    S_1p2(i) = trapz(y_1p2(1:length(x_1p2) + 1 - i))/length(S_1p2);
end
int_1p2 = trapz(x_1p2, y_1p2);
S_1p2 = S_1p2 * int_1p2/max(S_1p2);
scatter(x_1p2, S_1p2, 'm');
plot(x_1p2, S_1p2, 'm');
plot(x_1p2, y_1p2, 'm');
title('f and S of each isostere');
hold off;

%% 4/15 - Bit later in the day, doublechecking what's going on with trapz(). 
%The coverage vs temperature one's good since that already is scaled
%according to the LEED patterns vs dosage. 

zzzzz  = trapz(x_0p4 - min(x_0p4), y_0p4); %4/16 basically confirms that we're good to move on 
zzzzzz = trapz(x_0p4_interp - min(x_0p4_interp) , y_0p4_interp);



S_0p4 = zeros(1,length(x_0p4));



for i = 1:length(x_0p4)
    S_0p4(i) = trapz(y_0p4(1:length(x_0p4) +1 - i));
end
int_0p8 = trapz(x_0p4, y_0p4);
S_0p4 = S_0p4 * int_0p8/max(S_0p4);

for i = 1:length(x_0p4_interp)
    S_0p4_interp(i) = trapz(y_0p4_interp(1:length(x_0p4_interp) + 1 - i));
end

int_interp_0p4 = trapz(x_0p4_interp, y_0p4_interp);
S_0p4_interp = S_0p4_interp * int_interp_0p4/max(S_0p4_interp);


figure(6);
hold on; 
plot(x_0p4_interp, S_0p4_interp); 
plot(x_0p4, S_0p4);
hold off; %yep turn out to be basically the same after adjusting scaling like that huh whatever 


%% 4/15 - Interpolating values from previous plot to get more accurate 'f' and 'S'.  

figure(5);
hold on;
%0.4
interp = 10000;
x_0p4_interp = linspace(min(x_0p4), max(x_0p4), interp);
y_0p4_interp = interp1(x_0p4, y_0p4, x_0p4_interp, 'pchip');
%scatter(x_0p4_interp, y_0p4_interp);
S_0p4_interp = zeros(1,length(x_0p4_interp));


for i = 1:length(x_0p4_interp)
    S_0p4_interp(i) = trapz(y_0p4_interp(1:length(x_0p4_interp) + 1 - i));
end

int_interp_0p4 = trapz(x_0p4_interp, y_0p4_interp);
S_0p4_interp = S_0p4_interp * int_interp_0p4/max(S_0p4_interp);
%scatter(x_0p4_interp, S_0p4_interp, 'red');
plot(x_0p4_interp, S_0p4_interp, 'b');
plot(x_0p4_interp, y_0p4_interp,'b');


%0.8 
x_0p8_interp = linspace(min(x_0p8), max(x_0p8), interp);
y_0p8_interp = interp1(x_0p8, y_0p8, x_0p8_interp, 'pchip');
%scatter(x_0p4_interp, y_0p4_interp);
S_0p8_interp = zeros(1,length(x_0p8_interp));

for i = 1:length(x_0p8_interp)
    S_0p8_interp(i) = trapz(y_0p8_interp(1:length(x_0p8_interp) + 1 - i));
end


int_interp_0p8 = trapz(x_0p8_interp, y_0p8_interp);
S_0p8_interp = S_0p8_interp * int_interp_0p8/max(S_0p8_interp);

%scatter(x_0p4_interp, S_0p4_interp, 'red');
plot(x_0p8_interp, S_0p8_interp, 'r');
plot(x_0p8_interp, y_0p8_interp,'r');

x_1p2_interp = linspace(min(x_1p2), max(x_1p2), interp);
y_1p2_interp = interp1(x_1p2, y_1p2, x_1p2_interp, 'pchip');
%scatter(x_0p4_interp, y_0p4_interp);
S_1p2_interp = zeros(1,length(x_1p2_interp));
for i = 1:length(x_1p2_interp)
    S_1p2_interp(i) = trapz(y_1p2_interp(1:length(x_1p2_interp) + 1 - i));
end


int_interp_1p2 = trapz(x_1p2_interp, y_1p2_interp);
S_1p2_interp = S_1p2_interp * int_interp_1p2/max(S_1p2_interp);
%scatter(x_0p4_interp, S_0p4_interp, 'red');
plot(x_1p2_interp, S_1p2_interp, 'm');
plot(x_1p2_interp, y_1p2_interp,'m');




hold off;
%% 4/18 - Restricting temperature range to get linear relationship between isostere integrals. 
[~, I_0p4_choke] = max(y_0p4);
x_0p4_choke = x_0p4(1:I_0p4_choke);
y_0p4_choke = y_0p4(1:I_0p4_choke);
%scatter(x_0p4_choke, y_0p4_choke);
int_0p4_choke = trapz(x_0p4_choke, y_0p4_choke);

[~, I_0p8_choke] = max(y_0p8);
x_0p8_choke = x_0p8(1:I_0p8_choke);
y_0p8_choke = y_0p8(1:I_0p8_choke);
%scatter(x_0p8_choke, y_0p8_choke);
int_0p8_choke = trapz(x_0p8_choke, y_0p8_choke);

[~, I_1p2_choke] = max(y_1p2);
x_1p2_choke = x_1p2(1:I_1p2_choke);
y_1p2_choke = y_1p2(1:I_1p2_choke);
%scatter(x_1p2_choke, y_1p2_choke);
int_1p2_choke = trapz(x_1p2_choke, y_1p2_choke);

%Nope, still seems to double with each increment of 0.4 x 10^19. 

%% 5/19 - Isolating 'f' from isostere plots

figure(6);
hold on;
plot(x_0p4_interp, y_0p4_interp,'b');
plot(x_0p8_interp, y_0p8_interp,'r');
plot(x_1p2_interp, y_1p2_interp,'m');
hold off;


%% 4/16 Now doing analysis with interpolated plots. 
%0.078 coverage. 0.4 to 1.2 are at temperatures ~1065, 1162, and 1230 K. 
%Hold on - interpolated vs raw f and S plots line up, but a slight
%disagreement with coverage vs temperature plots. Integral over all signal
%should be proportional to the initial coverage, yeah? 
%4/18 - maybe only take stuff before peak? 




arrh_x = [1/1065, 1/1162, 1/1230];
arrh_y = [log(50*0.0382806/24.51), log(50*0.336007/26.39), log(50*0.457/33.6)];



% Linear regression with stats
X = [ones(length(arrh_x), 1), arrh_x];  % Design matrix
[b, b_int, ~, ~, stats] = regress(arrh_y, X);  % b = [intercept; slope]

% Extract results
m_0p78 = b(2);
b_0p78 = b(1);
R2 = stats(1);
slope_CI = b_int(2, :);
intercept_CI = b_int(1, :);

%predicted fit 
arrh_y_fit = X * b;

figure(7);
hold on;
set(gca, 'Ydir', 'reverse');
scatter(arrh_x, arrh_y);
plot(arrh_x, arrh_y_fit, 'DisplayName', sprintf('Regression: y = %.4fx + %.4f', m_0p78, b_0p78));
hold off;
