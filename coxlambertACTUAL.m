%% 4/11: Cleaning up first attempt. First, making a plot of coverage vs temperature. 
%Run 'coxlambertscan.m' first. 
x_0p8_raw = fig5_0p8(:,1)';
y_0p8_raw = fig5_0p8(:,2)';

%scatter(x_0p8_raw, y_0p8_raw);
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
figure(13);
hold on;
plot(x_0p8_raw, bg_pre * ones(size(x_0p8_raw)), '--k', 'LineWidth', 1.5);  % dashed black line
plot(x_0p8_raw, bg_post * ones(size(x_0p8_raw)), '--k', 'LineWidth', 1.5);
plot(bg_span,bg(bg_span),'b');
plot(x_0p8_raw,y_0p8_raw,'b');
hold off;

figure(13);
hold on;
y_0p8 = y_0p8_raw(51-offset_0p8:85) - bg(bg_span);
scatter(x_0p8,y_0p8,'b');
hold off;

% Doing 0.4 




x_0p4_raw = fig5_0p4(:,1)';

y_0p4_raw = fig5_0p4(:,2)';
figure(13);
hold on;
%scatter(x_0p4_raw, y_0p4_raw,'r');
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
figure(13);
hold on;
plot(x_0p4_raw, bg_pre_04 * ones(size(x_0p4_raw)), '--k', 'LineWidth', 1.5);  % dashed black line
plot(x_0p4_raw, bg_post_04 * ones(size(x_0p4_raw)), '--k', 'LineWidth', 1.5);
plot(bg_span_04,bg_04(bg_span_04),'r');
plot(x_0p4_raw,y_0p4_raw,'r');
hold off;

figure(13);
hold on;
y_0p4 = y_0p4_raw(49:83) - bg_04(bg_span_04);
scatter(x_0p4,y_0p4,'r');
hold off;

% 4/11 1.2


x_1p2_raw = fig5_1p2(:,1)';
y_1p2_raw = fig5_1p2(:,2)';
figure(13);
hold on;
%scatter(x_1p2_raw, y_1p2_raw,'m');
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
figure(13);
hold on;
plot(x_1p2_raw, bg_pre_12 * ones(size(x_1p2_raw)), '--k', 'LineWidth', 1.5);  % dashed black line
plot(x_1p2_raw, bg_post_12 * ones(size(x_1p2_raw)), '--k', 'LineWidth', 1.5);
plot(bg_span_12,bg_12(bg_span_12),'m');
plot(x_1p2_raw,y_1p2_raw,'m');
hold off;

figure(13);
hold on;
y_1p2 = y_1p2_raw(51-offset_1p2:84+offset_1p2) - bg_12(bg_span_12);
scatter(x_1p2,y_1p2,'m');
hold off;

title('Linear Background Reduction');
ylim([0 2]);
% 5/21 - Messed up the axes scale on the original scans, but only fixed for the lower coverage ones. Fixing that for 1.6 x 10^19 here. 
%Should just redo all the scans? 

fig5_1p6_x_corrected = fig5_1p6(:,1) ./(15/12) + 300;
fig5_1p6_y_corrected = (fig5_1p6(:,2) * 2.5);
%clf(10);
%figure(10);
%hold on;
%scatter(fig5_0p4(:,1),fig5_0p4(:,2));
%scatter(fig5_0p8(:,1),fig5_0p8(:,2));
%scatter(fig5_1p2(:,1),fig5_1p2(:,2));
%scatter(fig5_1p6_x_corrected,fig5_1p6_y_corrected);
%scatter(fig5_2p0(:,1),fig5_2p0(:,2));
%scatter(fig5_2p8(:,1),fig5_2p8(:,2));
%hold off;

 


% 5/20, 1.6 



x_1p6_raw = fig5_1p6_x_corrected;
y_1p6_raw = fig5_1p6_y_corrected;
%clf(1);
figure(13);
hold on;
%scatter(x_1p6_raw, y_1p6_raw,'g');
hold off;
offset_1p6 = 12 * 3;
bg_pre_raw_16 = y_1p6_raw(1:51 * 3-offset_1p6); %51 = 1000 K, multiplied by three because of higher sampling rate
bg_post_raw_16 = y_1p6_raw(84 * 3: length(x_1p6_raw)); %84 = 1350 K 
templowBG_1p6 = x_1p6_raw(51 * 3 -offset_1p6);
temphighBG_1p6 = x_1p6_raw(84 *3);

bg_pre_16 = mean(bg_pre_raw_16);
bg_post_16 = mean(bg_post_raw_16);

% Linear background reduction 

bg_16 = @(x) bg_pre_16 + (x-templowBG_1p6)*((bg_post_16-bg_pre_16)/(temphighBG_1p6-templowBG_1p6)); %Temps that index 51 and 85 correspond to 

bg_span_16 = linspace(templowBG_1p6,temphighBG_1p6,84*3-51*3+1+offset_1p6); %Number of indices from 51 to 85 
bg_span_16 = bg_span_16';

x_1p6 = x_1p6_raw(51*3-offset_1p6:84*3); %restricting temperature range 
figure(13);
hold on;
plot(x_1p6_raw, bg_pre_16 * ones(size(x_1p6_raw)), '--k', 'LineWidth', 1.5);  % dashed black line
plot(x_1p6_raw, bg_post_16 * ones(size(x_1p6_raw)), '--k', 'LineWidth', 1.5);
plot(bg_span_16,bg_16(bg_span_16),'g');
plot(x_1p6_raw,y_1p6_raw,'g');
hold off;
%
fprintf('NEW INSTANCE \n');
figure(13);
hold on;
y_1p6 = y_1p6_raw(51*3-offset_1p6:84*3) - bg_16(bg_span_16);
scatter(x_1p6,y_1p6,'g');
hold off;

title('Linear Background Reduction');
ylim([0 1.3]);
xlim([800 1500]);

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
N0_1p6 = 0.3856; %5/21 
figure(4);
clf(4);
hold on;

N_0p8 = zeros(1,length(x_0p8));
for i = 1:length(x_0p8)
    N_0p8(i) = trapz(y_0p8(1:length(x_0p8) +1 - i));
end
[max_0p8, ~] = max(N_0p8);
N_0p8 = N_0p8 * (N0_0p8 / max_0p8);
plot(x_0p8, N_0p8, 'b');
scatter(x_0p8, N_0p8 , 'b');
plot(x_0p8, y_0p8 .* N0_0p8,'b');


N_0p4 = zeros(1,length(x_0p4));
for i = 1:length(x_0p4)
    N_0p4(i) = trapz(y_0p4(1:length(x_0p4) + 1 - i));
end
[max_0p4, ~] = max(N_0p4);
N_0p4 = N_0p4 * (N0_0p4 / max_0p4);
scatter(x_0p4, N_0p4, 'r');
plot(x_0p4 , N_0p4, 'r');
plot(x_0p4, y_0p4 .* N0_0p4,'r');

N_1p2 = zeros(1,length(x_1p2));
for i = 1:length(x_1p2)
    N_1p2(i) = trapz(y_1p2(1:length(x_1p2) + 1 - i));
end
[max_1p2, ~] = max(N_1p2);
N_1p2 = N_1p2 * (N0_1p2 / max_1p2);
scatter(x_1p2, N_1p2, 'm');
plot(x_1p2 , N_1p2 , 'm');
plot(x_1p2, y_1p2 .* N0_1p2,'m');

N_1p6 = zeros(1,length(x_1p2)); %yeah this is screwed 5/21 re-scan the plots. inflection point (peak temp) is screwed)  
for i = 1:length(x_1p6)
    N_1p6(i) = trapz(y_1p6(1:length(x_1p6) + 1 - i));
end
[max_1p6, ~] = max(N_1p6);
N_1p6 = N_1p6 * (N0_1p6 / max_1p6);
%scatter(x_1p6, N_1p6, 'g'); %commenting out these because screwed
%plot(x_1p6 , N_1p6 , 'g');


title('Coverage vs temperature');
ylim([0 0.5]);


%% 4/15  Taking background reduction from figure 1 and adding integrals of each isostere. 

n_trace = 8;
cmap = parula(n_trace);

figure(1); clf;
hold on;
plot(new_x_0p4, new_y_0p4, 'o-');
plot(new_x_0p8, new_y_0p8, 'o-');
plot(newnew_x_1p2, newnew_y_1p2, 'o-');
plot(newnew_x_1p6, newnew_y_1p6, 'o-');
plot(newnew_x_2p0, newnew_y_2p0, 'o-');
plot(new_x_2p8, new_y_2p8, 'o-');
plot(new_x_4p0, new_y_4p0, 'o-');
plot(new_x_8p0, new_y_8p0, 'o-');
set(gca, 'ColorOrder', cmap, 'NextPlot', 'replacechildren');
legend('0.4', '0.8', '1.2', '1.6', '2.0', '2.8', '4.0', '8.0');
xlabel('Temperature (K)');
ylabel('Spectrometer Signal');
title('Desorption Traces, Signal vs Temperature')
hold off;
%%
% 5/21 - Multiplying the signal by the initial coverage to force the signal
% and desorption rate proportional. Use these values for F and S in the
% Arrhenius plot. 
%5/21, later in the day: this sucks nevermind
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
    S_0p8(i) = trapz(N0_0p8 .* y_0p8(1:length(x_0p8) +1 - i))/length(S_0p8);
end
int_0p8 = trapz(x_0p8, y_0p8 .* N0_0p8);
S_0p8 = S_0p8 * int_0p8/max(S_0p8);
scatter(x_0p8, S_0p8,'b');
plot(x_0p8, S_0p8,'b');
plot(x_0p8, y_0p8 .* N0_0p8, 'b');


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

%% 5/21 - Properly scaled S and F plots

N0_1p2 = 0.2877;
N0_0p4 = 0.0908;
N0_0p8 = 0.2004;
N0_1p6 = 0.3856; %5/21 
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
plot(x_0p8, y_0p8,'b');


N_0p4 = zeros(1,length(x_0p4));
for i = 1:length(x_0p4)
    N_0p4(i) = trapz(y_0p4(1:length(x_0p4) + 1 - i));
end
[max_0p4, ~] = max(N_0p4);
N_0p4 = N_0p4 * (N0_0p4 / max_0p4);
scatter(x_0p4, N_0p4, 'r');
plot(x_0p4 , N_0p4, 'r');
plot(x_0p4, y_0p4, 'r');

N_1p2 = zeros(1,length(x_1p2));
for i = 1:length(x_1p2)
    N_1p2(i) = trapz(y_1p2(1:length(x_1p2) + 1 - i));
end
[max_1p2, ~] = max(N_1p2);
N_1p2 = N_1p2 * (N0_1p2 / max_1p2);
scatter(x_1p2, N_1p2, 'm');
plot(x_1p2 , N_1p2 , 'm');
plot(x_1p2, y_1p2, 'm');

N_1p6 = zeros(1,length(x_1p2)); %yeah this is screwed 5/21 re-scan the plots 
for i = 1:length(x_1p6)
    N_1p6(i) = trapz(y_1p6(1:length(x_1p6) + 1 - i));
end
[max_1p6, ~] = max(N_1p6);
N_1p6 = N_1p6 * (N0_1p6 / max_1p6);
%scatter(x_1p6, N_1p6, 'g');
%plot(x_1p6 , N_1p6 , 'g');


title('Coverage vs temperature');
ylim([0 1/3]);

%% 4/16 Now doing analysis with interpolated plots. 
%0.078 coverage. 0.4 to 1.2 are at temperatures ~1065, 1162, and 1230 K.
%(5/21 - used 0.8 x 10^19 in figure 7. Corresponds to ~0.2 coverage.)
%Hold on - interpolated vs raw f and S plots line up, but a slight
%disagreement with coverage vs temperature plots. Integral over all signal
%should be proportional to the initial coverage, yeah? 
%4/18 - maybe only take stuff before peak? 




%arrh_x = [1/1065, 1/1162, 1/1230];
%arrh_y = [log(50*0.0382806/24.51), log(50*0.336007/26.39), log(50*0.457/33.6)];
arrh_x = [1/950, 1/1150];
arrh_y = [log(0.004/0.2), log(0.1692/0.2)]; % ONLY MULTIPLY ARGUMENT BY BETA IF RATE IN TERMS OF TIME INSTEAD OF TEMPERATURE 

p = polyfit(arrh_x, arrh_y, 1);

% Extract slope and intercept
slope = p(1);
intercept = p(2);

% Generate fit line
x_fit = linspace(min(arrh_x), max(arrh_x), 100);
y_fit = polyval(p, x_fit);

% Plot
figure(12);
scatter(arrh_x * 10^3, arrh_y, 'filled');
hold on;
plot(x_fit * 10^3, y_fit, 'r-', 'LineWidth', 2);
xlabel('1/T (1/K)  *10^3');
ylabel('ln(beta * f /s)');
title('Arrhenius Plot');
grid on;
legend('Data', 'Linear Fit');
set(gca, 'YDir','reverse')
% Display the equation
fprintf('Fit equation: ln(beta * f/s) = %.4f*(1/T) + %.4f\n', slope, intercept); 
fprintf('Ea in kJ/mol for 0.8 x 10^19 dosage \n = %.4f\n', slope*(-8.314)/1000); % Actual: 160 +- 20 kJ
fprintf('Pre-exponential factor in s^-1 for 0.8 x 10^19 dosage \n = %.4e\n', exp(intercept));% Actual: 2 x 10^7 s^-1
fprintf('Log10() of pre-exponential \n = %.4e\n', log10(exp(intercept))); % Actual: 7.3 +- 0.5 

%% 5/29 TEST - REDUCING SAMPLING RATE 



%5/29 - Signals must be normalized to 1; divide all of the signal by the
%max value of the top initial coverage, in this case the maximum of
%new_y_8p0.

[max_new_y, ~] = max(new_y_8p0);

new_N_0p4 = zeros(1,length(new_x_0p4)); %defines empty array, will be filled in with the areas (the remaining coverage) 

for i = 1:length(new_x_0p4)
    new_N_0p4(i) = trapz(new_y_0p4(1:length(new_x_0p4) + 1 - i)); %gets the area underneath the curve past the signal (new_y) at index 'i'
end
[new_max_0p4, ~] = max(new_N_0p4); %gets maximum of areas for below normalization 
new_N_0p4 = new_N_0p4 * (N0_0p4 / new_max_0p4); %trapz assumes a spacing of 1, so if the actual spacing is lower the reported area's going to be higher
 %more on that, it makes the maximum for new_N_ (the first value)be equal to the initial coverage, N0_

fprintf('NEW INSTANCE \n')
figure(5);
hold on;
%legend('0.4', '0.8', '1.2', '1.6', '2.0', '2.8', '4.0', '8.0');
[max_new_0p4, ~] = max(new_y_0p4);
plot(new_x_0p4 , new_N_0p4, 'o-');
%plot(new_x_0p8, new_y_0p8 .* (N0_0p8/max_new_0p8));
plot(new_x_0p4, new_y_0p4 .* (N0_0p4/max_new_y),'o-');

set(gca, 'ColorOrder', cmap, 'NextPlot', 'replacechildren');

hold off;




newnew_x_0p8 = new_x_0p8(1:2:end);
newnew_y_0p8 = new_y_0p8(1:2:end);

newnew_N_0p8 = zeros(1,length(newnew_x_0p8));
for i = 1:length(newnew_x_0p8)
    newnew_N_0p8(i) = trapz(newnew_y_0p8(1:length(newnew_x_0p8) + 1 - i)); %gets the area underneath the curve past the signal (new_y) at index 'i'
end
[newnew_max_0p8, ~] = max(newnew_N_0p8); %gets maximum of areas for below normalization 
newnew_N_0p8 = newnew_N_0p8 * (N0_0p8 / newnew_max_0p8); %trapz assumes a spacing of 1, so if the actual spacing is lower the reported area's going to be higher
 %more on that, it makes the maximum for new_N_ (the first value)be equal to the initial coverage, N0_

fprintf('NEW INSTANCE \n')
figure(5); 
hold on;
%legend('0.4', '0.8', '1.2', '1.6', '2.0', '2.8', '4.0', '8.0');
[max_new_0p8, ~] = max(newnew_y_0p8);
plot(newnew_x_0p8 , newnew_N_0p8, 'o-');
plot(newnew_x_0p8, newnew_y_0p8 .* (N0_0p8/max_new_y),'o-');

%plot(new_x_0p8, new_y_0p8 .* (N0_0p8/max_new_0p8));

set(gca, 'ColorOrder', cmap, 'NextPlot', 'replacechildren');

hold off;





newnew_x_1p2 = newnew_x_1p2(1:2:end);
newnew_y_1p2 = newnew_y_1p2(1:2:end);
newnew_N_1p2 = zeros(1,length(newnew_x_1p2));
for i = 1:length(newnew_x_1p2)
    newnew_N_1p2(i) = trapz(newnew_y_1p2(1:length(newnew_x_1p2) + 1 - i)); %gets the area underneath the curve past the signal (new_y) at index 'i'
end
[newnew_max_1p2, ~] = max(newnew_N_1p2); %gets maximum of areas for below normalization 
newnew_N_1p2 = newnew_N_1p2 * (N0_1p2 / newnew_max_1p2); %trapz assumes a spacing of 1, so if the actual spacing is lower the reported area's going to be higher
 %more on that, it makes the maximum for new_N_ (the first value)be equal to the initial coverage, N0_

fprintf('NEW INSTANCE \n')
figure(5); 
hold on;
%legend('0.4', '0.8', '1.2', '1.6', '2.0', '2.8', '4.0', '8.0');
[max_new_1p2, ~] = max(newnew_y_1p2);
plot(newnew_x_1p2 , newnew_N_1p2, 'o-');
plot(newnew_x_1p2, newnew_y_1p2 .* (N0_1p2/max_new_y),'o-');

%plot(new_x_1p2, new_y_1p2 .* (N0_1p2/max_new_1p2));

set(gca, 'ColorOrder', cmap, 'NextPlot', 'replacechildren');

hold off;






newnew_x_1p6 = new_x_1p6(1:2:end);
newnew_y_1p6 = new_y_1p6(1:2:end);
newnew_N_1p6 = zeros(1,length(new_x_0p4));


for i = 1:length(newnew_x_1p6)
    newnew_N_1p6(i) = trapz(newnew_y_1p6(1:length(newnew_x_1p6) + 1 - i)); %gets the area underneath the curve past the signal (new_y) at index 'i'
end
[newnew_max_1p6, ~] = max(newnew_N_1p6); %gets maximum of areas for below normalization 
newnew_N_1p6 = newnew_N_1p6 * (N0_1p6 / newnew_max_1p6); %trapz assumes a spacing of 1, so if the actual spacing is lower the reported area's going to be higher
 %more on that, it makes the maximum for new_N_ (the first value)be equal to the initial coverage, N0_

fprintf('NEW INSTANCE \n')
figure(5); 
hold on;
%legend('0.4', '0.8', '1.2', '1.6', '2.0', '2.8', '4.0', '8.0');

plot(newnew_x_1p6 , newnew_N_1p6, 'o-');
plot(newnew_x_1p6, newnew_y_1p6 .* (N0_1p6/max_new_y),'o-');

set(gca, 'ColorOrder', cmap, 'NextPlot', 'replacechildren');

hold off;





newnew_x_2p0 = new_x_2p0(1:2:end);
newnew_y_2p0 = new_y_2p0(1:2:end);
new_N_2p0 = zeros(1,length(newnew_x_2p0));

for i = 1:length(newnew_x_2p0)
    new_N_2p0(i) = trapz(newnew_y_2p0(1:length(newnew_x_2p0) + 1 - i)); %gets the area underneath the curve past the signal (new_y) at index 'i'
end
[new_max_2p0, ~] = max(new_N_2p0); %gets maximum of areas for below normalization 
new_N_2p0 = new_N_2p0 * (N0_2p0 / new_max_2p0); %trapz assumes a spacing of 1, so if the actual spacing is lower the reported area's going to be higher
 %more on that, it makes the maximum for new_N_ (the first value)be equal to the initial coverage, N0_

fprintf('NEW INSTANCE \n')
figure(5); 
hold on;
%legend('0.4', '0.8', '1.2', '1.6', '2.0', '2.8', '4.0', '8.0');

plot(newnew_x_2p0 , new_N_2p0, 'o-');
plot(newnew_x_2p0, newnew_y_2p0 .* (N0_2p0/max_new_y), 'o-');

set(gca, 'ColorOrder', cmap, 'NextPlot', 'replacechildren');

hold off;


