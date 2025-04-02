%% (b) Plot DTFT of 512-sample Kaiser windows
%clear; clc; close all;

N = 512;                      % Window length
beta_values = [0 2 4 6 8];    % Kaiser window shape parameters
NFFT = 16*1024;               % Large FFT size (16k) to approximate DTFT

figure('Name','DTFT of Kaiser Windows','NumberTitle','off');
hold on; grid on; box on;
colors = lines(length(beta_values));  % Distinct plot colors
legend_entries = cell(size(beta_values));

for i = 1:length(beta_values)
    beta = beta_values(i);
    w = kaiser(N, beta);      % 512-sample Kaiser window
    
    % Compute the "DTFT" via an NFFT-point FFT
    Wf = fft(w, NFFT);        
    % Frequency axis in cycles/sample, from 0 to 1
    f = (0:NFFT-1)/NFFT;      
    
    % Convert magnitude to dB, adding a small epsilon to avoid log(0)
    Wf_dB = 20*log10(abs(Wf)+0.00000000000000001);
    
    % Plot only the low-frequency region: 0 <= f < 0.01
    fmax = 0.01; 
    idx = f < fmax;  % indices for frequencies up to 0.01
    plot(f(idx), Wf_dB(idx), 'LineWidth', 1.5, ...
         'Color', colors(i,:), ...
         'DisplayName', sprintf('\\beta = %g', beta));
    
    legend_entries{i} = sprintf('\\beta = %g', beta);
end

xlabel('Normalized Frequency (cycles/sample)');
ylabel('Magnitude (dB)');
title('DTFT (in dB) of 512-Sample Kaiser Windows');
legend('show','Location','Best');
xlim([0 0.01]);    % Zoom in to the specified frequency range
ylim([-100 5]);    % Adjust dB range as desired
hold off;

%% 
% Generate Kaiser windows and compute their DTFTs
betas = [0, 2, 4, 6, 8, 10];
N_FFT = 16384; % 16k FFT
f = (0:N_FFT-1)/N_FFT; % Frequency axis in cycles/sample

figure;
for i = 1:length(betas)
    w = kaiser(512, betas(i)); % Generate Kaiser window
    W = fft(w, N_FFT); % Compute FFT
    W_mag_dB = 20*log10(abs(W)); % Convert to dB
    
    subplot(2, 3, i);
    plot(f(1:N_FFT/200), W_mag_dB(1:N_FFT/200)); % Plot 0 < f < 0.01
    title(['DTFT, \beta = ', num2str(betas(i))]);
    xlabel('Frequency (cycles/sample)');
    ylabel('Magnitude (dB)');
    grid on;
    ylim([-100 20]); % Adjust dB range for clarity
end
%% 2b (All plots on one figure)
samp = 512;
beta = [0, 2, 4, 6, 8, 10];
NFFT = 16*1024; % 16384-point FFT
colors = lines(length(beta)); % Distinct colors for each beta

figure;
hold on;

for i = 1:length(beta)
    % Generate Kaiser window and compute FFT
    w = kaiser(samp, beta(i));
    w_fft = fft(w, NFFT);
    f = linspace(0, 1, NFFT); 
    
    % Convert to dB (avoid log(0))
    w_fft_dB = 20 * log10(abs(w_fft) + 1e-9); 
    
    % Plot for 0 < f < 0.01 cycles/sample
    idx = f <= 0.01; 
    plot(f(idx), w_fft_dB(idx), 'LineWidth', 1.5, ...
        'Color', colors(i,:), 'DisplayName', ['\beta = ', num2str(beta(i))]);
end

hold off;

% Formatting
xlim([0, 0.01]);
ylim([-100, 20]);
xlabel('Frequency (cycles/sample)');
ylabel('Magnitude (dB)');
title('DTFT of Kaiser Windows (0 < f < 0.01)');
legend('show', 'Location', 'best');
grid on;
%% Problem 3(a): DTFT of Windowed Cosine Signal
%clear; clc; close all;

% Parameters
A = 3.7;          % Amplitude (V)
f0 = 0.3308;      % Normalized frequency (cycles/sample)
N = 512;          % Number of samples
NFFT = 32768;     % FFT size for DTFT approximation

% Generate x_w[n] = A*cos(2πf0n) for n = 0,1,...,511
n = 0:N-1;        % Time indices
x_w = A * cos(2*pi*f0*n); % Windowed cosine signal

% Compute DTFT via 32768-point FFT
Xw = fft(x_w, NFFT); 
Xw_mag_dB = 20*log10(abs(Xw) + 1e-9); % Magnitude in dB (avoid log(0))

% Frequency axis (0 to 1 cycles/sample)
f = linspace(0, 1, NFFT);

% Plot
figure;
plot(f, Xw_mag_dB, 'LineWidth', 1.5);
xlabel('Discrete-Time Frequency, f (cycles/sample)');
ylabel('Magnitude (dB)');
title('DTFT of x_w[n] (NFFT = 32768)');
grid on;
xlim([0 1]);      % Full frequency range
ylim([-20 50]);   % Adjusted for visibility

% Highlight leakage around f0
hold on;
plot([f0 f0], ylim, '--r', 'LineWidth', 1.2, 'DisplayName', 'f_0 = 0.3308');
legend('show', 'Location', 'northeast');
hold off;
%%
% Define parameters
A = 3.7;          % Amplitude in volts
f0 = 0.3308;      % Normalized frequency
N = 512;          % Number of samples in x_w[n]
Nfft = 32768;     % FFT size for detailed spectrum

% Generate the signal x_w[n] for n = 0 to 511
n = 0:N-1;
x_w = A * cos(2 * pi * f0 * n);

% Compute the 32768-point FFT of x_w[n]
X_w = fft(x_w, Nfft);

% Compute the magnitude in decibels (dB)
mag_dB = 20 * log10(abs(X_w));

% Generate the frequency axis (normalized frequency)
f = (0:Nfft-1)/Nfft;

% Extract the positive frequency range (0 to 0.5)
f_plot = f(1:Nfft/2 + 1);
mag_dB_plot = mag_dB(1:Nfft/2 + 1);

% Create the plot
figure;
plot(f_plot, mag_dB_plot);
xlabel('Discrete-time frequency f');
ylabel('Magnitude (dB)');
title('Spectrum of x_w[n] using 32768-point FFT');
grid on;
%% 
%% Problem 3(b): DTFT vs. 512-point DFT Comparison
clear; clc; close all;

% Parameters
A = 3.7;          % Amplitude (V)
f0 = 0.3308;      % Normalized frequency (cycles/sample)
N = 512;          % Number of samples
NFFT = 32768;     % FFT size for DTFT approximation

% Generate x_w[n] = A*cos(2πf0n) for n = 0,1,...,511
n = 0:N-1;        
x_w = A * cos(2*pi*f0*n);

% Compute DTFT (32768-point FFT)
Xw_DTFT = fft(x_w, NFFT);
Xw_DTFT_dB = 20*log10(abs(Xw_DTFT) + 1e-9); % Magnitude in dB
f_DTFT = linspace(0, 1, NFFT);              % DTFT frequency axis

% Compute 512-point DFT
Xw_DFT = fft(x_w);                          % 512-point FFT
Xw_DFT_dB = 20*log10(abs(Xw_DFT) + 1e-9);  % Magnitude in dB
f_DFT = (0:N-1)/N;                         % DFT frequency axis

% Plot DTFT and DFT
figure;
plot(f_DTFT, Xw_DTFT_dB, 'LineWidth', 1.5, 'DisplayName', 'DTFT (32768-point FFT)');
hold on;
plot(f_DFT, Xw_DFT_dB, 'ro', 'MarkerSize', 5, 'DisplayName', '512-point DFT');
xlim([0.30 0.35]);                         % Zoom to 0.30 < f < 0.35
ylim([0 60]);                              % Adjusted for amplitude visibility

% Theoretical peak amplitude (A*N/2 in linear scale)
peak_expected_linear = A * N / 2;          % Expected peak magnitude
peak_expected_dB = 20*log10(peak_expected_linear); % ~59.53 dB

% Annotate actual amplitude and frequency
plot([f0 f0], ylim, '--k', 'LineWidth', 1.2, 'DisplayName', 'f_0 = 0.3308');
plot(xlim, [peak_expected_dB peak_expected_dB], '--m', 'LineWidth', 1.2, ...
    'DisplayName', 'Expected Peak (59.53 dB)');

xlabel('Discrete-Time Frequency, f (cycles/sample)');
ylabel('Magnitude (dB)');
title('DTFT vs. 512-point DFT of x_w[n] (0.30 < f < 0.35)');
legend('show', 'Location', 'northeast');
grid on;
hold off;
%%
%% Problem 3(a) and 3(b): DTFT and DFT Comparison
A = 3.7;          % Amplitude (V)
f0 = 0.3308;      % Normalized frequency (cycles/sample)
N = 512;          % Number of samples
NFFT = 32768;     % FFT size for DTFT approximation

% Generate x_w[n] = A*cos(2πf0n) for n = 0,1,...,511
n = linspace(0, N-1, N);
x_w = A * cos(2*pi*f0*n);

% Compute DTFT (32768-point FFT)
fft_x_w = fft(x_w, NFFT);
fft_x_w_dB = 20 * log10(abs(fft_x_w) + 1e-13); % Magnitude in dB (avoid log(0))

% Frequency axis for DTFT
f_DTFT = linspace(0, 1, NFFT);

% Compute 512-point DFT
fft_x_w_DFT = fft(x_w);                      % 512-point FFT
fft_x_w_DFT_dB = 20 * log10(abs(fft_x_w_DFT) + 1e-13); % Magnitude in dB
f_DFT = (0:N-1)/N;                           % DFT frequency axis

% Plot DTFT and DFT
figure;
hold on;

% Plot DTFT
plot(f_DTFT, fft_x_w_dB, 'LineWidth', 1.5, 'DisplayName', 'DTFT (32768-point FFT)');

% Plot DFT points
plot(f_DFT, fft_x_w_DFT_dB, 'ro', 'MarkerSize', 5, 'DisplayName', '512-point DFT');

% Highlight f0
plot([f0 f0], ylim, '--r', 'LineWidth', 1.2, 'DisplayName', 'f_0 = 0.3308');

% Theoretical peak amplitude (A*N/2 in linear scale)
peak_expected_linear = A * N / 2;            % Expected peak magnitude
peak_expected_dB = 20*log10(peak_expected_linear); % ~59.53 dB

% Annotate expected peak amplitude
plot(xlim, [peak_expected_dB peak_expected_dB], '--m', 'LineWidth', 1.2, ...
    'DisplayName', 'Expected Peak (59.53 dB)');

% Formatting
xlim([0.30 0.35]);                           % Zoom to 0.30 < f < 0.35
ylim([0 60]);                                % Adjusted for amplitude visibility
xlabel('Discrete-Time Frequency, f (cycles/sample)');
ylabel('Magnitude (dB)');
title('DTFT vs. 512-point DFT of x_w[n] (0.30 < f < 0.35)');
legend('show', 'Location', 'northeast');
grid on;
hold off;

%% 3b, copied 3a code
fprintf('NEW INSTANCE \n');

A = 3.7;
f0 = 0.3308;
N = 512; 
NFFT = 32768;
n = linspace(0, N-1, N);
x_w = A * cos(2*pi*f0*n);

% Compute DTFT (32768-point FFT)
fft_x_w = fft(x_w, NFFT);
fft_x_w_dB = 20 * log10(abs(fft_x_w) + 1e-13); % Magnitude in dB

% Frequency axis for DTFT
f_DTFT = linspace(0, 1, NFFT);

% Plot DTFT
figure;
hold on;
plot(f_DTFT, fft_x_w_dB, 'LineWidth', 1.5, 'DisplayName', 'DTFT (32768-point FFT)');
ylim([-20 60]);
xlim([0.30 0.35]);
plot([f0 f0], ylim, '--r', 'DisplayName', 'f_0 = 0.3308'); % Highlight f0

% Compute 512-point DFT
fft_x_w_DFT = fft(x_w); % 512-point FFT
fft_x_w_DFT_dB = 20 * log10(abs(fft_x_w_DFT) + 1e-13); % Magnitude in dB

% Frequency axis for DFT (using linspace)
f_DFT = linspace(0, 1, N); % Correct frequency axis for DFT

% Plot DFT points
plot(f_DFT, fft_x_w_DFT_dB, 'ro', 'MarkerSize', 6, 'DisplayName', '512-point DFT');

% Formatting
ylim([-20 60]);
xlim([0.30 0.35]);
xlabel('Discrete-Time Frequency, f (cycles/sample)');
ylabel('Magnitude (dB)');
title('DTFT vs. 512-point DFT of x_w[n] (0.30 < f < 0.35)');
legend('show', 'Location', 'northeast');
grid on;
hold off;

%% Problem 3(c): DTFT and DFT with Kaiser Window (β=8)
%clear; clc; close all;

% Parameters
A = 3.7;          % Amplitude (V)
f0 = 0.3308;      % Normalized frequency (cycles/sample)
N = 512;          % Number of samples
NFFT = 32768;     % FFT size for DTFT approximation
beta = 8;         % Kaiser window parameter

% Generate original signal x[n]
n = 0:N-1;
x = A * cos(2*pi*f0*n);

% Generate Kaiser window and apply it
w = kaiser(N, beta);        % 512-point Kaiser window (β=8)
x_w_windowed = x .* w';     % Windowed signal (element-wise multiplication)

% Compute DTFT (32768-point FFT)
Xw_DTFT = fft(x_w_windowed, NFFT);
Xw_DTFT_dB = 20*log10(abs(Xw_DTFT) + 1e-13); % Magnitude in dB

% Compute 512-point DFT
Xw_DFT = fft(x_w_windowed);                % 512-point FFT
Xw_DFT_dB = 20*log10(abs(Xw_DFT) + 1e-13); % Magnitude in dB

% Frequency axes
f_DTFT = linspace(0, 1, NFFT);            % DTFT frequency axis
f_DFT = (0:N-1)/N;                        % DFT frequency axis

% Calculate expected peak using coherent gain
coherent_gain = sum(w)/N;                 % Coherent gain of the window
peak_expected_linear = A * coherent_gain * N/2; % Expected peak magnitude
peak_expected_dB = 20*log10(peak_expected_linear);

% Plot
figure;
hold on;
plot(f_DTFT, Xw_DTFT_dB, 'LineWidth', 1.5, 'DisplayName', 'DTFT (32768-point FFT)');
plot(f_DFT, Xw_DFT_dB, 'ro', 'MarkerSize', 5, 'DisplayName', '512-point DFT');
plot([f0 f0], ylim, '--r', 'LineWidth', 1.2, 'DisplayName', 'f_0 = 0.3308');
plot(xlim, [peak_expected_dB peak_expected_dB], '--m', 'LineWidth', 1.2, ...
    'DisplayName', ['Expected Peak (', num2str(peak_expected_dB, '%.1f'), ' dB)']);

% Formatting
xlim([0.30 0.35]);
ylim([-100 60]); % Adjusted to show reduced side-lobes
xlabel('Discrete-Time Frequency, f (cycles/sample)');
ylabel('Magnitude (dB)');
title('DTFT vs. 512-point DFT of Windowed x_w[n] (β=8, 0.30 < f < 0.35)');
legend('show', 'Location', 'northeast');
grid on;
hold off;
%% Problem 4: Frequency and Amplitude Analysis
%clear; clc; close all;

% =============================================
% Step 1: Load Data from hw4data.txt
% =============================================
current_folder = 'C:\Users\vifro\OneDrive\Documents\MATLAB'; 
file_name = 'hw4data.txt';
file_path = fullfile(current_folder, file_name);

% Load data
if exist(file_path, 'file') == 2
    x = load(file_path); 
    x = x(:); % Ensure column vector
    fprintf('Data loaded: %d samples\n', length(x));
else
    error('File "%s" not found in folder: %s', file_name, current_folder);
end

% =============================================
% Step 2: Signal Parameters
% =============================================
Fs = 75000;           % Sampling frequency (75 ksps)
N = length(x);        % 512 samples
NFFT = 32768;         % Zero-padded FFT size for high resolution

% =============================================
% Step 3: Apply Window and Compute FFT
% =============================================
win = hann(N);        % Hann window to reduce leakage
x_win = x .* win;     % Windowed signal

% Compute FFT
X_win = fft(x_win, NFFT);
X_win_mag = abs(X_win); 
X_win_dB = 20*log10(X_win_mag + eps); % Avoid log(0)

% Frequency axis (Hz)
f = (0:NFFT-1)*(Fs/NFFT); 

% =============================================
% Step 4: Find Peaks with Interpolation
% =============================================
% --- Larger signal (9-10 kHz) ---
f_range1 = [9000, 10000];
idx1 = find(f >= f_range1(1) & f <= f_range1(2));
[~, max_idx1] = max(X_win_dB(idx1));
% Quadratic interpolation around the peak
if max_idx1 > 1 && max_idx1 < length(idx1)
    y = X_win_dB(idx1(max_idx1-1:max_idx1+1));
    p = polyfit((-1:1)', y, 2);
    delta = -p(2)/(2*p(1));
    f1 = f(idx1(max_idx1)) + delta*(Fs/NFFT);
else
    f1 = f(idx1(max_idx1));
end

% --- Smaller signal (11-12 kHz) ---
f_range2 = [11000, 12000];
idx2 = find(f >= f_range2(1) & f <= f_range2(2));
[~, max_idx2] = max(X_win_dB(idx2));
% Quadratic interpolation around the peak
if max_idx2 > 1 && max_idx2 < length(idx2)
    y = X_win_dB(idx2(max_idx2-1:max_idx2+1));
    p = polyfit((-1:1)', y, 2);
    delta = -p(2)/(2*p(1));
    f2 = f(idx2(max_idx2)) + delta*(Fs/NFFT);
else
    f2 = f(idx2(max_idx2));
end

% =============================================
% Step 5: Calculate Amplitudes (Volts Peak)
% =============================================
coherent_gain = sum(win); % Hann window compensation
A1 = 2 * interp1(f, X_win_mag, f1, 'spline') / coherent_gain; 
A2 = 2 * interp1(f, X_win_mag, f2, 'spline') / coherent_gain;

% =============================================
% Step 6: Plot and Results
% =============================================
figure;
plot(f, X_win_dB, 'LineWidth', 1.5);
xlim([8000 13000]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Spectrum of Windowed Signal (Hann Window)');
grid on;
hold on;
plot(f1, interp1(f, X_win_dB, f1, 'spline'), 'ro', 'MarkerSize', 8, 'DisplayName', 'Larger Signal');
plot(f2, interp1(f, X_win_dB, f2, 'spline'), 'go', 'MarkerSize', 8, 'DisplayName', 'Smaller Signal');
legend('show');

% Display results
fprintf('=============================================\n');
fprintf('Larger signal frequency: \t%.1f Hz (Target: 9-10 kHz)\n', f1);
fprintf('Smaller signal frequency: \t%.1f Hz (Target: 11-12 kHz)\n', f2);
fprintf('Larger signal amplitude: \t%.2f V (peak)\n', A1);
fprintf('Smaller signal amplitude: \t%.2f V (peak)\n', A2);
fprintf('=============================================\n');

%%

fprintf('NEW INSTANCE \n');
% =============================================
% Step 1: Load Data from hw4data.txt
% =============================================
current_folder = 'C:\Users\vifro\OneDrive\Documents\MATLAB'; 
file_name = 'hw4data.txt';
file_path = fullfile(current_folder, file_name);

% Load data
if exist(file_path, 'file') == 2
    x = load(file_path); 
    x = x(:); % Ensure column vector
    fprintf('Data loaded: %d samples\n', length(x));
else
    error('File "%s" not found in folder: %s', file_name, current_folder);
end


% =============================================
% Step 2: Signal Parameters
% =============================================

Fs = 75000;           % Sampling frequency (75 ksps)
N = length(x);        % 512 samples
NFFT = 32768;         % Zero-padded FFT size for high resolution

% =============================================
% Step 3: Apply Kaiser Window and Compute FFT
% =============================================
beta = 10;             % Kaiser window parameter (controls sidelobe suppression)
win = kaiser(N, beta); % Generate Kaiser window
x_win = x .* win;     % Windowed signal

% Compute FFT
X_win = fft(x_win, NFFT);
X_win_mag = abs(X_win); 
X_win_dB = 20*log10(X_win_mag + eps); % Avoid log(0)

% Frequency axis (Hz)
f = (0:NFFT-1)*(Fs/NFFT); 

% =============================================
% Step 4: Find Peaks with Interpolation
% =============================================
% --- Larger signal (9-10 kHz) ---
f_range1 = [9000, 10000];
idx1 = find(f >= f_range1(1) & f <= f_range1(2));
[~, max_idx1] = max(X_win_dB(idx1));
f1 = f(idx1(max_idx1)); % Directly use bin frequency

% --- Smaller signal (11-12 kHz) ---
f_range2 = [11000, 12000];
idx2 = find(f >= f_range2(1) & f <= f_range2(2));
[~, max_idx2] = max(X_win_dB(idx2));
f2 = f(idx2(max_idx2)); % Directly use bin frequency

% =============================================
% Step 5: Calculate Amplitudes (No Spline Interpolation)
% =============================================
coherent_gain = sum(win); % Window compensation factor
A1 = 2 * X_win_mag(idx1(max_idx1)) / coherent_gain; % Use peak bin magnitude
A2 = 2 * X_win_mag(idx2(max_idx2)) / coherent_gain;

% =============================================
% Step 6: Plot and Results
% =============================================
figure;
plot(f, X_win_dB, 'LineWidth', 1.5);
xlim([8000 13000]); ylim([-100 60]);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Spectrum (Kaiser Window, \beta=10)');
grid on; hold on;
plot(f1, X_win_dB(idx1(max_idx1)), 'ro', 'MarkerSize', 8, 'DisplayName', 'Larger Signal');
plot(f2, X_win_dB(idx2(max_idx2)), 'go', 'MarkerSize', 8, 'DisplayName', 'Smaller Signal');
legend('show', 'Location', 'northeast');

% Display results
fprintf('=============================================\n');
fprintf('Larger signal frequency: \t%.1f Hz (Accuracy: ±%.1f Hz)\n', f1, Fs/NFFT);
fprintf('Smaller signal frequency: \t%.1f Hz (Accuracy: ±%.1f Hz)\n', f2, Fs/NFFT);
fprintf('Larger signal amplitude: \t%.2f V (peak)\n', A1);
fprintf('Smaller signal amplitude: \t%.2f V (peak)\n', A2);
fprintf('=============================================\n');
%% Problem 4: Frequency and Amplitude Analysis with Multiple Kaiser Windows
clear; clc; close all;

% =============================================
% Step 1: Load Data from hw4data.txt
% =============================================
current_folder = 'C:\Users\vifro\OneDrive\Documents\MATLAB'; 
file_name = 'hw4data.txt';
file_path = fullfile(current_folder, file_name);

% Load data
if exist(file_path, 'file') == 2
    x = load(file_path); 
    x = x(:); % Ensure column vector
    fprintf('Data loaded: %d samples\n', length(x));
else
    error('File "%s" not found in folder: %s', file_name, current_folder);
end

% =============================================
% Step 2: Signal Parameters
% =============================================
Fs = 75000;           % Sampling frequency (75 ksps)
N = length(x);        % 512 samples
NFFT = 32768;         % Zero-padded FFT size for high resolution

% Frequency axis (Hz)
f = (0:NFFT-1)*(Fs/NFFT); 

% =============================================
% Step 3: Apply Kaiser Windows and Compute FFTs
% =============================================
betas = [0, 2, 4, 6, 8, 10]; % Kaiser window parameters
colors = lines(length(betas)); % Distinct colors for each beta

figure;
hold on;

for i = 1:length(betas)
    beta = betas(i);
    win = kaiser(N, beta); % Generate Kaiser window
    x_win = x .* win;     % Windowed signal

    % Compute FFT
    X_win = fft(x_win, NFFT);
    X_win_dB = 20*log10(abs(X_win) + eps); % Avoid log(0)

    % Plot spectrum
    plot(f, X_win_dB, 'LineWidth', 1.5, 'DisplayName', ['\beta = ', num2str(beta)]);
end

% =============================================
% Step 4: Format Plot
% =============================================
xlim([8000 13000]);
ylim([-100 60]); % Adjusted for visibility
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Spectrum of Windowed Signal (Kaiser Window, Multiple \beta)');
grid on;
legend('show', 'Location', 'northeast');
hold off;

%% Problem 4: Frequency and Amplitude Analysis with Averaging
clear; clc; close all;

% =============================================
% Step 1: Load Data from hw4data.txt
% =============================================
current_folder = 'C:\Users\vifro\OneDrive\Documents\MATLAB'; 
file_name = 'hw4data.txt';
file_path = fullfile(current_folder, file_name);

% Load data
if exist(file_path, 'file') == 2
    x = load(file_path); 
    x = x(:); % Ensure column vector
    fprintf('Data loaded: %d samples\n', length(x));
else
    error('File "%s" not found in folder: %s', file_name, current_folder);
end

% =============================================
% Step 2: Signal Parameters
% =============================================
Fs = 75000;           % Sampling frequency (75 ksps)
N = length(x);        % 512 samples
NFFT = 32768;         % Zero-padded FFT size for high resolution

% Frequency axis (Hz) using linspace
f = linspace(0, Fs, NFFT); 

% =============================================
% Step 3: Apply Kaiser Windows and Compute FFTs
% =============================================
betas = [2, 4, 6, 8, 10]; % Kaiser window parameters
colors = lines(length(betas)); % Distinct colors for each beta

% Arrays to store results
f1_results = zeros(1, length(betas)); % Larger signal frequencies
f2_results = zeros(1, length(betas)); % Smaller signal frequencies
A1_results = zeros(1, length(betas)); % Larger signal amplitudes
A2_results = zeros(1, length(betas)); % Smaller signal amplitudes

figure;
hold on;

for i = 1:length(betas)
    beta = betas(i);
    win = kaiser(N, beta); % Generate Kaiser window
    x_win = x .* win;     % Windowed signal

    % Compute FFT
    X_win = fft(x_win, NFFT);
    X_win_mag = abs(X_win); 
    X_win_dB = 20*log10(X_win_mag + eps); % Avoid log(0)

    % Plot spectrum
    plot(f, X_win_dB, 'LineWidth', 1.5, 'DisplayName', ['\beta = ', num2str(beta)]);

    % =============================================
    % Step 4: Find Peaks and Calculate Amplitudes
    % =============================================
    % --- Larger signal (9-10 kHz) ---
    f_range1 = [9000, 10000];
    idx1 = find(f >= f_range1(1) & f <= f_range1(2));
    [~, max_idx1] = max(X_win_dB(idx1));
    f1 = f(idx1(max_idx1));

    % --- Smaller signal (11-12 kHz) ---
    f_range2 = [11000, 12000];
    idx2 = find(f >= f_range2(1) & f <= f_range2(2));
    [~, max_idx2] = max(X_win_dB(idx2));
    f2 = f(idx2(max_idx2));

    % Calculate amplitudes (corrected for window gain)
    coherent_gain = sum(win); % Kaiser window compensation
    A1 = 2 * X_win_mag(idx1(max_idx1)) / coherent_gain; 
    A2 = 2 * X_win_mag(idx2(max_idx2)) / coherent_gain;

    % Store results
    f1_results(i) = f1;
    f2_results(i) = f2;
    A1_results(i) = A1;
    A2_results(i) = A2;

    % Display results for this beta
    fprintf('=============================================\n');
    fprintf('Kaiser Window (β = %d)\n', beta);
    fprintf('Larger signal frequency: \t%.1f Hz (Target: 9-10 kHz)\n', f1);
    fprintf('Smaller signal frequency: \t%.1f Hz (Target: 11-12 kHz)\n', f2);
    fprintf('Larger signal amplitude: \t%.2f V (peak)\n', A1);
    fprintf('Smaller signal amplitude: \t%.2f V (peak)\n', A2);
    fprintf('=============================================\n');
end

% =============================================
% Step 5: Calculate Averages
% =============================================
f1_avg = mean(f1_results); % Average larger signal frequency
f2_avg = mean(f2_results); % Average smaller signal frequency
A1_avg = mean(A1_results); % Average larger signal amplitude
A2_avg = mean(A2_results); % Average smaller signal amplitude

% Display averaged results
fprintf('=============================================\n');
fprintf('Averaged Results Across All β Values:\n');
fprintf('Larger signal frequency: \t%.1f Hz (Target: 9-10 kHz)\n', f1_avg);
fprintf('Smaller signal frequency: \t%.1f Hz (Target: 11-12 kHz)\n', f2_avg);
fprintf('Larger signal amplitude: \t%.2f V (peak)\n', A1_avg);
fprintf('Smaller signal amplitude: \t%.2f V (peak)\n', A2_avg);
fprintf('=============================================\n');

% =============================================
% Step 6: Format Plot
% =============================================
xlim([8000 13000]);
ylim([-100 60]); % Adjusted for visibility
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Spectrum of Windowed Signal (Kaiser Window, Multiple \beta)');
grid on;
legend('show', 'Location', 'northeast');
hold off;