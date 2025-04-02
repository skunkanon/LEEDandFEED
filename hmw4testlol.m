fprintf('NEW INSTANCE \n');

%getting data from text file
directory = 'C:\Users\vifro\OneDrive\Documents\MATLAB'; 
file_name = 'hw4data.txt';
file_path = fullfile(directory, file_name);
x = load(file_path); 
x = x(:); 



%%
Fs = 75 * 1000; %75 ksps = sampling frequency 
N = length(x); %is 512 

NFFT = 32768;

f = linspace(0, Fs, NFFT);

beta = 10;

figure;
hold on;


    win = kaiser(N,beta);
    x_win = x .* win;

    fft_x_win = fft(x_win,NFFT);
    fft_x_win_dB = 20 * log10( abs(fft_x_win) + 0.0000000000001);

    
    plot(f, fft_x_win_dB);




xlim([8000 13000]);
ylim([-100 60]);

hold off;

xlabel('Frequency (Hz)');
ylabel('Magnitude (dB');
title('hw4data.txt Kaiser Window, Beta = 10');

%%calculate where gain is maximum, getting frequency 

ind_large = find(f >9000 & f < 10000);
[M_large, I_large] = max(fft_x_win_dB(ind_large));
f_large = f(ind_large(I_large));



ind_small = find(f > 11000 & f < 12000);
[M_small, I_small] = max(fft_x_win_dB(ind_small));
f_small = f(ind_small(I_small));

%% calculate amplitude in V

fprintf('NEW INSTANCE \n');

%multiply by 2 because symmetric
%divide by window size to account for tapering by window 

A_large = 2 * abs(fft_x_win(ind_large(I_large))) / sum(win);

A_small = 2 * abs(fft_x_win(ind_small(I_small))) / sum(win);

%%

fprintf('Large Component Frequency, Hz = %.4f\n',f_large);
fprintf('Large Component Amplitude, V = %.4f\n', A_large);

fprintf('Small Component Frequency, Hz = %.4f\n',f_small);
fprintf('Small Component Amplitude, V = %.4f\n', A_small);


%% Testing out different beta = 0, see if same result (it's the same, mostly) 

%getting data from text file
directory = 'C:\Users\vifro\OneDrive\Documents\MATLAB'; 
file_name = 'hw4data.txt';
file_path = fullfile(directory, file_name);
x = load(file_path); 
x = x(:); 




Fs = 75 * 1000; %75 ksps = sampling frequency 
N = length(x); %is 512 

NFFT = 32768;

f = linspace(0, Fs, NFFT);

beta = 0;

figure;
hold on;


    win = kaiser(N,beta);
    x_win = x .* win;

    fft_x_win = fft(x_win,NFFT);
    fft_x_win_dB = 20 * log10( abs(fft_x_win) + 0.0000000000001);

    
    plot(f, fft_x_win_dB);




xlim([8000 13000]);
ylim([-100 60]);

hold off;

xlabel('Frequency (Hz)');
ylabel('Magnitude (dB');
title('hw4data.txt Kaiser Window, Beta = 0');

%calculate where gain is maximum, getting frequency 

ind_large = find(f >9000 & f < 10000);
[M_large, I_large] = max(fft_x_win_dB(ind_large));
f_large = f(ind_large(I_large));



ind_small = find(f > 11000 & f < 12000);
[M_small, I_small] = max(fft_x_win_dB(ind_small));
f_small = f(ind_small(I_small));

% calculate amplitude in V

fprintf('NEW INSTANCE \n');

%multiply by 2 because symmetric
%divide by window size to account for tapering by window 

A_large = 2 * abs(fft_x_win(ind_large(I_large))) / sum(win);

A_small = 2 * abs(fft_x_win(ind_small(I_small))) / sum(win);



fprintf('Beta = 0, Large Component Frequency, Hz = %.4f\n',f_large);
fprintf('Beta = 0, Large Component Amplitude, V = %.4f\n', A_large);

fprintf('Beta = 0, Small Component Frequency, Hz = %.4f\n',f_small);
fprintf('Beta = 0, Small Component Amplitude, V = %.4f\n', A_small);
%%
% Task 2: Front-end Implementation - Filter Design for Receiver Front-End
% This filter is intended to extract the baseband (after mixing) signal.
% The original band is 11 kHz ±1.3 kHz, and after mixing it shifts to baseband (~ -1.3 to +1.3 kHz).
% We design a lowpass filter with the following specifications:
%   - Sampling Rate: 50 ksps
%   - Passband edge: 1.5 kHz (to comfortably include the desired signal)
%   - Stopband edge: 2.5 kHz (to reject signals that would alias after decimation by 10)
%   - Passband ripple: ±0.5 dB maximum
%   - Stopband attenuation: at least 65 dB

clear; close all; clc;

% Parameters
fs = 50e3;         % Sampling rate: 50 ksps
fp = 1.5e3;        % Passband edge (Hz)
fsb = 2.5e3;       % Stopband edge (Hz)
Rp = 0.5;          % Passband ripple (dB)
Ast = 65;          % Stopband attenuation (dB)

% Design a lowpass IIR filter using the elliptic design method
lpFilt = designfilt('lowpassiir', ...
    'PassbandFrequency', fp, ...
    'StopbandFrequency', fsb, ...
    'PassbandRipple', Rp, ...
    'StopbandAttenuation', Ast, ...
    'DesignMethod', 'ellip', ...
    'SampleRate', fs);

% Display the designed filter coefficients and order information
disp(lpFilt);

% Plot the magnitude response using fvtool for an interactive view
fvtool(lpFilt, 'Fs', fs, 'MagnitudeDisplay', 'Magnitude (dB)');
title('Lowpass Filter Magnitude Response (Elliptic Design)');

% Alternatively, plot using freqz for a static plot:
Nfft = 2^14; % FFT length for smooth curve
[H, f] = freqz(lpFilt, Nfft, fs);
figure;
plot(f, 20*log10(abs(H)), 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Designed Lowpass Filter Frequency Response');
grid on;
xlim([0 5e3]); % Focus on frequencies up to 5 kHz
