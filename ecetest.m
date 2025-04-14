%part 1 spec 
fp = 0.07;        % passband edge
fs = 0.11;        % stopband edge
rip_max_db = 0.1; % max ripple (dB)
min_att_db = 80;  % min attenuation (dB) 

% from log scale to linear 
rip_max = (10^(rip_max_db/20) - 1)/(10^(rip_max_db/20) + 1);
min_att = 10^(-min_att_db/20);
A = -20*log10(min(rip_max, min_att));
df = fs - fp;

% alternatively, can hardcode values for M and beta 
beta = 8;
M = 151;



fprintf("Using M = %d and beta = %.2f\n", M, beta);

% Create ideal impulse response (lowpass)
n = -(M-1)/2 : (M-1)/2;
fc = (fp + fs)/2; % cutoff in middle of transition
hd = 2 * fc * sinc(2 * fc * n);

% Apply Kaiser window
w = kaiser(M, beta)';
h = hd .* w;

% Frequency response
N = 1024;
H = fft(h, N);
f = linspace(0, 0.5, N/2);

% Plot magnitude response
figure(1); clf;
plot(f, 20*log10(abs(H(1:N/2))));
grid on;
xlabel('Frequency (cycles/sample)');
ylabel('Magnitude (dB)');
title('FIR Low-Pass Filter (Window Method)');
ylim([-100 5]);
xline(fp, '--g', 'Passband edge');
xline(fs, '--r', 'Stopband edge');

%%
% Calculate passband ripple and stopband attenuation
H_mag = abs(H(1:N/2));
passband = f <= fp;
stopband = f >= fs;

ripple_db = 20 * log10(max(H_mag(passband)) / min(H_mag(passband)));
attenuation_db = -20 * log10(max(H_mag(stopband)));

fprintf('Passband ripple: %.4f dB\n', ripple_db);
fprintf('Stopband attenuation: %.2f dB\n', attenuation_db);
%%
num_str = sprintf('%.15f, ', Num);               
num_str = ['{', num_str(1:end-2), '}'];         
disp(num_str);                               
clipboard('copy', num_str);  % Copy to clipboard (Windows/macOS)