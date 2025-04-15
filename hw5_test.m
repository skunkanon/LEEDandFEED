% Rebuilt FIR Low-Pass Filter Design (Window Method, Kaiser)
% Based on example code structure — FIXED

N = 1024;         % FFT length (power of 2, fine)
M = 127;          % Filter length (must be odd for symmetry)
beta = 7.85;      % Kaiser beta for ~80 dB attenuation

% Draw spec boxes (visual only)
figure(1); clf;
patch([0 0.07 0.07 0], [-100 -100 -0.1 -0.1], 0.9*[1 1 1]);  % passband
hold on; grid on;
patch([0.11 0.5 0.5 0.11], [-100 -100 -80 -80], 0.9*[1 1 1]); % stopband
axis([0 0.5 -90 10]);
xlabel('Frequency (cycles/sample)');
ylabel('Filter Gain (dB)');

% Frequency axis: [-0.5, 0.5)
f = (0:N-1)/N;
f(N/2+2:end) = f(N/2+2:end) - 1.0;

% Build desired frequency response (low-pass with cosine transition)
Hd = zeros(1, N);
Hd(abs(f) < 0.07) = 1;  % passband
trans = (abs(f) >= 0.07 & abs(f) <= 0.11);
Hd(trans) = 0.5 * (1 + cos(pi * (abs(f(trans)) - 0.07) / 0.04));

% Linear phase centering
Hd = Hd .* exp(1j * 2 * pi * f * ((M - 1)/2));

% IFFT to get time domain response
hd = real(ifft(Hd));       % impulse response
hd = ifftshift(hd);        % center it around 0

% Extract centered M-tap segment
center = N/2;
start = round(center - floor(M/2));
stop = start + M - 1;
hd_windowed = hd(start:stop);  % this fixes the colon error

% Apply Kaiser window
w = kaiser(M, beta)';
h = hd_windowed .* w;

% Plot impulse response
figure(2); clf;
stem(0:M-1, h, '.');
title('Windowed FIR Impulse Response');
xlabel('n'); ylabel('Amplitude');

% Compute frequency response of resulting filter
H = fft(h, N);

% Plot actual frequency response
figure(1);
plot(f, 20*log10(abs(H) + eps), 'b');
legend('Specs', 'Ideal', 'Actual');

% Diagnostics
passband = abs(f) < 0.07;
stopband = abs(f) > 0.11;

max_pass = max(20*log10(abs(H(passband)) + eps));
min_pass = min(20*log10(abs(H(passband)) + eps));
ripple = max_pass - min_pass;

attenuation = -max(20*log10(abs(H(stopband)) + eps));

fprintf('Passband ripple: %.3f dB\n', ripple);
fprintf('Max passband gain: %.3f dB\n', max_pass);
fprintf('Min passband gain: %.3f dB\n', min_pass);
fprintf('Stopband attenuation: %.2f dB\n', attenuation);
