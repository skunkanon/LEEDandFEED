%%
N = 8 * 1024;
M = 100;
beta = 8; %higher lowers ripple but widens transition band


f = (0:N-1) / N;
f(N/2+1+1:N) = f(N/2+1+1:N)-1;

figure(1); clf; %make filter mask 
Fp_end = 0.07;
Fp_start = 0; 
Fs_start = 0.11;
Fs_end = 0.5;
Fp_dB_start = 0;
Fp_dB_end = 0;
Fp_dB_rip = -0.1;
Fs_dB_att = -80;


%patch([Fp_start Fp_end Fp_end Fp_start] ,[-100 -100 Fs_dB_att Fs_dB_att], 0.9 * [1 1 1]);  %comment out if not bandpass
patch([Fs_start Fs_end Fs_end Fs_start] , [-100 -100 Fs_dB_att Fs_dB_att], 0.9 * [1 1 1]);
patch([Fp_start Fp_end Fp_end Fp_start], ...
    [(Fp_dB_start - Fp_dB_rip) (Fp_dB_end - Fp_dB_rip) (Fp_dB_end + Fp_dB_rip) (Fp_dB_start + Fp_dB_rip)],...
    0.9 * [1 1 1]);
grid on;
hold on;
axis([0 0.5 -90 10]);

Hd = (abs(f) > 0 & abs(f) < 0.08); %Numbers respectively are close to start and end of passband 
Hd = Hd .* exp(-j*2*pi*f*(M-1)/2);

%if you want to have a sloped passband; comment out Hd statement if not  
slope = (Fp_dB_start-Fp_dB_end) / (Fp_start - Fp_end); 
%slope = 0.1; %alternatively, hardcode it 
Fp_height = 0.4;
%Hd = Hd .* (10 .^ ((10/ slope) .* (abs(f)- Fp_height)/20 )); 


%plot(f, 20*log10(abs(Hd)), 'r');

hd = ifft(Hd);
figure(2); clf;
n = 0:N-1;
stem(n,hd,'.');
hold on; 

win = kaiser(M,beta)';

h = hd(1:M) .* win; %coefficients 
plot(0:M-1, h, 'r.');

H = fft(h, N );
figure(1);
plot ( f(1:N/2), 20 * log10(abs(H(1:N/2))));
