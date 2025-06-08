
%%
N = 64;
n = 0:N-1;
h = [1 2 -3 -3 2 1];
x = cos(2*pi*.1*n);
Y = fft(x,N) .* fft(h,N);
y = ifft(Y);
%%
N = 8 * 1024;
M = 200;
beta = 9; %higher lowers ripple but widens transition band


f = (0:N-1) / N;
f(N/2+1+1:N) = f(N/2+1+1:N)-1;

figure(3); clf; %make filter mask 
Fp_end = 8000/50000;
Fp_start = 4000/50000; 
Fs_start = 0;
Fs_end = 2000/50000;
Fp_dB_start = 6;
Fp_dB_end = 6;
Fp_dB_rip = -1;
Fs_dB_att = -80;


patch([(9/50) .5 .5 (9/50)] ,[-100 -100 Fs_dB_att Fs_dB_att], 0.9 * [1 1 1]);  %comment out if not bandpass
patch([Fs_start Fs_end Fs_end Fs_start] , [-100 -100 Fs_dB_att Fs_dB_att], 0.9 * [1 1 1]);
patch([Fp_start Fp_end Fp_end Fp_start], ...
    [(Fp_dB_start - Fp_dB_rip) (Fp_dB_end - Fp_dB_rip) (Fp_dB_end + Fp_dB_rip) (Fp_dB_start + Fp_dB_rip)],...
    0.9 * [1 1 1]);
grid on;
hold on;
axis([0 0.5 -90 10]);

devi = 0.007; %positive value = grow
Hd = 2.2*(abs(f) > 0.08 - 0.023 & abs(f) < 0.16 + 0.0048); %Numbers respectively are close to start and end of passband, %CHANGE SCALAR TO ADJUST HEIGHT OF PASSBAND 
Hd = Hd .* exp(-j*2*pi*f*(M-1)/2);

%if you want to have a sloped passband; comment out Hd statement if not  
slope = (Fp_dB_start-Fp_dB_end) / (Fp_start - Fp_end); 
%slope = 0.1; %alternatively, hardcode it 
Fp_height = 0.4;
%Hd = Hd .* (10 .^ ((10/ slope) .* (abs(f)- Fp_height)/20 )); 


%plot(f, 20*log10(abs(Hd)), 'r');

hd = ifft(Hd);
figure(8); clf;
n = 0:N-1;
stem(n,hd,'.');
hold on; 

win = kaiser(M,beta)';

h = hd(1:M) .* win; %coefficients 
plot(0:M-1, h, 'r.');

H = fft(h, N );
figure(3);
plot ( f(1:N/2), 20 * log10(abs(H(1:N/2))));
ylim([-100, 20]);

%% START HERE, MAY 2ND 
N = 512;
w = kaiser(N, 7.3);
s = sum(w);
x = x(:); %VERY IMPORTANT 
%plot(0:N-1, abs(fft(x .* w, N)));
%%

N = 64*2048;
plot((0:N-1)*Fs/N, abs(fft(x .* w, N)/(s*2)));
%22.312 kHz, 37.99 V
%21.352 kHz, 0.366866 V

%CODE IN ORDER 
%N = 2048
%plot(0:N-1, abs(fft(x))
%N = 64 * 2048
%w = kaiser(2048, 7.3) 
%stem(w)
%plot(0:N-1, abs(fft(x .* w, N)) 
%s = sum(w) 
%plot((0:N-1)*Fs/N, abs(x .* w, N)/(s*2)) % ZOOM IN