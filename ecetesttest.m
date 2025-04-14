N = 1024; 
M =31; 
beta = 5;

figure(1); clf; %freq response
patch([0 .2 .2 0], [-100, -100, -80, -80], 0.9*[1 1 1]);
hold on; grid on; 
patch([0.45 0.5 0.5 0.5] , [-100 -100 -80 -80], 0.9*[1 1 1]);
patch([0.25 0.35 0.35 0.25], [-11, -1, +1 -9], 0.9*[1 1 1]);
axis([0, 0.5, -90, 10]);
xlabel('Frequency (cycles/sample)');
ylabel('Filter Gain (dB)');

f = (0:N-1)/N;
f(N/2+1+1:end) = f(N/2+1+1:end) - 1.0;

Hd = (abs(f)>0.225 & abs(f)<0.4) .* 10.^( (10/0.1)*(abs(f)-0.35)/20);
Hd = Hd .* exp(-j*2*pi*f*(M-1)/2);

figure(1);
plot(f, 20 * log10(abs(Hd)), 'r');

hd = ifft(Hd);
figure(2); clf;
n = 0:N-1;
stem(n,hd,'.');
hold on;

win = kaiser(M,beta);
h = hd(1:M) .* win;
stem(0:M-1,h,'r.');

H = fft(h,N);
figure(1);
plot(f, 20*log10(abs(H)));