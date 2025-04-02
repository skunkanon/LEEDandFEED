fprintf('NEW INSTANCE \n');
%% 1) POLES AND ZEROES


N = 1024;
w = 2*pi()*[0:N-1]/N;
z = exp(j*w);
minz = -50; maxz = 50;
iir_poles = [.8*exp(j*2*pi()*0.22*[-1,1])
             .8*exp(j*2*pi()*0.18*[-1,1])
             .6*exp(j*2*pi()*0.080*[-1,1])
             .6*exp(j*2*pi()*0.080*[-1,1])];

iir_zeros = [-1*exp(j*2*pi()*0.350*[-1,1])
             -1*exp(j*2*pi()*0.350*[-1,1])
             1*exp(j*2*pi()*0.250*[-1,1])
             1*exp(j*2*pi()*0.750*[-1,1])];

figure(1); clf
plot(z);
hold on
plot(real(iir_poles), imag(iir_poles), 'x');
plot(real(iir_zeros),imag(iir_zeros),'o');
axis square
axis([-1.5 1.5 -1.5 1.5]);
grid 
xlabel('Re(z)');
ylabel('Im(z)');
title('Unit Circle with Poles and Zeroes');
%% 2) FREQUENCY RESPONSE 

iir_poles = iir_poles(:);
iir_zeros = iir_zeros(:);

dem = poly(iir_poles);
num = poly(iir_zeros);

figure(2); clf
f= linspace(-0.5,0.5, N);
w = 2*pi()*f;
z = exp(j*w);
H = polyval(num, z) ./ polyval(dem, z);
Hdb = 20 * log10(abs(H));
plot(f,Hdb);
grid
xlabel('Normalized Frequency f(cycles/sample)');
ylabel('Magnitude');
title('Frequency Response');
axis([0 0.5 -100 10]);
%% 3) 3D PLOT 

x = linspace(-1, 1, 101); 
y = linspace(-1, 1, 101); 
[real_x, imag_y] = meshgrid(x, y); 
z_grid = real_x + j*imag_y; 

H_grid = polyval(num, z_grid(:)) ./ polyval(dem, z_grid(:)); 
H_grid = reshape(H_grid, size(real_x)); 
Hdb_grid = 20 * log10(abs(H_grid)); 
Hdb_grid(Hdb_grid < minz) = nan;

figure(3); clf;
surf(real_x, imag_y, Hdb_grid, 'EdgeColor', 'none');
xlabel('Re(z)');
ylabel('Im(z)');
zlabel('Magnitude (dB)');
title('3D Frequency Response');
axis([-1 1 -1 1 minz maxz]);
colorbar;
grid on;
%% 4) CREATE.SOS INTEGRATION

%some error around here

num = poly(iir_zeros);
den = poly(iir_poles);


max_gain = filternorm(num,den,inf);
k = 1/max_gain;
createsos(iir_zeros,iir_poles,k);

