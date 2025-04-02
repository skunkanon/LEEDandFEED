function [sos,g] = createsos(iir_zeros,iir_poles,k)
zz = iir_zeros;
pp = iir_poles;

[sos, g] = Zp2sos(zz,pp,k,'up','inf');

sections= size(sos,1);
f = linspace(0,0.5,10000);
w = 2*pi*f;
z = exp(j*w);
desired_num = poly(zz)*k;
desired_den = poly(pp);
figure; clf;
desired_H = polyval(desired_num,z) ./ polyval(desired_den,z);
poly(f,20*log10(abs(desired_H)));

%print coefficients for C somewhere here

end
