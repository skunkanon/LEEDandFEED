function [town,car] = equationfit_test_feb25(crown,vic)

town = crown * 1997;
car = vic * 2011; 

%{
x = linspace(0,10,100); %input, 0 to 100
shifttest = rand();
y_raw = cos(shifttest*x); %raw cos(x)
y_actual = []; %noisy cos(x)
for i = 1:length(y_raw)
    y_actual = [y_actual,y_raw(i) - rand() + rand()]; %fake noise on y_raw

end



model = @(params,x) params(1) * cos(params(2) * x + params(3));

initial_params = [1,1,0];

params = lsqcurvefit(model, initial_params, x, y_actual);

disp('fitted parameters');
disp(params);

fitted_actual = model(params,x);
plot(x,y_actual,x,fitted_actual,'b--o');


A =[1,2,30,4];
[M,I] = max(A);
%}

end
%% 
g = @(C) (C/C_0)*(C_0 - C + f*C_0);

C_0 = 1;
f = 0.9;
c_span = linspace(0,1);
g_span = [];

for i = 1:length(c_span)
    g_span = [g_span,g(c_span(i))];
end
plot(flip(c_span),g_span);
%% R = 8.314; %gas constant in J/(K*mol)
%N_0 = 0.6; %initial surface 
kcal_to_j = 4184; %conversion factor

Ea_0 = 32.6 * kcal_to_j;
%Ea_0 = 60 * 10^3;
v = 8.5 * 10^15;

%coverage-dependent activation energy
y_E = 4 * kcal_to_j; %4kcal/mol into j/mol, eq) 14 
Ea = @(C) Ea_0 - C*y_E;



%defining experimental conditions
beta = 9; % Kelvin/s
init_tmp = 300; %initial temperature
ds_time = 20; %define desorption duration in seconds
t_span = linspace(0,ds_time,10000); %array of times from start through whole duration

%temperature as function of time 
tmp = @(t) beta*t + init_tmp; %function declaration
tmp_span = tmp(t_span); %independent variable for spectra spanning duration 

%C_0 = 1;
f = 0.99;
g = @(C) (C/C_0)*(C_0 - C + f*C_0);

dCdt = @(t,C) -v * g(C) * exp(-Ea(C)/(R*tmp(t))); % = -N, the rate of desorption
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9); % Adjust tolerances
[t,C] = ode23s(dCdt, t_span,C_0,options);


dCdt_span = arrayfun(@(t,C) dCdt(t,C), t,C); %IS NEGATIVE, returned multiplied by -1
%Ea_span = arrayfun(@(C) Ea(C),C);
%g_span = arrayfun(@(C) g(C), C);

rate_span = -dCdt_span;

tmp_C_span = tmp(t_span) - 273;
%figure
%plot(tmp_C_span, -dCdt_span);
%legend('coverage', 'rate');
%% 

town = {1,2,3
        4,5,6
        7,8,9};
car = town{1,2};
%% 

f = [0    0.1   0.1   0.35  0.35  0.5];  

H_dB = [-45     -45    0    0    -55   -55];

figure;
plot(f, H_dB, 'r-o','LineWidth',2); 
grid on;
xlabel('Normalized Frequency  f = F / F_s');
ylabel('Magnitude (dB)');
title('Ideal Filter');
axis([0 0.5 -60 5]);

%% 
%% Automated Bandpass Filter Design Using fdesign
% We define the filter so that the passband is between 0.1 and 0.35 (normalized)
% and the stopband requirements are met in the lower and upper stopbands.
% Note: Frequencies here are defined with Fs = 1 (so Nyquist = 0.5).

Fs = 1;  % Sampling frequency (normalized so that 0.5 is Nyquist)

% Define filter specification parameters:
%   Fst1: end of lower stopband
%   Fp1:  beginning of passband
%   Fp2:  end of passband
%   Fst2: beginning of upper stopband
%   Ast1: minimum stopband attenuation in lower stopband (in dB)
%   Ap:   maximum passband ripple (in dB)
%   Ast2: minimum stopband attenuation in upper stopband (in dB)
%
% Our verbal spec is that the filter should be “off” (deep attenuation) for
% f < 0.1 and f > 0.35. For instance, we set:
Fst1 = 0.05;  % Lower stopband edge (before passband begins)
Fp1  = 0.1;   % Lower passband edge (transition here from stopband to passband)
Fp2  = 0.35;  % Upper passband edge (transition from passband to stopband)
Fst2 = 0.5;   % Upper stopband edge (Nyquist)
Ast1 = 40;    % At least 40 dB attenuation in the lower stopband
Ap   = 1;     % Allow up to 1 dB ripple in the passband
Ast2 = 50;    % At least 50 dB attenuation in the upper stopband

% Create the fdesign object. Here the frequencies are in Hz (with Fs = 1).
d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
    Fst1, Fp1, Fp2, Fst2, Ast1, Ap, Ast2, Fs);

% Design the filter using the elliptic method (which typically gives the minimum order).
Hd = design(d, 'ellip');

%% Frequency Response and Specification Mask Plot

% Compute the frequency response of the designed filter.
Nfft = 1024;
[fResp, H] = freqz(Hd, Nfft, Fs);  % fResp in Hz (0 to 0.5)
HdB = 20*log10(abs(H));

figure; clf;
plot(fResp, HdB, 'b-', 'LineWidth', 1.5);
grid on; hold on;
xlabel('Normalized Frequency (F/Fs)');
ylabel('Magnitude (dB)');
title('Bandpass Filter Response with Spec Mask');
axis([0 0.5 -100 10]);

% Create the specification mask.
% Here we define a piecewise mask that, for example, is intended to show that:
% - In the lower stopband (0 to 0.1) the filter should be below about -45 dB.
% - In the passband (0.1 to 0.35) the filter should be near 0 dB.
% - In the upper stopband (0.35 to 0.5) the filter should be below about -55 dB.
%
% (This mask gives a little margin relative to the verbal spec of -40 dB at 0.1
% and -50 dB at 0.35.)
f_mask = [0    0.1   0.1   0.35  0.35   0.5];  
H_mask = [-45  -45    0     0    -55   -55];

plot(f_mask, H_mask, 'r-o', 'LineWidth',2);
legend('Designed Filter','Spec Mask','Location','best');
%% 
fprintf('NEW INSTANCE \n');

% Specifications (frequencies are normalized so that 1 corresponds to Nyquist)
Ap   = 1;     % Maximum passband ripple in dB
Ast1 = 40;    % Minimum attenuation in the lower stopband in dB
Ast2 = 50;    % Minimum attenuation in the upper stopband in dB

% Frequency edge choices:
Fst1 = 0.08;  % Lower stopband edge (for f < 0.1, we require at least 40 dB attenuation)
Fp1  = 0.1;   % Lower passband edge
Fp2  = 0.35;  % Upper passband edge
Fst2 = 0.37;  % Upper stopband edge (for f > 0.35, we require at least 50 dB attenuation)

% Create the filter design specification object.
d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
                     Fst1, Fp1, Fp2, Fst2, Ast1, Ap, Ast2);

% Design the filter using the elliptic method.
Hd = design(d, 'ellip');

% Visualize the magnitude response.
fvtool(Hd);

%% 

Rp = 1;              % +- this in the passband
As = 50;             % min attenuation in stopband
Wp = [0.15, 0.25];    % passband edge
Ws = [0.1, 0.35];   % stopband edge 

[n, Wn] = ellipord(Wp, Ws, Rp, As);

[z, p,k] = ellip(n, Rp, As, Wn, 'bandpass');

%fvtool(b, a);

[sos,g] = zp2sos(z,p,k);


%% 
fprintf('NEW INSTANCE \n');

fprintf('ZEROS \n');

for i=1:length(z)
    %fprintf('%.4f\n',abs(z(i)));
    fprintf('%.4f',abs(z(i)));
    fprintf(',');
    fprintf('%.4f\n',angle(z(i))/(2*pi()));
end

fprintf('POLES \n');

for i=1:length(p)
    fprintf('%.4f',abs(p(i)));
    fprintf(',');
    fprintf('%.4f\n',angle(p(i)/(2*pi())));
end


%%

for i=1:length(Np_index)
    fprintf('%.4',Np_index(i))
end

