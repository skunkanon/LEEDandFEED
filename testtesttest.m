
%%

%gui specs 
Fs = 50000;
Fpass = 2000; %centered after mixing, limit is bandwidth 
Fstop = 2333.3333; %3.33 kHz - 1 kHz 
Ap = 0.2;     % dB
As = 62;      % dB

% normalize frequencies
f = [Fpass Fstop]/(Fs/2);
a = [1 0];  %for 0 dB passband, negative infinity dB in stop 
dev = [(10^(Ap/20)-1)/(10^(Ap/20)+1), 10^(-As/20)]; %passband ripple, then stopband attenuation 

% kaiser params 
[N, Wn, beta, ftype] = kaiserord(f, a, dev);

%%

w = kaiser(2048,beta);
wvtool(w);

%%
num_str = sprintf('%.15f, ', Num);               
num_str = ['{', num_str(1:end-2), '}'];         
disp(num_str);                                
%%
clipboard('copy', num_str);  % Copy to clipboard (Windows/macOS)
%%
thing = sum(w);
thingthing = sprintf('%.15f, ', thing);
