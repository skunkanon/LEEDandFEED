f = linspace(-0.5,0.5,100); % set f from -0.5 to 0.5
z = exp(2*pi*f*1i); %make independent variable
%%
%filter 1

a1 = -1.5731;
a2 = 0.9025;
b0 = 1;
b1 = 1.6180;
b2 = 1;

H = [];

for i = 1:length(z)

    H = [H,(b0+b1/z(i)+b2/z(i)^2)/(1+a1/z(i)+a2/z(i)^2)];
end

plot (f,H);

clear H

%%


%filter 2

a1 = 0;
a2 = -0.81;
b0 = 1;
b1 = -0.6180;
b2 = 1;

H = [];

for i = 1:length(z)

    H = [H,(b0+b1/z(i)+b2/z(i)^2)/(1+a1/z(i)+a2/z(i)^2)];
end

plot (f,H);

clear H