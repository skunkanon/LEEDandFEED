b = fir1(20,[0.35 0.65]);
freqz(b,0.5,512)


%% 


function [y] = conv(b,x) end


M = length(b)-1;
b = reshape(b,1,M+1);
Lx = length(x);
x = [zeros(K,1);x(:)];
Ly = Lx+M; 
y = zeros(1,Ly);
for n = K+1:Ly
    y(n) = b*x(n:-1:n-M);
end
y = y(K+1:Ly);

