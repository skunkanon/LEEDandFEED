n = linspace(0,20,21);

function y = u(x)
    y = zeros(size(x));  
    y(x > 0) = 1;      
    y(x == 0) = 1;       
end

x_1 = [];

for i = 1:length(n)

    x_1 = [x_1,2*sqrt(67)/3*(0.8)^n(i)*cos(0.723*n(i)-1.32)*u(n(i))];
end

x_2 = [];

for i = 1:length(n)

    x_2 = [x_2,(4/3)*((0.8)^n(i)*cos(0.723*n(i)))*u(n(i))+2*sqrt(7)*((0.8)^n(i)*sin(0.723*n(i)))*u(n(i))];
end

x_2 = [];

for i = 1:length(n)

    x_2 = [x_2,(4/3)*((0.8)^n(i)*cos(0.723*n(i))*u(i))+2*sqrt(7)*((0.8)^n(i)*sin(0.723*n(i)))*u(n(i))];
end

x_3 = [];


for i = 1:length(n)

    x_3 = [x_3,2.52*(0.8^(n(i)+1)*sin(0.723*n(i)+1)*u(n(i)+1))-1.2*(0.8^(n(i))*sin(0.723*n(i))*u(n(i)))+2*sqrt(7)*(0.8^n(i)*sin(0.723*n(i))*(u(n(i))))];
end




plot (n,x_1,n,x_2,'b--o',n,x_3,'k*');