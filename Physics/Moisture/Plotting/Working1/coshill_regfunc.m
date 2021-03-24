close all
clear
clc

D = 0.5;
x = 0:0.01:1;
d = 0.2;
n = 2;

f = zeros(1,length(x));

for j = 1:length(x)
    if x(j) <= D+d && x(j) >= D-d
        f(j) = (cos(pi*(x(j)-D)/(2*d)))^n;
    end
end
    


plot(x,f)