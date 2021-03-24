close all
clear
clc


x = 0:0.5:27;

for i = 1:length(x)

    y(i) = x(i) + rand/2 ;

    if x(i) <= 3
        
        y1(i) = y(i);
        
    else
        
        y1(i) = y(i)+2 + rand;
        
        y(i) = y(i) - 0.3*x(i);
        
    end
    
end
    
plot(x,y)
hold on
plot(x,y1)