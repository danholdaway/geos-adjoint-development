close all
clear
clc

theta(1) = 300;

theta_p(1) = 300+0.4;

theta_tl(1) = 0.01;

t = 1:20;

for i = 2:length(t)
    
    theta(i) =  theta(i-1) + 150*sin(i*pi/3)+300;
    
    theta_p(i) =  theta_p(i-1) + theta_p(i-1)*sind(i);
   
    theta_tl(i) = theta_tl(i-1) - sind(theta(i-1))*theta_tl(i-1);
    
end


plot(t,theta)
hold on
% plot(t,theta_p,'r')