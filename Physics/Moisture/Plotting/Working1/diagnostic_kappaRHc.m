close all
clear
clc


sigma = 0:0.01:1;

kappa = zeros(1,length(sigma));
kappa1 = zeros(1,length(sigma));
RHc = zeros(1,length(sigma));
RHc1 = zeros(1,length(sigma));

for i = 1:length(sigma)
   
    if sigma(i) >= 0.2
        
        kappa(i)  = 0.9* ( sigma(i) - 0.2 )^0.2;
        kappa1(i) = 0.90*( sigma(i) - 0.2 )^0.10;
%         kappa1(i) = 0.80*( sigma(i) - 0.2 )^0.10;

    else
        
        kappa(i) = 0.0;
        kappa1(i) = 0.0;
        
    end
       
    RHc(i) =  0.85 - 0.7* sigma(i).*(1 - sigma(i))*(1.85 + 0.95*(sigma(i) - 0.5));
    RHc1(i) = 0.85 - 0.40*sigma(i).*(1 - sigma(i))*(1.85 + 0.95*(sigma(i) - 0.5));

end


plot(kappa,sigma)
hold on
plot(kappa1,sigma,'r')
set(gca,'YDir','reverse')
xlim([0 1])

figure
plot(RHc,sigma)
hold on
plot(RHc1,sigma,'r')
set(gca,'YDir','reverse')
xlim([0 1])
