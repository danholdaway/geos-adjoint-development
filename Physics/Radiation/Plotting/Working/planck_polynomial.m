close all
clear
clc
mydir = pwd;

cb = [ 5.3443e+0,  -2.0617e-1,   2.5333e-3, ...
      -6.8633e-6,   1.0115e-8,  -6.2672e-12, ...
       2.7148e+1,  -5.4038e-1,   2.9501e-3, ...
       2.7228e-7,  -9.3384e-9,   9.9677e-12, ...
      -3.4860e+1,   1.1132e+0,  -1.3006e-2, ...
       6.4955e-5,  -1.1815e-7,   8.0424e-11, ...
      -6.0513e+1,   1.4087e+0,  -1.2077e-2, ...
       4.4050e-5,  -5.6735e-8,   2.5660e-11, ...
      -2.6689e+1,   5.2828e-1,  -3.4453e-3, ...
       6.0715e-6,   1.2523e-8,  -2.1550e-11, ...
      -6.7274e+0,   4.2256e-2,   1.0441e-3, ...
      -1.2917e-5,   4.7396e-8,  -4.4855e-11, ...
       1.8786e+1,  -5.8359e-1,   6.9674e-3, ...
      -3.9391e-5,   1.0120e-7,  -8.2301e-11, ...
       1.0344e+2,  -2.5134e+0,   2.3748e-2, ...
      -1.0692e-4,   2.1841e-7,  -1.3704e-10, ...
      -1.0482e+1,   3.8213e-1,  -5.2267e-3, ...
       3.4412e-5,  -1.1075e-7,   1.4092e-10, ...
       1.6769e+0,   6.5397e-2,  -1.8125e-3, ...
       1.2912e-5,  -2.6715e-8,   1.9792e-11];
   
cb = reshape(cb,6,10);
   
T = 250:0.01:300;

i = length(T);
   
xlayer = zeros(length(T),10);
xlayerp = zeros(length(T),10);
xlayer_lin = zeros(length(T),10);

for ibn = 1:10

    for i = 1:length(T)

        xlayer(i,ibn) = T(i)*(T(i)*(T(i)*(T(i)*(T(i)*cb(6, ibn)+cb(5, ibn))+cb(4, ibn))+cb(3, ibn))+ cb(2, ibn)) + cb(1, ibn);
        xlayerp(i,ibn) = 1*(T(i)*(T(i)*(T(i)*(T(i)*cb(6, ibn)+cb(5, ibn))+cb(4, ibn))+cb(3, ibn)) ...
                         +cb(2, ibn)) + T(i)*(1*(T(i)*(T(i)*(T(i)*cb(6, ibn)+cb(5, ibn))+cb(4, ibn))+cb(3, ibn)) ...
                        +T(i)*(1*(T(i)*(T(i)*cb(6, ibn)+cb(5, ibn))+cb(4, ibn))+T(i)*(1*(T(i)*cb (6, ibn)+cb(5, ibn))+T(i)*cb(6, ibn)*1)));
        
    end
    
end

for ibn = 1:10

    for i = 1:length(T)

        m = (max(xlayer(:,ibn)) - min(xlayer(:,ibn)))/(max(T)-min(T));
        xlayer_lin(i,ibn) = m * T(i) - m * T(1) + min(xlayer(:,ibn));

    end
    
end

for ibn = 1:1

    figure
    plot(T,xlayer(:,ibn))
    hold on
    plot(T,xlayer_lin(:,ibn),'r')
    
    figure
    plot(T,xlayerp(:,ibn))
%     hold on
%     plot(T,xlayer_lin(:,ibn),'r')
end
