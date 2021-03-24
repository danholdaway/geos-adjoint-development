close all
clear
clc

T = 240:0.1:280;
qi_frac = zeros(1,length(T));

for i = 1:length(T)
    
    if T(i) <= 253.15 
        
        qi_frac(i) = 1;
        
    elseif T(i) > 253.15 && T(i) <= 273.15
        
        qi_frac(i) = -0.05*(T(i) - 273.15);
        
    else
        
        qi_frac(i) = 0.0;
        
    end
    
end

plot(T,qi_frac)