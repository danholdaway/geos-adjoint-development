close all
clear
clc

for i = 0:23

    for j = 0:3
        
        v = [2013, 1, 7, 0, i, 15*j];
        
        a = datestr(v);
        
        if length(a) < 20
            h = '0000';
        else
            h = [a(16:17) a(19:20)]'
        end
        
    end
    
end

% datestr(v)