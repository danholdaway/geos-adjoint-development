close all
clear
clc

q = 0:10^-6:10^-3;

qs = 3*10^-4;
c = 1.1;

qscrit = c*qs;

a = 1.1;

C1 = zeros(1,length(q));
C2 = zeros(1,length(q));

b = (0.5*(c+1)*qs ) ;

for i = 1:length(q)
    
    if q(i) < qs 
        
        C1(i) = 0;
        
    elseif qs <= q(i) && q(i) < qscrit
        
        C1(i) = q(i) * 1/((a-1)*qs) - 1/(a-1);

    elseif q(i) >= qscrit
        
        C1(i) = 1;
                
    end
    
    
    C2(i) = 0.5 + 0.5*tanh( 60000*(q(i) - b));
    
end


plot(q,C1)
hold on
plot(q,C2,'r')