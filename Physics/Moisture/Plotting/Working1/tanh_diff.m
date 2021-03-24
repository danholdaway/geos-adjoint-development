close all
clear
clc


xmin = 0.8;
xmax = 1.2;
qt = xmin:0.0001:xmax;

qoff = 0.0665;
q_min = 1-qoff;
q_max = 1+qoff;
grad = 1/(q_max - q_min);
q1 = 3*grad;

CF = 0.5 + 0.5*tanh( q1*(qt - 1.0)); 

dCFdqt = 0.5 * q1 * sech(q1*(qt - 1.0)).^2



a = q_max/q_min;
b = (0.5*(a+1)*q_min ) ;

CFlin = zeros(1,length(qt));
dCFlindqt = zeros(1,length(qt));
for i = 1:length(qt)
    
    if qt(i) < q_min 
        CFlin(i) = 0;
        dCFlindqt(i) = 0;
    elseif q_min <= qt(i) && qt(i) < q_max
        CFlin(i) = qt(i) * 1/((a-1)*q_min) - 1/(a-1);
        dCFlindqt(i) = 1/((a-1)*q_min);
    elseif qt(i) >= q_max
        CFlin(i) = 1;
        dCFlindqt(i) = 0;
    end
        
end


reg1 = 0.66*( cosh(10*(qt-1.0)).^-2);


dCFdqt_reg = dCFdqt .* reg1;

figure
plot(qt,CF)
hold on
plot(qt,CFlin,'r')

figure
plot(qt,dCFlindqt)
hold on
plot(qt,dCFdqt,'r')
plot(qt,reg1,'g')
plot(qt,dCFdqt_reg,'k')

