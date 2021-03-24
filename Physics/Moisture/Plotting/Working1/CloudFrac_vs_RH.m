close all
clear
clc

cd /discover/nobackup/drholdaw/ExperimentData/BacmeisterPaper/

load CloudFrac_vs_RH
qsx_cf = qsx;
qt_cf = qt_in;

load LSCloud_vs_RH
qsx_ql = qsx;
qt_ql = qt_in;

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

xmin = 0.8;
xmax = 1.2;

qt = xmin:0.0001:xmax;



C2 = zeros(1,length(qt));

q_min = 0.93;
q_max = 1.07;
a = q_max/q_min;
b = (0.5*(a+1)*q_min ) ;

for i = 1:length(qt)
    
    if qt(i) < q_min 
        
        C2(i) = 0;
        
    elseif q_min <= qt(i) && qt(i) < q_max
        
        C2(i) = qt(i) * 1/((a-1)*q_min) - 1/(a-1);

    elseif qt(i) >= q_max
        
        C2(i) = 1;
                
    end
        
end


C1 = 0.5 + 0.5*tanh( 15*(qt - 1.0015)); 



scatter(qt_cf(1:100:end)./qsx_cf(1:100:end),CF(1:100:end),'kx')
xlim([xmin xmax])
hold on
plot(qt,C1,'r','LineWidth',2.5)
plot(qt,C2,'g','LineWidth',2.5)


qdiff = -2e-3:1e-6:4e-3;
q1 = zeros(1,length(qdiff));

for i = 1:length(qdiff)
    
    if qdiff(i) < 0 
        
        q1(i) = 0;
        
    else
        
        q1(i) = 1.0*qdiff(i);
                
    end
        
end




figure
scatter(qt_ql(1:100:end)-qsx_ql(1:100:end),qc(1:100:end),'k.')
hold on
plot(qdiff,q1,'r','LineWidth',1.5)
% plot(qdiff,q2,'g','LineWidth',1.5)