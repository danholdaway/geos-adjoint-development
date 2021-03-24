close all
clear
clc

cd /discover/nobackup/drholdaw/tmp.22292/sens.20130117.000000/
% fg
Qprogress_NoCloudUpdate

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/


[a,b] = size(QIC);

for i = 1:b

    QIL(a+1,i) = max(QIL(:,i));
    QLL(a+1,i) = max(QLL(:,i));
    QIC(a+1,i) = max(QIC(:,i));
    QLC(a+1,i) = max(QLC(:,i));
    
end

% plot(QIL(b+1,:))
hold on
plot(QLL(b+1,:),'r')
% plot(QIC(b+1,:),'g--')
plot(QLC(b+1,:),'k')

figure

for j = 1:96

% plot(QIL(j,:))
hold on
plot(QLL(j,:),'r')
% plot(QIC(j,:),'g--')
plot(QLC(j,:),'k')

end

xlabel('Time Step')
xlim([1 96])