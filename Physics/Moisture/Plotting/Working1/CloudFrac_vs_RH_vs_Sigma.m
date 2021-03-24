close all
clear
clc

cd /discover/nobackup/drholdaw/ExperimentData/BacmeisterPaper/

% cfsig
% save CF_vs_RH_vs_Sig.mat

load CF_vs_RH_vs_Sig

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/


CF_array = [cf' rh' s'];

[sig_sort, i_sort] = sort(CF_array(:,3),'ascend');

CF_array_sort = CF_array(i_sort,:);

x_min = 0.9;
x_max = 1.1;
q_min = 0.93;
q_max = 1.07;

a = q_max/q_min;
b = (0.5*(a+1)*q_min ) ;

qt = x_min:0.001:x_max;

C2 = zeros(1,length(qt));
for i = 1:length(qt)
    
    if qt(i) < q_min 
        C2(i) = 0;
    elseif q_min <= qt(i) && qt(i) < q_max
        C2(i) = qt(i) * 1/((a-1)*q_min) - 1/(a-1);
    elseif qt(i) >= q_max
        C2(i) = 1;
    end
        
end
C1 = 0.5 + 0.5*tanh( 20*(qt - 1.00)); 

i1  = find(0.9 <= CF_array_sort(:,3) & CF_array_sort(:,3) <= 1.0);
i2  = find(0.8 <= CF_array_sort(:,3) & CF_array_sort(:,3) <  0.9);
i3  = find(0.7 <= CF_array_sort(:,3) & CF_array_sort(:,3) <  0.8);
i4  = find(0.6 <= CF_array_sort(:,3) & CF_array_sort(:,3) <  0.7);
i5  = find(0.5 <= CF_array_sort(:,3) & CF_array_sort(:,3) <  0.6);
i6  = find(0.4 <= CF_array_sort(:,3) & CF_array_sort(:,3) <  0.5);
i7  = find(0.3 <= CF_array_sort(:,3) & CF_array_sort(:,3) <  0.4);
i8  = find(0.2 <= CF_array_sort(:,3) & CF_array_sort(:,3) <  0.3);
i9  = find(0.1 <= CF_array_sort(:,3) & CF_array_sort(:,3) <  0.2);
i10 = find(0.0 <= CF_array_sort(:,3) & CF_array_sort(:,3) <  0.1);

msize = 0.1;

figure
set(gcf,'position',[1281 1 1280 943])
h = scatter(CF_array(i1,2),CF_array(i1,1),'.');
hold on
plot(qt,C1,'r','LineWidth',2.0)
plot(qt,C2,'g','LineWidth',2.0)
xlim([x_min x_max])
hChildren = get(h, 'Children');
set(hChildren, 'Markersize', msize)

asdad

figure
set(gcf,'position',[1281 1 1280 943])
h = scatter(CF_array(i2,2),CF_array(i2,1),'.');
hold on
plot(qt,C1,'r','LineWidth',2.0)
plot(qt,C2,'g','LineWidth',2.0)
xlim([x_min x_max])
hChildren = get(h, 'Children');
set(hChildren, 'Markersize', msize)

asd

figure
set(gcf,'position',[1281 1 1280 943])
h = scatter(CF_array(i3,2),CF_array(i3,1),'.');
hold on
plot(qt,C1,'r','LineWidth',2.0)
plot(qt,C2,'g','LineWidth',2.0)
xlim([x_min x_max])
hChildren = get(h, 'Children');
set(hChildren, 'Markersize', msize)

figure
set(gcf,'position',[1281 1 1280 943])
h = scatter(CF_array(i4,2),CF_array(i4,1),'.');
hold on
plot(qt,C1,'r','LineWidth',2.0)
plot(qt,C2,'g','LineWidth',2.0)
xlim([x_min x_max])
hChildren = get(h, 'Children');
set(hChildren, 'Markersize', msize)

figure
set(gcf,'position',[1281 1 1280 943])
h = scatter(CF_array(i5,2),CF_array(i5,1),'.');
hold on
plot(qt,C1,'r','LineWidth',2.0)
plot(qt,C2,'g','LineWidth',2.0)
xlim([x_min x_max])
hChildren = get(h, 'Children');
set(hChildren, 'Markersize', msize)

figure
set(gcf,'position',[1281 1 1280 943])
h = scatter(CF_array(i6,2),CF_array(i6,1),'.');
hold on
plot(qt,C1,'r','LineWidth',2.0)
plot(qt,C2,'g','LineWidth',2.0)
xlim([x_min x_max])
hChildren = get(h, 'Children');
set(hChildren, 'Markersize', msize)

figure
set(gcf,'position',[1281 1 1280 943])
h = scatter(CF_array(i7,2),CF_array(i7,1),'.');
hold on
plot(qt,C1,'r','LineWidth',2.0)
plot(qt,C2,'g','LineWidth',2.0)
xlim([x_min x_max])
hChildren = get(h, 'Children');
set(hChildren, 'Markersize', msize)

figure
set(gcf,'position',[1281 1 1280 943])
h = scatter(CF_array(i8,2),CF_array(i8,1),'.');
hold on
plot(qt,C1,'r','LineWidth',2.0)
plot(qt,C2,'g','LineWidth',2.0)
xlim([x_min x_max])
hChildren = get(h, 'Children');
set(hChildren, 'Markersize', msize)

figure
set(gcf,'position',[1281 1 1280 943])
h = scatter(CF_array(i9,2),CF_array(i9,1),'.');
hold on
plot(qt,C1,'r','LineWidth',2.0)
plot(qt,C2,'g','LineWidth',2.0)
xlim([x_min x_max])
hChildren = get(h, 'Children');
set(hChildren, 'Markersize', msize)

figure
set(gcf,'position',[1281 1 1280 943])
h = scatter(CF_array(i10,2),CF_array(i10,1),'.');
hold on
plot(qt,C1,'r','LineWidth',2.0)
plot(qt,C2,'g','LineWidth',2.0)
xlim([x_min x_max])
hChildren = get(h, 'Children');
set(hChildren, 'Markersize', msize)