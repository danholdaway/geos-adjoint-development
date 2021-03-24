close all
clear
clc

cd /discover/nobackup/drholdaw/ExperimentData/BacmeisterPaper/

% qcrh1
% save LSCloud_vs_QTmQs.mat

% load LSCloud_vs_RH.mat
load LSCloud_vs_QTmQs.mat

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

QC_array = [qc' rh' s'];

[sig_sort, i_sort] = sort(QC_array(:,3),'ascend');

QC_array_sort = QC_array(i_sort,:);

i1  = find(0.9 <= QC_array_sort(:,3) & QC_array_sort(:,3) <= 1.0);
i2  = find(0.8 <= QC_array_sort(:,3) & QC_array_sort(:,3) <  0.9);
i3  = find(0.7 <= QC_array_sort(:,3) & QC_array_sort(:,3) <  0.8);
i4  = find(0.6 <= QC_array_sort(:,3) & QC_array_sort(:,3) <  0.7);
i5  = find(0.5 <= QC_array_sort(:,3) & QC_array_sort(:,3) <  0.6);
i6  = find(0.4 <= QC_array_sort(:,3) & QC_array_sort(:,3) <  0.5);
i7  = find(0.3 <= QC_array_sort(:,3) & QC_array_sort(:,3) <  0.4);
i8  = find(0.2 <= QC_array_sort(:,3) & QC_array_sort(:,3) <  0.3);
i9  = find(0.1 <= QC_array_sort(:,3) & QC_array_sort(:,3) <  0.2);
i10 = find(0.0 <= QC_array_sort(:,3) & QC_array_sort(:,3) <  0.1);

msize = 0.1;
x_min = -1.3e-3;
x_max = 0.5e-3;

y_min = 0;
y_max = 0.5e-3;

figure
set(gcf,'position',[1281 1 1280 943])
h = scatter(QC_array(i1,2),QC_array(i1,1),'.');
xlim([x_min x_max])
ylim([y_min y_max])
hChildren = get(h, 'Children');
set(hChildren, 'Markersize', msize)

figure
set(gcf,'position',[1281 1 1280 943])
h = scatter(QC_array(i2,2),QC_array(i2,1),'.');
xlim([x_min x_max])
ylim([y_min y_max])
hChildren = get(h, 'Children');
set(hChildren, 'Markersize', msize)

figure
set(gcf,'position',[1281 1 1280 943])
h = scatter(QC_array(i3,2),QC_array(i3,1),'.');
xlim([x_min x_max])
ylim([y_min y_max])
hChildren = get(h, 'Children');
set(hChildren, 'Markersize', msize)

figure
set(gcf,'position',[1281 1 1280 943])
h = scatter(QC_array(i4,2),QC_array(i4,1),'.');
xlim([x_min x_max])
ylim([y_min y_max])
hChildren = get(h, 'Children');
set(hChildren, 'Markersize', msize)

figure
set(gcf,'position',[1281 1 1280 943])
h = scatter(QC_array(i5,2),QC_array(i5,1),'.');
xlim([x_min x_max])
ylim([y_min y_max])
hChildren = get(h, 'Children');
set(hChildren, 'Markersize', msize)

figure
set(gcf,'position',[1281 1 1280 943])
h = scatter(QC_array(i6,2),QC_array(i6,1),'.');
xlim([x_min x_max])
ylim([y_min y_max])
hChildren = get(h, 'Children');
set(hChildren, 'Markersize', msize)

figure
set(gcf,'position',[1281 1 1280 943])
h = scatter(QC_array(i7,2),QC_array(i7,1),'.');
xlim([x_min x_max])
ylim([y_min y_max])
hChildren = get(h, 'Children');
set(hChildren, 'Markersize', msize)

figure
set(gcf,'position',[1281 1 1280 943])
h = scatter(QC_array(i8,2),QC_array(i8,1),'.');
xlim([x_min x_max])
ylim([y_min y_max])
hChildren = get(h, 'Children');
set(hChildren, 'Markersize', msize)

figure
set(gcf,'position',[1281 1 1280 943])
h = scatter(QC_array(i9,2),QC_array(i9,1),'.');
xlim([x_min x_max])
ylim([y_min y_max])
hChildren = get(h, 'Children');
set(hChildren, 'Markersize', msize)

figure
set(gcf,'position',[1281 1 1280 943])
h = scatter(QC_array(i10,2),QC_array(i10,1),'.');
xlim([x_min x_max])
ylim([y_min y_max])
hChildren = get(h, 'Children');
set(hChildren, 'Markersize', msize)
