close all
clear
clc

cd /discover/nobackup/drholdaw/ExperimentData/BacmeisterPaper/

% Variables
% Jacobian
% clear i
% save Jacobian.mat
% save Variables.mat

load Jacobian.mat
load Variables.mat

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

im = length(V);

im90 = round(im/100 * 90);
im80 = round(im/100 * 80);
im70 = round(im/100 * 70);
im60 = round(im/100 * 60);
im50 = round(im/100 * 50);

percentile = zeros(5,8);
sortJ = zeros(im,8);

% figure
for i = 1:8

    [sortJ(:,i), sorti] = sort(J(:,i),'ascend');
    sortV = V(sorti,:);
    percentile(1,i) = sortJ(im90);
    percentile(2,i) = sortJ(im80);
    percentile(3,i) = sortJ(im70);
    percentile(4,i) = sortJ(im60);
    percentile(5,i) = sortJ(im50);
    
%     subplot(2,4,i)
%     plot(sortJ(:,i))
end


plot_vj = 5;
plot_jac = 1;
plot_inc = 1:1:im;

msize = 0.1;

figure
set(gcf,'position',[1281 1 1280 943])
h = scatter(V(:,plot_vj),J(:,plot_vj),'.');
xlims = get(gca,'xlim');
% hold on
% plot(-5:5,ones(1,11)*percentile(1,plot_vj),'r')
% plot(-5:5,ones(1,11)*percentile(2,plot_vj),'g')
% plot(-5:5,ones(1,11)*percentile(3,plot_vj),'k')
% plot(-5:5,ones(1,11)*percentile(4,plot_vj),'c')
% plot(-5:5,ones(1,11)*percentile(5,plot_vj),'m')
xlim([xlims(1) xlims(2)])
% hChildren = get(h, 'Children');
% set(hChildren, 'Markersize', msize)


