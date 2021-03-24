close all
clear
clc

heating_vs_heatingpert_container_6hr

scrsz = get(0,'ScreenSize');

fontsize = 22;

%%%%%%%%%%%%
figure('visible','on','Position',[0 scrsz(4) 0.5*scrsz(3) scrsz(4)]) 
subplot(2,4,1)
scatter(R6(:,1), R6(:,3))
set(gca,'YScale','log')
set(gca,'XScale','log')
title('\partial H/\partial \theta\prime','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Jacobian vs H\prime','FontSize',fontsize,'FontName','TimesNewRoman')
box on

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(2,4,2)
scatter(R6(:,1), R6(:,4))
set(gca,'YScale','log')
set(gca,'XScale','log')
title('\partial M/\partial \theta\prime','FontSize',fontsize,'FontName','TimesNewRoman')
box on

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(2,4,3)
scatter(R6(:,1), P6(:,3))
set(gca,'YScale','log')
set(gca,'XScale','log')
title('J1')
title('\partial H/\partial q\prime','FontSize',fontsize,'FontName','TimesNewRoman')
box on

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(2,4,4)
scatter(R6(:,1), P6(:,4))
set(gca,'YScale','log')
set(gca,'XScale','log')
title('J1')
title('\partial M/\partial q\prime','FontSize',fontsize,'FontName','TimesNewRoman')
box on

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(2,4,5)
scatter(R6(:,2), R6(:,3))
set(gca,'YScale','log')
set(gca,'XScale','log')
xlim([10e-15 10e-4])
ylabel('Jacobian vs M\prime','FontSize',fontsize,'FontName','TimesNewRoman')
box on

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(2,4,6)
scatter(R6(:,2), R6(:,4))
set(gca,'YScale','log')
set(gca,'XScale','log')
xlim([10e-15 10e-4])
box on

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(2,4,7)
scatter(R6(:,2), P6(:,3))
set(gca,'YScale','log')
set(gca,'XScale','log')
xlim([10e-15 10e-4])
box on

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(2,4,8)
scatter(R6(:,2), P6(:,4))
set(gca,'YScale','log')
set(gca,'XScale','log')
xlim([10e-15 10e-4])
box on

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

close

k = 10000;

J1_max = zeros(length(R6));
J2_max = zeros(length(R6));
J3_max = zeros(length(R6));
J4_max = zeros(length(R6));

J1_max(1:length(R6)) = R6(:,3);
J2_max(1:length(R6)) = R6(:,4);
J3_max(1:length(R6)) = P6(:,3);
J4_max(1:length(R6)) = P6(:,4);

R_keep = zeros(length(R6),4);
R_disc = zeros(length(R6),4);

count_J1 = 0;
count_J2 = 0;
count_J3 = 0;
count_J4 = 0;

count_J4_keep = 0;

for i = 1:length(R6)

    if J4_max(i) > 10^(-3.5)

        count_J4 = count_J4 + 1;
        R_disc(count_J4,:) = R6(count_J4,:);
        
    elseif J4_max(i) >= 0

        count_J4_keep = count_J4_keep + 1;
        R_keep(count_J4_keep,:) = R6(count_J4_keep,:);

    end

end

%%%%%%%%%%%%
figure('visible','on','Position',[0 scrsz(4) 0.5*scrsz(3) scrsz(4)]) 
subplot(2,4,1)
scatter(R_keep(:,1), R_keep(:,3))
hold on
scatter(R_disc(:,1), R_disc(:,3),'r')
set(gca,'YScale','log')
set(gca,'XScale','log')
title('\partial H/\partial \theta\prime','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Jacobian vs H\prime','FontSize',fontsize,'FontName','TimesNewRoman')
box on

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(2,4,2)
scatter(R_keep(:,1), R_keep(:,4))
hold on
scatter(R_disc(:,1), R_disc(:,4),'r')
set(gca,'YScale','log')
set(gca,'XScale','log')
title('\partial M/\partial \theta\prime','FontSize',fontsize,'FontName','TimesNewRoman')
box on

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

% subplot(2,4,3)
% scatter(R6(:,1), P6(:,3))
% set(gca,'YScale','log')
% set(gca,'XScale','log')
% title('J1')
% title('\partial H/\partial q\prime','FontSize',fontsize,'FontName','TimesNewRoman')
% box on
% 
% set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% 
% subplot(2,4,4)
% scatter(R6(:,1), P6(:,4))
% set(gca,'YScale','log')
% set(gca,'XScale','log')
% title('J1')
% title('\partial M/\partial q\prime','FontSize',fontsize,'FontName','TimesNewRoman')
% box on

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(2,4,5)
scatter(R_keep(:,2), R_keep(:,3))
hold on
scatter(R_disc(:,2), R_disc(:,3),'r')
set(gca,'YScale','log')
set(gca,'XScale','log')
xlim([10e-15 10e-4])
ylabel('Jacobian vs M\prime','FontSize',fontsize,'FontName','TimesNewRoman')
box on

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(2,4,6)
scatter(R_keep(:,2), R_keep(:,4))
hold on
scatter(R_disc(:,2), R_disc(:,4),'r')
set(gca,'YScale','log')
set(gca,'XScale','log')
xlim([10e-15 10e-4])
box on

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

% subplot(2,4,7)
% scatter(R6(:,2), P6(:,3))
% set(gca,'YScale','log')
% set(gca,'XScale','log')
% xlim([10e-15 10e-4])
% box on
% 
% set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% 
% subplot(2,4,8)
% scatter(R6(:,2), P6(:,4))
% set(gca,'YScale','log')
% set(gca,'XScale','log')
% xlim([10e-15 10e-4])
% box on

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
