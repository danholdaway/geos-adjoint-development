close all
clear
clc

cd /discover/nobackup/drholdaw/ExperimentData/Journal_Articles/moist_dev_mwr/

% A = ncread('jacobian_vs_heating_container.nc4','A');
% B = ncread('jacobian_vs_heating_container.nc4','B');

A  =   ncread('jVr_6h_colkcbl_filteron.nc4','A');
B  =   ncread('jVr_6h_colkcbl_filteron.nc4','B');
ivec = ncread('jVr_6h_colkcbl_filteron.nc4','ivec');

J = ncread('JACOBIAN_DEEP_SERIES.nc4','J');
kcbl = ncread('JACOBIAN_DEEP_SERIES.nc4','KCBL');
top = ncread('JACOBIAN_DEEP_SERIES.nc4','TOP');
depth = ncread('JACOBIAN_DEEP_SERIES.nc4','DEPTH');

cd /home/drholdaw/Lin_Moist_Physics/Moist_Dev_MWR/

clc

ivec = ivec - 1;

Nk = 72;
Nk_vec = 1:72;


line_wid_cont = 0.8;
line_wid_det = 0.6;
fontsize = 12;

%%%%%%%%%%%%%%%%%%%
J_Auu = zeros(Nk,Nk);
J_Aut = zeros(Nk,Nk);
J_Auq = zeros(Nk,Nk);

J_Avv = zeros(Nk,Nk);
J_Avt = zeros(Nk,Nk);
J_Avq = zeros(Nk,Nk);

J_Ht  = zeros(Nk,Nk);
J_Hq  = zeros(Nk,Nk);

J_Mt  = zeros(Nk,Nk);
J_Mq  = zeros(Nk,Nk);

J_Auu(:,:) = J(1,0*Nk+1:1*Nk,0*Nk+1:1*Nk);
J_Aut(:,:) = J(1,0*Nk+1:1*Nk,2*Nk+1:3*Nk);
J_Auq(:,:) = J(1,0*Nk+1:1*Nk,3*Nk+1:4*Nk);

J_Avv(:,:) = J(1,1*Nk+1:2*Nk,1*Nk+1:2*Nk);
J_Avt(:,:) = J(1,1*Nk+1:2*Nk,2*Nk+1:3*Nk);
J_Avq(:,:) = J(1,1*Nk+1:2*Nk,3*Nk+1:4*Nk);

J_Ht(:,:)  = J(1,2*Nk+1:3*Nk,2*Nk+1:3*Nk);
J_Hq(:,:)  = J(1,2*Nk+1:3*Nk,3*Nk+1:4*Nk);

J_Mt(:,:)  = J(1,3*Nk+1:4*Nk,2*Nk+1:3*Nk);
J_Mq(:,:)  = J(1,3*Nk+1:4*Nk,3*Nk+1:4*Nk);
%%%%%%%%%%%%%%%%%%%
topm = 2;


figure('Position',[247 376 1132 500])

subplot(2,4,1)
[C,h] = contour(Nk_vec(top-topm:end),Nk_vec(top-topm:end),J_Ht(top-topm:end,top-topm:end),'LineWidth',line_wid_cont);
climmax = max(abs(caxis)); caxis([-climmax climmax])
set(gca,'YDir','reverse')
title('(a) \partial H /\partial \theta (s^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Output Index (H)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Perturbation Index (\theta)','FontSize',fontsize,'FontName','TimesNewRoman')
colormap jet
pos1 = get(gca,'Position');
colorbar
set(gca,'Position',[pos1(1)-0.08 pos1(2) pos1(3) pos1(4)])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
box on
axis square

levlist = get(h,'LevelList');
indzero = find(levlist==0);
levlist(indzero) = [];
set(h,'LevelList',levlist)

% line([70.5 70.5],[1 72],'LineWidth',line_wid_det','Color',[0 0 0])
% line([69.5 69.5],[1 72],'LineWidth',line_wid_det','Color',[0 0 0])
line([45.5 45.5],[1 72],'LineWidth',line_wid_det','Color',[0 0 0])
line([44.5 44.5],[1 72],'LineWidth',line_wid_det','Color',[0 0 0])

get(gca,'Position')

subplot(2,4,2)
[C,h] = contour(Nk_vec(top-topm:end),Nk_vec(top-topm:end),J_Hq(top-topm:end,top-topm:end),'LineWidth',line_wid_cont);
climmax = max(abs(caxis)); caxis([-climmax climmax])
set(gca,'YDir','reverse')
title('(b) \partial H /\partial q (Ks^{-1}kg^{-1}kg)','FontSize',fontsize,'FontName','TimesNewRoman')
%ylabel('Output Index (H)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Perturbation Index (q)','FontSize',fontsize,'FontName','TimesNewRoman')
colormap jet
pos2 = get(gca,'Position');
colorbar
set(gca,'Position',[pos2(1)-0.04 pos2(2) pos2(3) pos2(4)])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
box on
axis square

    levlist = get(h,'LevelList');
    indzero = find(levlist==0);
    levlist(indzero) = [];
    set(h,'LevelList',levlist)
    
line([71.5 71.5],[1 72],'LineWidth',line_wid_det','Color',[0 0 0])
line([70.5 70.5],[1 72],'LineWidth',line_wid_det','Color',[0 0 0])


subplot(2,4,5)
[C,h] = contour(Nk_vec(top-topm:end),Nk_vec(top-topm:end),J_Mt(top-topm:end,top-topm:end),'LineWidth',line_wid_cont);
climmax = max(abs(caxis)); caxis([-climmax climmax])
set(gca,'YDir','reverse')
title('(c) \partial Q /\partial \theta (kgkg^{-1}s^{-1}K^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Output Index (Q)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Perturbation Index (\theta)','FontSize',fontsize,'FontName','TimesNewRoman')
colormap jet
pos3 = get(gca,'Position');
colorbar
set(gca,'Position',[pos3(1)-0.08 pos3(2) pos3(3) pos3(4)])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
box on
axis square

    levlist = get(h,'LevelList');
    indzero = find(levlist==0);
    levlist(indzero) = [];
    set(h,'LevelList',levlist)

% line([70.5 70.5],[1 72],'LineWidth',line_wid_det','Color',[0 0 0])
% line([69.5 69.5],[1 72],'LineWidth',line_wid_det','Color',[0 0 0])
line([45.5 45.5],[1 72],'LineWidth',line_wid_det','Color',[0 0 0])
line([44.5 44.5],[1 72],'LineWidth',line_wid_det','Color',[0 0 0])


subplot(2,4,6)
[C,h] = contour(Nk_vec(top-topm:end),Nk_vec(top-topm:end),J_Mq(top-topm:end,top-topm:end),'LineWidth',line_wid_cont);
climmax = max(abs(caxis)); caxis([-climmax climmax])
set(gca,'YDir','reverse')
title('(d) \partial Q /\partial q  (s^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
%ylabel('Output Index (Q)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Perturbation Index (q)','FontSize',fontsize,'FontName','TimesNewRoman')
colormap jet
pos4 = get(gca,'Position');
colorbar
set(gca,'Position',[pos4(1)-0.04 pos4(2) pos4(3) pos4(4)])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
box on
axis square

    levlist = get(h,'LevelList');
    indzero = find(levlist==0);
    levlist(indzero) = [];
    set(h,'LevelList',levlist)

line([71.5 71.5],[1 72],'LineWidth',line_wid_det','Color',[0 0 0])
line([70.5 70.5],[1 72],'LineWidth',line_wid_det','Color',[0 0 0])



% figure

i = 13;

A_tmp = zeros(ivec(i),4);
A_tmp(:,:) = A(i,1:ivec(i),:);
B_tmp = zeros(ivec(i),4);
B_tmp(:,:) = B(i,1:ivec(i),:);

[J1_sort,ind_J1] = sort(A_tmp(:,3),'ascend');
[J2_sort,ind_J2] = sort(A_tmp(:,4),'ascend');
[J3_sort,ind_J3] = sort(B_tmp(:,3),'ascend');
[J4_sort,ind_J4] = sort(B_tmp(:,4),'ascend');

num_points = length(A_tmp);
points_90 = round(90*num_points/100);
points_95 = round(95*num_points/100);
points_98 = round(98*num_points/100);

% fprintf('J1_95 = %g \n', J1_sort(points_95))
% fprintf('J2_95 = %g \n', J2_sort(points_95))
% fprintf('J3_95 = %g \n', J3_sort(points_95))
% fprintf('J4_95 = %g \n\n\n', J4_sort(points_95))


subplot(2,4,3)
scatter(A(i,ind_J1(1:points_95),1),J1_sort(1:points_95),'b.')
hold on
scatter(A(i,ind_J1(points_90+1:points_95),1),J1_sort(points_90+1:points_95),'g.')
scatter(A(i,ind_J1(points_95+1:points_98),1),J1_sort(points_95+1:points_98),'CData',[1.0 0.50 0],'Marker','.')
scatter(A(i,ind_J1(points_98+1:end),1),J1_sort(points_98+1:end),'r.')
set(gca,'YScale','log')
set(gca,'XScale','log')
ylim([10e-8 10e-3])
xlim([10e-8 10e-1])
ylabel('max(|\partial H/\partial \theta|)_{k}','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('max(|H\prime|)','FontSize',fontsize,'FontName','TimesNewRoman')
title('(e)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
box on
axis square
set(gca,'YTick',[10e-8 10e-7 10e-6 10e-5 10e-4 10e-3])
set(gca,'XTick',[10e-8  10e-6  10e-4  10e-2 ])

pos5 = get(gca,'Position');
set(gca,'Position',[pos5(1)+0.04 pos5(2) pos4(3) pos4(4)])

subplot(2,4,7)
scatter(A(i,ind_J1(1:points_95),2),J2_sort(1:points_95),'b.')
hold on
scatter(A(i,ind_J1(points_90+1:points_95),2),J2_sort(points_90+1:points_95),'g.')
scatter(A(i,ind_J1(points_95+1:points_98),2),J2_sort(points_95+1:points_98),'CData',[1.0 0.50 0],'Marker','.')
scatter(A(i,ind_J1(points_98+1:end),2),J2_sort(points_98+1:end),'r.')
set(gca,'YScale','log')
set(gca,'XScale','log')
xlim([10e-11 10e-4])
ylim([10e-11 10e-6])
ylabel('max(|\partial Q/\partial \theta|)_{k}','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('max(|Q\prime|)','FontSize',fontsize,'FontName','TimesNewRoman')
title('(g)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
box on
axis square
set(gca,'YTick',[10e-11 10e-10 10e-9 10e-8 10e-7 10e-6])
set(gca,'XTick',[10e-11  10e-9  10e-7  10e-5 ])

pos6 = get(gca,'Position');
set(gca,'Position',[pos6(1)+0.04 pos6(2) pos4(3) pos4(4)])

subplot(2,4,4)
scatter(A(i,ind_J1(1:points_95),1),J3_sort(1:points_95),'b.')
hold on
scatter(A(i,ind_J1(points_90+1:points_95),1),J3_sort(points_90+1:points_95),'g.')
scatter(A(i,ind_J1(points_95+1:points_98),1),J3_sort(points_95+1:points_98),'CData',[1.0 0.50 0],'Marker','.')
scatter(A(i,ind_J1(points_98+1:end),1),J3_sort(points_98+1:end),'r.')
set(gca,'YScale','log')
set(gca,'XScale','log')
xlim([10e-8 10e-1])
ylim([10e-5 10e-0])
ylabel('max(|\partial H/\partial q|)_{k}','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('max(|H\prime|)','FontSize',fontsize,'FontName','TimesNewRoman')
title('(f)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
box on
axis square
set(gca,'YTick',[10e-5 10e-4 10e-3 10e-2 10e-1 10e-0])
set(gca,'XTick',[10e-8  10e-6  10e-4  10e-2 ])
legend('<90%','>90%','>95%','>98%','Location','SouthEast')

pos7 = get(gca,'Position');
set(gca,'Position',[pos7(1)+0.06 pos7(2) pos4(3) pos4(4)])

subplot(2,4,8)
scatter(A(i,ind_J1(1:points_95),2),J4_sort(1:points_95),'b.')
hold on
scatter(A(i,ind_J1(points_90+1:points_95),2),J4_sort(points_90+1:points_95),'g.')
scatter(A(i,ind_J1(points_95+1:points_98),2),J4_sort(points_95+1:points_98),'CData',[1.0 0.50 0],'Marker','.')
scatter(A(i,ind_J1(points_98+1:end),2),J4_sort(points_98+1:end),'r.')
set(gca,'YScale','log')
set(gca,'XScale','log')
xlim([10e-11 10e-4])
ylim([10e-8  10e-2])
ylabel('max(|\partial Q/\partial q|)_{k}','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('max(|Q\prime|)','FontSize',fontsize,'FontName','TimesNewRoman')
title('(h)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
box on
axis square
set(gca,'YTick',[10e-8 10e-7 10e-6 10e-5 10e-4 10e-3 10e-2])
set(gca,'XTick',[10e-11  10e-9  10e-7  10e-5 ])

pos8 = get(gca,'Position');
set(gca,'Position',[pos8(1)+0.06 pos8(2) pos4(3) pos4(4)])
set(gca,'XMinorTick','on')
    
saveas(gcf,'filter_scatter.eps', 'psc2')



