close all
clear
clc

cd /discover/nobackup/drholdaw/ExperimentData/jacobianVSheatingrate/

% A = ncread('jacobian_vs_heating_container.nc4','A');
% B = ncread('jacobian_vs_heating_container.nc4','B');

A  =   ncread('jVr_6h_colkcbl_filteron.nc4','A');
B  =   ncread('jVr_6h_colkcbl_filteron.nc4','B');
ivec = ncread('jVr_6h_colkcbl_filteron.nc4','ivec');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

clc

ivec = ivec - 1;

scrsz = get(0,'ScreenSize');

line_wid = 1;
fontsize = 10;

for i = 13%:18

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
    
    fprintf('J1_95 = %g \n', J1_sort(points_95))
    fprintf('J2_95 = %g \n', J2_sort(points_95))
    fprintf('J3_95 = %g \n', J3_sort(points_95))
    fprintf('J4_95 = %g \n\n\n', J4_sort(points_95))
    
    
%     figure('visible','on','Position',[0 scrsz(4) 0.5*scrsz(3) scrsz(4)]) 
%     subplot(2,4,1)
%     scatter(A(i,ind_J1(1:points_90),1),J1_sort(1:points_90),'m')
%     hold on
%     scatter(A(i,ind_J1(points_90+1:points_95),1),J1_sort(points_90+1:points_95),'b')
%     scatter(A(i,ind_J1(points_95+1:points_98),1),J1_sort(points_95+1:points_98),'r')
%     scatter(A(i,ind_J1(points_98+1:end),1),J1_sort(points_98+1:end),'g')
%     set(gca,'YScale','log')
%     set(gca,'XScale','log')
%     
%     subplot(2,4,2)
%     scatter(A(i,ind_J1(1:points_95),1),J2_sort(1:points_95),'m')
%     hold on
%     scatter(A(i,ind_J1(points_90+1:points_95),1),J2_sort(points_90+1:points_95),'b')
%     scatter(A(i,ind_J1(points_95+1:points_98),1),J2_sort(points_95+1:points_98),'r')
%     scatter(A(i,ind_J1(points_98+1:end),1),J2_sort(points_98+1:end),'g')
%     set(gca,'YScale','log')
%     set(gca,'XScale','log')
% 
%     subplot(2,4,3)
%     scatter(A(i,ind_J1(1:points_95),1),J3_sort(1:points_95),'m')
%     hold on
%     scatter(A(i,ind_J1(points_90+1:points_95),1),J3_sort(points_90+1:points_95),'b')
%     scatter(A(i,ind_J1(points_95+1:points_98),1),J3_sort(points_95+1:points_98),'r')
%     scatter(A(i,ind_J1(points_98+1:end),1),J3_sort(points_98+1:end),'g')
%     set(gca,'YScale','log')
%     set(gca,'XScale','log')
% 
%     subplot(2,4,4)
%     scatter(A(i,ind_J1(1:points_95),1),J4_sort(1:points_95),'m')
%     hold on
%     scatter(A(i,ind_J1(points_90+1:points_95),1),J4_sort(points_90+1:points_95),'b')
%     scatter(A(i,ind_J1(points_95+1:points_98),1),J4_sort(points_95+1:points_98),'r')
%     scatter(A(i,ind_J1(points_98+1:end),1),J4_sort(points_98+1:end),'g')
%     set(gca,'YScale','log')
%     set(gca,'XScale','log')

    subplot(2,2,1)
    scatter(A(i,ind_J1(1:points_95),1),J1_sort(1:points_95),'m')
    hold on
    scatter(A(i,ind_J1(points_90+1:points_95),1),J1_sort(points_90+1:points_95),'b')
    scatter(A(i,ind_J1(points_95+1:points_98),1),J1_sort(points_95+1:points_98),'r')
    scatter(A(i,ind_J1(points_98+1:end),1),J1_sort(points_98+1:end),'g')
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    ylim([10e-8 10e-3])
    xlim([10e-8 10e-1])
    ylabel('Max(|J1|)','FontSize',fontsize,'FontName','TimesNewRoman')
    xlabel('Max(|H^\prime|)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
box on
   
    subplot(2,2,2)
    scatter(A(i,ind_J1(1:points_95),1),J2_sort(1:points_95),'m')
    hold on
    scatter(A(i,ind_J1(points_90+1:points_95),1),J2_sort(points_90+1:points_95),'b')
    scatter(A(i,ind_J1(points_95+1:points_98),1),J2_sort(points_95+1:points_98),'r')
    scatter(A(i,ind_J1(points_98+1:end),1),J2_sort(points_98+1:end),'g')
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    xlim([10e-8 10e-1])
    ylabel('Max(|J2|)','FontSize',fontsize,'FontName','TimesNewRoman')
    xlabel('Max(|H^\prime|)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
box on

    subplot(2,2,3)
    scatter(A(i,ind_J1(1:points_95),2),J3_sort(1:points_95),'m')
    hold on
    scatter(A(i,ind_J1(points_90+1:points_95),2),J3_sort(points_90+1:points_95),'b')
    scatter(A(i,ind_J1(points_95+1:points_98),2),J3_sort(points_95+1:points_98),'r')
    scatter(A(i,ind_J1(points_98+1:end),2),J3_sort(points_98+1:end),'g')
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    xlim([10e-11 10e-4])
    ylabel('Max(|J3|)','FontSize',fontsize,'FontName','TimesNewRoman')
    xlabel('Max(|Q^\prime|)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
box on

    subplot(2,2,4)
    scatter(A(i,ind_J1(1:points_95),2),J4_sort(1:points_95),'m')
    hold on
    scatter(A(i,ind_J1(points_90+1:points_95),2),J4_sort(points_90+1:points_95),'b')
    scatter(A(i,ind_J1(points_95+1:points_98),2),J4_sort(points_95+1:points_98),'r')
    scatter(A(i,ind_J1(points_98+1:end),2),J4_sort(points_98+1:end),'g')
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    xlim([10e-11 10e-4])
    ylim([10e-8  10e-2])
    ylabel('Max(|J4|)','FontSize',fontsize,'FontName','TimesNewRoman')
    xlabel('Max(|Q^\prime|)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
box on
    
%     subplot(2,4,9)
%     scatter(A(i,:,1),C(i,:,3))
%     set(gca,'YScale','log')
%     set(gca,'XScale','log')
% 
%     subplot(2,4,10)
%     scatter(A(i,:,1),C(i,:,4))
%     set(gca,'YScale','log')
%     set(gca,'XScale','log')
% 
%     subplot(2,4,11)
%     scatter(A(i,:,1),D(i,:,3))
%     set(gca,'XScale','log')
% 
%     subplot(2,4,12)
%     scatter(A(i,:,1),D(i,:,4))
%     set(gca,'YScale','log')
%     set(gca,'XScale','log')
% 
%     subplot(2,4,13)
%     scatter(A(i,:,2),C(i,:,3))
%     set(gca,'YScale','log')
%     set(gca,'XScale','log')
% 
%     subplot(2,4,14)
%     scatter(A(i,:,2),C(i,:,4))
%     set(gca,'YScale','log')
%     set(gca,'XScale','log')
% 
%     subplot(2,4,15)
%     scatter(A(i,:,2),D(i,:,3))
%     set(gca,'XScale','log')
% 
%     subplot(2,4,16)
%     scatter(A(i,:,2),D(i,:,4))
%     set(gca,'YScale','log')
%     set(gca,'XScale','log')
    
end
    