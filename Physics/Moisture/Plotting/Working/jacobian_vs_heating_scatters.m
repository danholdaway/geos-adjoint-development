close all
clear
clc

addpath(cd)

cd /discover/nobackup/drholdaw/tmp.9748/sens.20130112.000000/
% asd
H_M_J_tkcbl_qkcblm1_66steps


cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/


ivec = i_vec - 1;

scrsz = get(0,'ScreenSize');

imax = 60;

J1_CUT = zeros(4,imax);
J2_CUT = zeros(4,imax);
J3_CUT = zeros(4,imax);
J4_CUT = zeros(4,imax);

cut1 = 95;
cut2 = 96;
cut3 = 97;
cut4 = 98;

cut_vec = [cut1 cut2 cut3 cut4]';

for i = 1:imax

    A_tmp = zeros(ivec(i),4);
    A_tmp(:,:) = A(i,1:ivec(i),:);
    B_tmp = zeros(ivec(i),4);
    B_tmp(:,:) = B(i,1:ivec(i),:);
        
    [J1_sort,ind_J1] = sort(A_tmp(:,3),'ascend');
    H_byJ1 = A_tmp(ind_J1,1);
    M_byJ1 = A_tmp(ind_J1,2);
    [J2_sort,ind_J2] = sort(A_tmp(:,4),'ascend');
    H_byJ2 = A_tmp(ind_J2,1);
    M_byJ2 = A_tmp(ind_J2,2);
    [J3_sort,ind_J3] = sort(B_tmp(:,3),'ascend');
    H_byJ3 = A_tmp(ind_J3,1);
    M_byJ3 = A_tmp(ind_J3,2);
    [J4_sort,ind_J4] = sort(B_tmp(:,4),'ascend');
    H_byJ4 = A_tmp(ind_J4,1);
    M_byJ4 = A_tmp(ind_J4,2);
    
    num_points = length(A_tmp);
    points_1 = round(cut1*num_points/100);
    points_2 = round(cut2*num_points/100);
    points_3 = round(cut3*num_points/100);
    points_4 = round(cut4*num_points/100);
        
    %FIND THE THRESHOLDS BY Js FOR EACH PERCENTAGE CUT OFF
    J1_CUT(1,i) = J1_sort(points_1);
    J2_CUT(1,i) = J2_sort(points_1);
    J3_CUT(1,i) = J3_sort(points_1);
    J4_CUT(1,i) = J4_sort(points_1);
    
    J1_CUT(2,i) = J1_sort(points_2);
    J2_CUT(2,i) = J2_sort(points_2);
    J3_CUT(2,i) = J3_sort(points_2);
    J4_CUT(2,i) = J4_sort(points_2);
    
    J1_CUT(3,i) = J1_sort(points_3);
    J2_CUT(3,i) = J2_sort(points_3);
    J3_CUT(3,i) = J3_sort(points_3);
    J4_CUT(3,i) = J4_sort(points_3);
    
    J1_CUT(4,i) = J1_sort(points_4);
    J2_CUT(4,i) = J2_sort(points_4);
    J3_CUT(4,i) = J3_sort(points_4);
    J4_CUT(4,i) = J4_sort(points_4);
 
    
%     figure('visible','on','Position',[0 scrsz(4) 0.5*scrsz(3) scrsz(4)]) 
%     
%     subplot(2,4,1)
%     scatter(H_byJ1(1:points_1),J1_sort(1:points_1),'m')
%     hold on
%     scatter(H_byJ1(points_1+1:points_2),J1_sort(points_1+1:points_2),'b')
%     scatter(H_byJ1(points_2+1:points_3),J1_sort(points_2+1:points_3),'r')
%     scatter(H_byJ1(points_3+1:points_4),J1_sort(points_3+1:points_4),'g')
%     scatter(H_byJ1(points_4+1:end),      J1_sort(points_4+1:end),      'c')
%     
%     ylim([min(J1_sort) 2*max(J1_sort)])
%     xlim([min(H_byJ1) 2*max(H_byJ1)])
%     set(gca,'YScale','log')
%     set(gca,'XScale','log')
%     
%     subplot(2,4,5)
%     scatter(M_byJ1(1:points_1),J1_sort(1:points_1),'m')
%     hold on
%     scatter(M_byJ1(points_1+1:points_2),J1_sort(points_1+1:points_2),'b')
%     scatter(M_byJ1(points_2+1:points_3),J1_sort(points_2+1:points_3),'r')
%     scatter(M_byJ1(points_3+1:points_4),J1_sort(points_3+1:points_4),'g')
%     scatter(M_byJ1(points_4+1:end),      J1_sort(points_4+1:end),      'c')
%     
%     
%     ylim([min(J1_sort) 2*max(J1_sort)])
%     xlim([min(M_byJ1) 2*max(M_byJ1)])
%     set(gca,'YScale','log')
%     set(gca,'XScale','log')
%     
%     
%     
%     subplot(2,4,2)
%     scatter(H_byJ2(1:points_1),J2_sort(1:points_1),'m')
%     hold on
%     scatter(H_byJ2(points_1+1:points_2),J2_sort(points_1+1:points_2),'b')
%     scatter(H_byJ2(points_2+1:points_3),J2_sort(points_2+1:points_3),'r')
%     scatter(H_byJ2(points_3+1:points_4),J2_sort(points_3+1:points_4),'g')
%     scatter(H_byJ2(points_4+1:end),      J2_sort(points_4+1:end),      'c')
%     
%     ylim([min(J2_sort) 2*max(J2_sort)])
%     xlim([min(H_byJ2) 2*max(H_byJ2)])
%     set(gca,'YScale','log')
%     set(gca,'XScale','log')
%     
%     subplot(2,4,6)
%     scatter(M_byJ2(1:points_1),J2_sort(1:points_1),'m')
%     hold on
%     scatter(M_byJ2(points_1+1:points_2),J2_sort(points_1+1:points_2),'b')
%     scatter(M_byJ2(points_2+1:points_3),J2_sort(points_2+1:points_3),'r')
%     scatter(M_byJ2(points_3+1:points_4),J2_sort(points_3+1:points_4),'g')
%     scatter(M_byJ2(points_4+1:end),      J2_sort(points_4+1:end),      'c')
%     
%     
%     ylim([min(J2_sort) 2*max(J2_sort)])
%     xlim([min(M_byJ2) 2*max(M_byJ2)])
%     set(gca,'YScale','log')
%     set(gca,'XScale','log')
%     
%     
%     
%     subplot(2,4,3)
%     scatter(H_byJ3(1:points_1),J3_sort(1:points_1),'m')
%     hold on
%     scatter(H_byJ3(points_1+1:points_2),J3_sort(points_1+1:points_2),'b')
%     scatter(H_byJ3(points_2+1:points_3),J3_sort(points_2+1:points_3),'r')
%     scatter(H_byJ3(points_3+1:points_4),J3_sort(points_3+1:points_4),'g')
%     scatter(H_byJ3(points_4+1:end),      J3_sort(points_4+1:end),      'c')
%     
%     ylim([min(J3_sort) 2*max(J3_sort)])
%     xlim([min(H_byJ3) 2*max(H_byJ3)])
%     set(gca,'YScale','log')
%     set(gca,'XScale','log')
%     
%     subplot(2,4,7)
%     scatter(M_byJ3(1:points_1),J3_sort(1:points_1),'m')
%     hold on
%     scatter(M_byJ3(points_1+1:points_2),J3_sort(points_1+1:points_2),'b')
%     scatter(M_byJ3(points_2+1:points_3),J3_sort(points_2+1:points_3),'r')
%     scatter(M_byJ3(points_3+1:points_4),J3_sort(points_3+1:points_4),'g')
%     scatter(M_byJ3(points_4+1:end),      J3_sort(points_4+1:end),      'c')
%     
%     
%     ylim([min(J3_sort) 2*max(J3_sort)])
%     xlim([min(M_byJ3) 2*max(M_byJ3)])
%     set(gca,'YScale','log')
%     set(gca,'XScale','log')
%     
%     
%     
%     subplot(2,4,4)
%     scatter(H_byJ4(1:points_1),J4_sort(1:points_1),'m')
%     hold on
%     scatter(H_byJ4(points_1+1:points_2),J4_sort(points_1+1:points_2),'b')
%     scatter(H_byJ4(points_2+1:points_3),J4_sort(points_2+1:points_3),'r')
%     scatter(H_byJ4(points_3+1:points_4),J4_sort(points_3+1:points_4),'g')
%     scatter(H_byJ4(points_4+1:end),      J4_sort(points_4+1:end),      'c')
%     
%     ylim([min(J4_sort) 2*max(J4_sort)])
%     xlim([min(H_byJ4) 2*max(H_byJ4)])
%     set(gca,'YScale','log')
%     set(gca,'XScale','log')
%     
%     subplot(2,4,8)
%     scatter(M_byJ4(1:points_1),J4_sort(1:points_1),'m')
%     hold on
%     scatter(M_byJ4(points_1+1:points_2),J4_sort(points_1+1:points_2),'b')
%     scatter(M_byJ4(points_2+1:points_3),J4_sort(points_2+1:points_3),'r')
%     scatter(M_byJ4(points_3+1:points_4),J4_sort(points_3+1:points_4),'g')
%     scatter(M_byJ4(points_4+1:end),      J4_sort(points_4+1:end),      'c')
%     
%     
%     ylim([min(J4_sort) 2*max(J4_sort)])
%     xlim([min(M_byJ4) 2*max(M_byJ4)])
%     set(gca,'YScale','log')
%     set(gca,'XScale','log')
    
end
  
for j = 1:4
    
    J1_CUT(j,imax+1) = mean(J1_CUT(j,1:imax));
    J2_CUT(j,imax+1) = mean(J2_CUT(j,1:imax));
    J3_CUT(j,imax+1) = mean(J3_CUT(j,1:imax));
    J4_CUT(j,imax+1) = mean(J4_CUT(j,1:imax));
    
    J1_CUT(j,imax+2) = max(J1_CUT(j,1:imax));
    J2_CUT(j,imax+2) = max(J2_CUT(j,1:imax));
    J3_CUT(j,imax+2) = max(J3_CUT(j,1:imax));
    J4_CUT(j,imax+2) = max(J4_CUT(j,1:imax));
    
end

J1_CUT = [cut_vec J1_CUT];
J2_CUT = [cut_vec J2_CUT];
J3_CUT = [cut_vec J3_CUT];
J4_CUT = [cut_vec J4_CUT];


fprintf('MEANS \n')
fprintf('   !Filter %g percent \n', 100-cut1)
fprintf('   J1_cut = %g \n', J1_CUT(1,imax+2))
fprintf('   J2_cut = %g \n', J3_CUT(1,imax+2))
fprintf('   J3_cut = %g \n', J2_CUT(1,imax+2))
fprintf('   J4_cut = %g \n \n', J4_CUT(1,imax+2))

fprintf('   !Filter %g percent \n', 100-cut2)
fprintf('   J1_cut = %g \n', J1_CUT(2,imax+2))
fprintf('   J2_cut = %g \n', J3_CUT(2,imax+2))
fprintf('   J3_cut = %g \n', J2_CUT(2,imax+2))
fprintf('   J4_cut = %g \n \n', J4_CUT(2,imax+2))

fprintf('   !Filter %g percent \n', 100-cut3)
fprintf('   J1_cut = %g \n', J1_CUT(3,imax+2))
fprintf('   J2_cut = %g \n', J3_CUT(3,imax+2))
fprintf('   J3_cut = %g \n', J2_CUT(3,imax+2))
fprintf('   J4_cut = %g \n \n', J4_CUT(3,imax+2))

fprintf('   !Filter %g percent \n', 100-cut4)
fprintf('   J1_cut = %g \n', J1_CUT(4,imax+2))
fprintf('   J2_cut = %g \n', J3_CUT(4,imax+2))
fprintf('   J3_cut = %g \n', J2_CUT(4,imax+2))
fprintf('   J4_cut = %g \n \n', J4_CUT(4,imax+2))


fprintf('MAXIMUMS \n')
fprintf('   !Filter %g percent \n', 100-cut1)
fprintf('   !J1_cut = %g \n', J1_CUT(1,imax+3))
fprintf('   !J2_cut = %g \n', J3_CUT(1,imax+3))
fprintf('   !J3_cut = %g \n', J2_CUT(1,imax+3))
fprintf('   !J4_cut = %g \n \n', J4_CUT(1,imax+3))

fprintf('   !Filter %g percent \n', 100-cut2)
fprintf('   !J1_cut = %g \n', J1_CUT(2,imax+3))
fprintf('   !J2_cut = %g \n', J3_CUT(2,imax+3))
fprintf('   !J3_cut = %g \n', J2_CUT(2,imax+3))
fprintf('   !J4_cut = %g \n \n', J4_CUT(2,imax+3))

fprintf('   !Filter %g percent \n', 100-cut3)
fprintf('   !J1_cut = %g \n', J1_CUT(3,imax+3))
fprintf('   !J2_cut = %g \n', J3_CUT(3,imax+3))
fprintf('   !J3_cut = %g \n', J2_CUT(3,imax+3))
fprintf('   !J4_cut = %g \n \n', J4_CUT(3,imax+3))

fprintf('   !Filter %g percent \n', 100-cut4)
fprintf('   J1_cut = %g \n', J1_CUT(4,imax+3))
fprintf('   J2_cut = %g \n', J3_CUT(4,imax+3))
fprintf('   J3_cut = %g \n', J2_CUT(4,imax+3))
fprintf('   J4_cut = %g \n \n', J4_CUT(4,imax+3))


