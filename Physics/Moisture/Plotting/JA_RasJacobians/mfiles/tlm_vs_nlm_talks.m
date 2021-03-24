close all
clear
clc

cd /home/drholdaw/Lin_Moist_Physics/Jacobian_Exact/

J_tlm = ncread('JACOBIAN_tlm.nc4','J_TL');

k = 72;

J_Ht_tlm = J_tlm(2*k+1:3*k,2*k+1:3*k);
J_Hq_tlm = J_tlm(2*k+1:3*k,3*k+1:4*k);

J_Mt_tlm = J_tlm(3*k+1:4*k,2*k+1:3*k);
J_Mq_tlm = J_tlm(3*k+1:4*k,3*k+1:4*k);

cd /home/drholdaw/Lin_Moist_Physics/Jacobian_Pert/

J_nlm1 = ncread('JACOBIAN_nlm.nc4','J');

cd /home/drholdaw/Lin_Moist_Physics/RAS_Jacobian_Paper/mfiles/

Nk = 72;

J_nlm = zeros(360,360);
J_nlm(:,:) = J_nlm1(1,:,:);

J_Ht_nlm  = J_nlm(2*Nk+1:3*Nk,2*Nk+1:3*Nk);
J_Hq_nlm  = J_nlm(2*Nk+1:3*Nk,3*Nk+1:4*Nk);

J_Mt_nlm  = J_nlm(3*Nk+1:4*Nk,2*Nk+1:3*Nk);
J_Mq_nlm  = J_nlm(3*Nk+1:4*Nk,3*Nk+1:4*Nk);

line_wid = 1;
fontsize = 10;

top = 40;
topm = 7;

Nk_vec = 1:72;

figure
subplot(1,2,1)
contour(Nk_vec(top-topm:end),Nk_vec(top-topm:end),J_Ht_nlm(top-topm:end,top-topm:end))
colorbar
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Output Index (H)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Perturbation Index (\theta)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YDir','reverse')
title('\partial H /\partial \theta (s^{-1}) Nonlinear model','FontSize',fontsize,'FontName','TimesNewRoman')


subplot(1,2,2)
contour(Nk_vec(top-topm:end),Nk_vec(top-topm:end),J_Ht_tlm(top-topm:end,top-topm:end))
colorbar
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% ylabel('Output Index (H)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Perturbation Index (\theta)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'YDir','reverse')
title('\partial H /\partial \theta (s^{-1}) Tangent linear model','FontSize',fontsize,'FontName','TimesNewRoman')
