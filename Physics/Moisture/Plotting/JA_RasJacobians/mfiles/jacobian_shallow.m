close all
close all
clear
clc


J = ncread('JACOBIAN_SHALLOW.nc4','J');
kcbl = ncread('JACOBIAN_SHALLOW.nc4','KCBL');
top = ncread('JACOBIAN_SHALLOW.nc4','TOP');
depth = ncread('JACOBIAN_SHALLOW.nc4','DEPTH');


Nk = 72;

Nk_vec = 1:72;

vec1 = [top kcbl];
vec2 = ones(1,2);

% x = zeros(4*Nk,1);
% x(0*Nk+1:1*Nk) = U(150,100,:);
% x(1*Nk+1:2*Nk) = V(150,100,:);
% x(2*Nk+1:3*Nk) = TV(150,100,:);
% x(3*Nk+1:4*Nk) = QV(150,100,:);
% J = jtimesx(J(1:4*Nk,1:4*Nk),x);

%%%%%%%%%%%%%%%%%%%
J_Auu = J(0*Nk+1:1*Nk,0*Nk+1:1*Nk);
J_Aut = J(0*Nk+1:1*Nk,2*Nk+1:3*Nk);
J_Auq = J(0*Nk+1:1*Nk,3*Nk+1:4*Nk);

J_Avv = J(1*Nk+1:2*Nk,1*Nk+1:2*Nk);
J_Avt = J(1*Nk+1:2*Nk,2*Nk+1:3*Nk);
J_Avq = J(1*Nk+1:2*Nk,3*Nk+1:4*Nk);

J_Ht  = J(2*Nk+1:3*Nk,2*Nk+1:3*Nk);
J_Hq  = J(2*Nk+1:3*Nk,3*Nk+1:4*Nk);

J_Mt  = J(3*Nk+1:4*Nk,2*Nk+1:3*Nk);
J_Mq  = J(3*Nk+1:4*Nk,3*Nk+1:4*Nk);
%%%%%%%%%%%%%%%%%%%

percent_cut = 0.01;

%REMOVE SMALL ELEMENT OF J_Ht
J_Ht = remove0(J_Ht,percent_cut);
J_Hq = remove0(J_Hq,percent_cut);
J_Mt = remove0(J_Mt,percent_cut);
J_Mq = remove0(J_Mq,percent_cut);

J_Auu = remove0(J_Auu,percent_cut);
J_Aut = remove0(J_Aut,percent_cut);
J_Auq = remove0(J_Auq,percent_cut);

J_Avv = remove0(J_Avv,percent_cut);
J_Avt = remove0(J_Avt,percent_cut);
J_Avq = remove0(J_Avq,percent_cut);



line_wid = 1;
fontsize = 8;

scrsz = get(0,'ScreenSize');
figure('visible','on','Position',[scrsz(1) scrsz(2) (0.22)*scrsz(3) (0.4)*scrsz(4)])

subplot(2,2,1)
[C,h] = contour(Nk_vec(top-5:end),Nk_vec(top-5:end),J_Ht(top-5:end,top-5:end));
climmax = max(abs(caxis)); caxis([-climmax climmax])
set(h,'LineWidth',line_wid);
set(gca,'YDir','reverse')
title('\partial H /\partial \theta','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Output Index (H)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Perturbation Index (\theta)','FontSize',fontsize,'FontName','TimesNewRoman')
colormap jet
colorbar
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(2,2,2)
[C,h] = contour(Nk_vec(top-5:end),Nk_vec(top-5:end),J_Hq(top-5:end,top-5:end));
climmax = max(abs(caxis)); caxis([-climmax climmax])
set(h,'LineWidth',line_wid);
set(gca,'YDir','reverse')
title('\partial H /\partial q','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Output Index (H)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Perturbation Index (q)','FontSize',fontsize,'FontName','TimesNewRoman')
colormap jet
colorbar
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(2,2,3)
[C,h] = contour(Nk_vec(top-5:end),Nk_vec(top-5:end),J_Mt(top-5:end,top-5:end));
climmax = max(abs(caxis)); caxis([-climmax climmax])
set(h,'LineWidth',line_wid);
set(gca,'YDir','reverse')
title('\partial Q /\partial \theta','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Output Index (Q)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Perturbation Index (\theta)','FontSize',fontsize,'FontName','TimesNewRoman')
colormap jet
colorbar
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(2,2,4)
[C,h] = contour(Nk_vec(top-5:end),Nk_vec(top-5:end),J_Mq(top-5:end,top-5:end));
climmax = max(abs(caxis)); caxis([-climmax climmax])
set(h,'LineWidth',line_wid);
set(gca,'YDir','reverse')
title('\partial Q /\partial q','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Output Index (Q)','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Perturbation Index (q)','FontSize',fontsize,'FontName','TimesNewRoman')
colormap jet
colorbar
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

saveas(gcf,'jacobian_shallow_hq.eps', 'psc2')


% figure('visible','on','Position',[scrsz(1) scrsz(2) (0.3)*scrsz(3) (0.5)*scrsz(4)])
% 
% subplot(2,3,1)
% [C,h] = contour(Nk_vec(top-5:end),Nk_vec(top-5:end),J_Auu(top-5:end,top-5:end));
% climmax = max(abs(caxis)); caxis([-climmax climmax])
% set(h,'LineWidth',line_wid);
% set(gca,'YDir','reverse')
% title('\partial A_u /\partial u','FontSize',fontsize,'FontName','TimesNewRoman')
% ylabel('Output Index (A_u)','FontSize',fontsize,'FontName','TimesNewRoman')
% xlabel('Perturbation Index (u)','FontSize',fontsize,'FontName','TimesNewRoman')
% colormap jet
% colorbar
% set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% 
% subplot(2,3,2)
% [C,h] = contour(Nk_vec(top-5:end),Nk_vec(top-5:end),J_Aut(top-5:end,top-5:end));
% climmax = max(abs(caxis)); caxis([-climmax climmax])
% set(h,'LineWidth',line_wid);
% set(gca,'YDir','reverse')
% title('\partial A_u /\partial \theta','FontSize',fontsize,'FontName','TimesNewRoman')
% ylabel('Output Index (A_u)','FontSize',fontsize,'FontName','TimesNewRoman')
% xlabel('Perturbation Index (\theta)','FontSize',fontsize,'FontName','TimesNewRoman')
% colormap jet
% colorbar
% set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% 
% subplot(2,3,3)
% [C,h] = contour(Nk_vec(top-5:end),Nk_vec(top-5:end),J_Auq(top-5:end,top-5:end));
% climmax = max(abs(caxis)); caxis([-climmax climmax])
% set(h,'LineWidth',line_wid);
% set(gca,'YDir','reverse')
% title('\partial A_u /\partial q','FontSize',fontsize,'FontName','TimesNewRoman')
% ylabel('Output Index (A_u)','FontSize',fontsize,'FontName','TimesNewRoman')
% xlabel('Perturbation Index (q)','FontSize',fontsize,'FontName','TimesNewRoman')
% colormap jet
% colorbar
% set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% 
% subplot(2,3,4)
% [C,h] = contour(Nk_vec(top-5:end),Nk_vec(top-5:end),J_Avv(top-5:end,top-5:end));
% climmax = max(abs(caxis)); caxis([-climmax climmax])
% set(h,'LineWidth',line_wid);
% set(gca,'YDir','reverse')
% title('\partial A_v /\partial v','FontSize',fontsize,'FontName','TimesNewRoman')
% ylabel('Output Index (A_v)','FontSize',fontsize,'FontName','TimesNewRoman')
% xlabel('Perturbation Index (v)','FontSize',fontsize,'FontName','TimesNewRoman')
% colormap jet
% colorbar
% set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% 
% subplot(2,3,5)
% [C,h] = contour(Nk_vec(top-5:end),Nk_vec(top-5:end),J_Avt(top-5:end,top-5:end));
% climmax = max(abs(caxis)); caxis([-climmax climmax])
% set(h,'LineWidth',line_wid);
% set(gca,'YDir','reverse')
% title('\partial A_v /\partial \theta','FontSize',fontsize,'FontName','TimesNewRoman')
% ylabel('Output Index (A_v)','FontSize',fontsize,'FontName','TimesNewRoman')
% xlabel('Perturbation Index (\theta)','FontSize',fontsize,'FontName','TimesNewRoman')
% colormap jet
% colorbar
% set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% 
% subplot(2,3,6)
% [C,h] = contour(Nk_vec(top-5:end),Nk_vec(top-5:end),J_Avq(top-5:end,top-5:end));
% climmax = max(abs(caxis)); caxis([-climmax climmax])
% set(h,'LineWidth',line_wid);
% set(gca,'YDir','reverse')
% title('\partial A_v /\partial q','FontSize',fontsize,'FontName','TimesNewRoman')
% ylabel('Output Index (A_v)','FontSize',fontsize,'FontName','TimesNewRoman')
% xlabel('Perturbation Index (q)','FontSize',fontsize,'FontName','TimesNewRoman')
% colormap jet
% colorbar
% set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% 
% saveas(gcf,'jacobian_shallow_a.eps', 'psc2')





