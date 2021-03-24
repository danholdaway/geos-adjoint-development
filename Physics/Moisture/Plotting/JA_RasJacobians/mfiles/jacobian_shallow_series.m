close all
close all
clear
clc


Jac = ncread('JACOBIAN_SHALLOW_SERIES.nc4','J');
kcbl = ncread('JACOBIAN_SHALLOW_SERIES.nc4','KCBL');
top = ncread('JACOBIAN_SHALLOW_SERIES.nc4','TOP');
depth = ncread('JACOBIAN_SHALLOW_SERIES.nc4','DEPTH');

Nk = 72;
Nk_vec = 1:72;

vec1 = [top kcbl];
vec2 = ones(1,2);

for k = 2:7
    n = 0;
    sum2 = 0;
    for i = 1:4*Nk
        for j = 1:4*Nk

            n = n+1;
            sum2 = sum2 + (Jac(1,i,j) - Jac(k,i,j))^2;

        end
    end
    RMS = sqrt(sum2/n)
end

% asd

scrsz = get(0,'ScreenSize');
figure('visible','on','Position',[scrsz(1) scrsz(2) (0.22)*scrsz(3) (0.4)*scrsz(4)])

imax = 4;

J_vec = [1 4 6 7];
Jac_new = Jac(J_vec,:,:);

J_max = max(max(max(Jac_new(:,2*Nk+1:3*Nk,2*Nk+1:3*Nk))));
J_min = min(min(min(Jac_new(:,2*Nk+1:3*Nk,2*Nk+1:3*Nk))));
climmax1 = max([J_max abs(J_min)]);

J_max = max(max(max(Jac_new(:,2*Nk+1:3*Nk,3*Nk+1:4*Nk))));
J_min = min(min(min(Jac_new(:,2*Nk+1:3*Nk,3*Nk+1:4*Nk))));
climmax2 = max([J_max abs(J_min)]);

J_max = max(max(max(Jac_new(:,3*Nk+1:4*Nk,2*Nk+1:3*Nk))));
J_min = min(min(min(Jac_new(:,3*Nk+1:4*Nk,2*Nk+1:3*Nk))));
climmax3 = max([J_max abs(J_min)]);

J_max = max(max(max(Jac_new(:,3*Nk+1:4*Nk,3*Nk+1:4*Nk))));
J_min = min(min(min(Jac_new(:,3*Nk+1:4*Nk,3*Nk+1:4*Nk))));
climmax4 = max([J_max abs(J_min)]);

pert_size = [0.0001 0.1 1.0 2.0 ];

xoffset = 0.04;

for i = 1:imax
    
    j = J_vec(i);
    
    J = zeros(5*Nk,5*Nk);
    J(:,:) = Jac(j,:,:);
    
    %%%%%%%%%%%%%%%%%%%
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


    line_wid = 1;
    fontsize = 8;
    
    subplot(2,imax,i)
    [C,h] = contour(Nk_vec(top-5:end),Nk_vec(top-5:end),J_Ht(top-5:end,top-5:end));
    gcapos = get(gca,'position');
    if i == 1; 
        ylabel({'\partial H /\partial \theta';'Output Index'},'FontSize',fontsize,'FontName','TimesNewRoman')
    end
    caxis([-climmax1 climmax1]);
    set(h,'LineWidth',line_wid);
    set(gca,'YDir','reverse')
    title(['\delta = ',num2str(pert_size(i))],'FontSize',fontsize,'FontName','TimesNewRoman')
    xlabel('Perturbation Index (\theta)','FontSize',fontsize,'FontName','TimesNewRoman')
    colormap jet
    if i == imax
        colorbar
    end
    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
    set(gca,'position',[gcapos(1)-xoffset gcapos(2) gcapos(3) gcapos(4)])

    levlist = get(h,'LevelList');
    indzero = find(levlist==0)
    levlist(indzero) = [];
    set(h,'LevelList',levlist)

    subplot(2,imax,imax+i)
    [C,h] = contour(Nk_vec(top-5:end),Nk_vec(top-5:end),J_Hq(top-5:end,top-5:end));
    gcapos = get(gca,'position');
    if i == 1; 
        ylabel({'\partial H /\partial q';'Output Index'},'FontSize',fontsize,'FontName','TimesNewRoman')
    end
    caxis([-climmax2 climmax2]);
    set(h,'LineWidth',line_wid);
    set(gca,'YDir','reverse')
    xlabel('Perturbation Index (q)','FontSize',fontsize,'FontName','TimesNewRoman')
    colormap jet
    if i == imax
        colorbar
    end
    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
    set(gca,'position',[gcapos(1)-xoffset gcapos(2) gcapos(3) gcapos(4)])

    levlist = get(h,'LevelList');
    indzero = find(levlist==0)
    levlist(indzero) = [];
    set(h,'LevelList',levlist)

%     subplot(2,imax,2*imax+i)
%     [C,h] = contour(Nk_vec(top-5:end),Nk_vec(top-5:end),J_Mt(top-5:end,top-5:end));
%     gcapos = get(gca,'position');
%     if i == 1; 
%         ylabel({'\partial Q /\partial \theta';'Output Index'},'FontSize',fontsize,'FontName','TimesNewRoman')
%     end
%     caxis([-climmax3 climmax3]);
%     set(h,'LineWidth',line_wid);
%     set(gca,'YDir','reverse')
%     xlabel('Perturbation Index (\theta)','FontSize',fontsize,'FontName','TimesNewRoman')
%     colormap jet
%     if i == imax
%         colorbar
%     end
%     set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
%     set(gca,'position',[gcapos(1)-xoffset gcapos(2) gcapos(3) gcapos(4)])
% 
%     subplot(2,imax,3*imax+i)
%     [C,h] = contour(Nk_vec(top-5:end),Nk_vec(top-5:end),J_Mq(top-5:end,top-5:end));
%     gcapos = get(gca,'position');
%     if i == 1; 
%         ylabel({'\partial Q /\partial q';'Output Index'},'FontSize',fontsize,'FontName','TimesNewRoman')
%     end
%     caxis([-climmax4 climmax4]);
%     set(h,'LineWidth',line_wid);
%     set(gca,'YDir','reverse')
%     xlabel('Perturbation Index (q)','FontSize',fontsize,'FontName','TimesNewRoman')
%     colormap jet
%     if i == imax
%         colorbar
%     end
%     set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
%     set(gca,'position',[gcapos(1)-xoffset gcapos(2) gcapos(3) gcapos(4)])

end

saveas(gcf,'jacobian_shallow_series_hq.eps', 'psc2')


