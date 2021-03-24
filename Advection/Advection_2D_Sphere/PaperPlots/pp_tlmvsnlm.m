close all
clear
clc
mydir = pwd;

figvis = 'on';

rehash
mo_Figure4end
cd(mydir)

lat = -90:181/90:90;
lon = -179:2:180;

posmult = 1.155555;
loweroffset = 0.00;
rightoffset = 0.00;
fontsize = 12;

qnlm = qreplutlim'-qfreeutlim';
qtlm = qtlmutlim';

qmax = max(abs([qnlm(:); qtlm(:)]));
qmin = -qmax;
qstep = qmax/10;


coli = 255;
colvec = 0:1/(coli-1):1;
fac = 1;

blue = horzcat(colvec.^fac',colvec.^fac',ones(coli,1));
red = horzcat(ones(coli,1),colvec.^fac',colvec.^fac');

cmap_sens_grey_p = colormap(gray);
cmap_sens_grey_n = colormap(flipud(gray));
cmap_sens_grey = [cmap_sens_grey_p;cmap_sens_grey_n];

cmap_sens_colo_p = blue;
cmap_sens_colo_n = flipud(red);
cmap_sens_colo = [cmap_sens_colo_p;cmap_sens_colo_n];

close

cints = 10;

for i = 1:2

    if i == 1
        C1 = figure('visible',figvis);
        colormap(cmap_sens_grey)
    else
        C2 = figure('visible',figvis);
        colormap(cmap_sens_colo)
    end
    set(gcf,'position',[18 133 600 566])

    subplot(2,1,1)
    set(gcf,'renderer','painters');
    contourf(lon,lat,qnlm,'LevelList',qstep:qstep:qmax,'LineStyle','none')
    hold on
    contour(lon,lat,qnlm,'LevelList',qmin:qstep:-qstep)
    caxis([qmin qmax])
    axis equal
    ylim([-90 90])
    xlim([-179 180])
    colorbar
    pos = get(gca,'position');
    set(gca,'position',[pos(1)-0.02 pos(2)-0.02 posmult*pos(3) posmult*pos(4)])
    ylabel('Latitude','FontSize',fontsize,'FontName','TimesNewRoman')
    %xlabel('Longitude','FontSize',fontsize,'FontName','TimesNewRoman')
    title('(a) Nonlinear perturbation','FontSize',fontsize,'FontName','TimesNewRoman')
    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

    subplot(2,1,2)
    set(gcf,'renderer','painters');
    contourf(lon,lat,qtlm,'LevelList',qstep:qstep:qmax,'LineStyle','none')
    hold on
    contour(lon,lat,qtlm,'LevelList',qmin:qstep:-qstep,'LineStyle','--')
    caxis([qmin qmax])
    axis equal
    ylim([-90 90])
    xlim([-179 180])
    colorbar
    pos = get(gca,'position');
    set(gca,'position',[pos(1)-0.02 pos(2)-0.02 posmult*pos(3) posmult*pos(4)])
    ylabel('Latitude','FontSize',fontsize,'FontName','TimesNewRoman')
    xlabel('Longitude','FontSize',fontsize,'FontName','TimesNewRoman')
    title('(b) Tangent linear perturbation','FontSize',fontsize,'FontName','TimesNewRoman')
    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')


end

cd(mydir)

clearvars -except C1 C2
pause(0.1)
 
set(C1,'Units','Inches');
pos = get(C1,'Position');
set(C1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(C1,'tlmvsnlm.eps','-painters','-depsc2','-r300')

set(C2,'Units','Inches');
pos = get(C2,'Position');
set(C2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(C2,'tlmvsnlm_col.eps','-painters','-depsc2','-r300')
