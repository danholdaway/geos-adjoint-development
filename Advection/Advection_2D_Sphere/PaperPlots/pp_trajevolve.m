close all
clear
clc
mydir = pwd;

figvis = 'off';

rehash
mo_Figure2half
qfreeuthalf = qfreeut; clear qfreeut
qtlmuthalf = qtlmut; clear qtlmut
mo_Figure2end
qfreeutend = qfreeut; clear qfreeut
qreplutend = qreplut; clear qreplut
qtlmutend = qtlmut; clear qtlmut
cd(mydir)

lat = -90:181/90:90;
lon = -179:2:180;

posmult = 1.155555;
loweroffset = 0.00;
rightoffset = 0.00;
fontsize = 12;

coli = 255;
colvec = 0:1/(coli-1):1;
fac = 1;
blue = horzcat(colvec.^fac',colvec.^fac',ones(coli,1));
red = horzcat(ones(coli,1),colvec.^fac',colvec.^fac');

for i = 1:2

    if i == 1
        C1 = figure('visible',figvis);
        colormap(flipud(gray))
    else
        C2 = figure('visible',figvis);
        colormap(flipud(red))
    end
    set(gcf,'position',[18 133 600 566])

    subplot(2,1,1)
    contourf(lon,lat,qfreeuthalf','LineStyle','none','LevelStep',0.05)
    axis equal
    ylim([-90 90])
    xlim([-179 180])
    caxis([0 1])
    colorbar
    pos = get(gca,'position');
    set(gca,'position',[pos(1)-0.02 pos(2)-0.02 posmult*pos(3) posmult*pos(4)])
    ylabel('Latitude','FontSize',fontsize,'FontName','TimesNewRoman')
    %xlabel('Longitude','FontSize',fontsize,'FontName','TimesNewRoman')
    title('(a) Solution at t = T/2','FontSize',fontsize,'FontName','TimesNewRoman')
    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

    subplot(2,1,2)
    contourf(lon,lat,qfreeutend','LineStyle','none','LevelStep',0.05)
    axis equal
    ylim([-90 90])
    xlim([-179 180])
    caxis([0 1])
    colorbar
    pos = get(gca,'position');
    set(gca,'position',[pos(1)-0.02 pos(2)-0.02 posmult*pos(3) posmult*pos(4)])
    ylabel('Latitude','FontSize',fontsize,'FontName','TimesNewRoman')
    xlabel('Longitude','FontSize',fontsize,'FontName','TimesNewRoman')
    title('(b) Solution at t = T','FontSize',fontsize,'FontName','TimesNewRoman')
    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')


end

cd(mydir)

clearvars -except C1 C2
pause(0.1)
 
set(C1,'Units','Inches');
pos = get(C1,'Position');
set(C1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(C1,'trajevolve.eps','-painters','-depsc2','-r300')

set(C2,'Units','Inches');
pos = get(C2,'Position');
set(C2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(C2,'trajevolve_col.eps','-painters','-depsc2','-r300')
