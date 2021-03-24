close all
clear
clc
mydir = pwd;

figvis = 'off';

rehash
mo_gaussgaussIC
mo_cylindstepIC
cd(mydir)

lat = -90:181/90:90;
lon = -179:2:180;

posmult = 1.155555;
loweroffset = 0.03;
rightoffset = 0.03;
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
    set(gcf,'position',[18 133 1200 566])

    subplot(2,2,1)
    contourf(lon,lat,qfreeicgauss','LineStyle','none','LevelStep',0.05)
    axis equal
    ylim([-90 90])
    xlim([-179 180])
    caxis([0 1])
    colorbar
    pos = get(gca,'position');
    set(gca,'position',[pos(1) pos(2) posmult*pos(3) posmult*pos(4)])
    ylabel('Latitude','FontSize',fontsize,'FontName','TimesNewRoman')
    %xlabel('Longitude','FontSize',fontsize,'FontName','TimesNewRoman')
    title('(a) Gaussian reference initital condition','FontSize',fontsize,'FontName','TimesNewRoman')
    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

    subplot(2,2,3)
    contourf(lon,lat,qfreeiccylind','LineStyle','none','LevelStep',0.05)
    axis equal
    ylim([-90 90])
    xlim([-179 180])
    caxis([0 1])
    colorbar
    pos = get(gca,'position');
    set(gca,'position',[pos(1) pos(2)+loweroffset posmult*pos(3) posmult*pos(4)])
    ylabel('Latitude','FontSize',fontsize,'FontName','TimesNewRoman')
    xlabel('Longitude','FontSize',fontsize,'FontName','TimesNewRoman')
    title('(b) Slotted cylinder reference initital condition','FontSize',fontsize,'FontName','TimesNewRoman')
    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

    subplot(2,2,2)
    contourf(lon,lat,qperticgauss','LineStyle','none','LevelStep',0.05e-4)
    caxis([0 1e-4])
    axis equal
    ylim([-90 90])
    xlim([-179 180])
    colorbar
    pos = get(gca,'position');
    set(gca,'position',[pos(1)-rightoffset pos(2) posmult*pos(3) posmult*pos(4)])
    %ylabel('Latitude','FontSize',fontsize,'FontName','TimesNewRoman')
    %xlabel('Longitude','FontSize',fontsize,'FontName','TimesNewRoman')
    title('(c) Gaussian perturbation initital condition','FontSize',fontsize,'FontName','TimesNewRoman')
    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

    subplot(2,2,4)
    contourf(lon,lat,qperticstep','LineStyle','none','LevelStep',0.05e-4)
    caxis([0 1e-4])
    axis equal
    ylim([-90 90])
    xlim([-179 180])
    colorbar
    pos = get(gca,'position');
    set(gca,'position',[pos(1)-rightoffset pos(2)+loweroffset posmult*pos(3) posmult*pos(4)])
    %ylabel('Latitude','FontSize',fontsize,'FontName','TimesNewRoman')
    xlabel('Longitude','FontSize',fontsize,'FontName','TimesNewRoman')
    title('(d) Step perturbation initital condition','FontSize',fontsize,'FontName','TimesNewRoman')
    set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

end

cd(mydir)
    
clearvars -except C1 C2
pause(0.1)
 
set(C1,'Units','Inches');
pos = get(C1,'Position');
set(C1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(C1,'pertic.eps','-painters','-depsc2','-r300')

set(C2,'Units','Inches');
pos = get(C2,'Position');
set(C2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(C2,'pertic_col.eps','-painters','-depsc2','-r300')
