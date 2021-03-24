close all
clear
clc
mydir = pwd;
load coast
coast_lat = lat; clear lat
coast_lon = long; clear long
grey = 0.75;
fontsize = 14;



plot_level_ql = 70;
plot_level_qi = 37;

file_prog  = 'v000_C180.prog.eta.20140201_00z.nc4';


cd /discover/nobackup/drholdaw/btmp.25956/prog/prog_free/

lon = ncread(file_prog,'lon'); im = length(lon);
lat = ncread(file_prog,'lat'); jm = length(lat);
lev = ncread(file_prog,'lev'); lm = length(lev);

ql_free = ncread(file_prog,'qltot');
qi_free = ncread(file_prog,'qitot');

cd /discover/nobackup/drholdaw/btmp.25956/

ql_tanh = ncread(file_prog,'qltot');
qi_tanh = ncread(file_prog,'qitot');

cd(mydir)

ql_free_plot = ql_free(:,:,plot_level_ql);
ql_tanh_plot = ql_tanh(:,:,plot_level_ql);

max_ql_plot = max([ql_free_plot(:); ql_tanh_plot(:) ]);
max_ql_step = max_ql_plot/20;


qi_free_plot = qi_free(:,:,plot_level_qi);
qi_tanh_plot = qi_tanh(:,:,plot_level_qi);

max_qi_plot = max([qi_free_plot(:); qi_tanh_plot(:) ]);
max_qi_step = max_qi_plot/20;

mean_ql_pdf  = zeros(1,lm);
mean_ql_tanh = zeros(1,lm);
mean_qi_pdf  = zeros(1,lm);
mean_qi_tanh = zeros(1,lm);

for k = 1:lm

    mean_ql_pdf(k)  = mean(mean(ql_free(:,:,k)));
    mean_ql_tanh(k) = mean(mean(ql_tanh(:,:,k)));
    mean_qi_pdf(k)  = mean(mean(qi_free(:,:,k)));
    mean_qi_tanh(k) = mean(mean(qi_tanh(:,:,k)));
    
end

figure
subplot(1,2,1)
plot(mean_ql_pdf,1:lm,'b')
hold on
plot(mean_ql_tanh,1:lm,'r')
set(gca,'Ydir','reverse')
subplot(1,2,2)
plot(mean_ql_pdf./mean_ql_tanh,1:lm,'b')
set(gca,'Ydir','reverse')

figure
subplot(1,2,1)
plot(mean_qi_pdf,1:lm,'b')
hold on
plot(mean_qi_tanh,1:lm,'r')
set(gca,'Ydir','reverse')
subplot(1,2,2)
plot(mean_qi_pdf./mean_ql_tanh,1:lm,'b')
set(gca,'Ydir','reverse')




figure
set(gcf,'position',[1305 44 1000 800])
set(gcf, 'Color', 'None')
set(gcf,'Renderer','Zbuffer')

subplot(3,1,1)
contour(lon,lat,ql_free_plot','LevelStep',max_ql_step)
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([0 max_ql_plot])
colorbar
title('Specific humidity (kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Longitude','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Latitude','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(3,1,2)
contour(lon,lat,ql_tanh_plot','LevelStep',max_ql_step)
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([0 max_ql_plot])
colorbar
title('Specific humidity (kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Longitude','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Latitude','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(3,1,3)
contour(lon,lat,ql_free_plot'-ql_tanh_plot','LevelStep',max_ql_step)
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-max_ql_plot max_ql_plot])
colorbar
title('Specific humidity (kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Longitude','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Latitude','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')


figure
set(gcf,'position',[1305 44 1000 800])
set(gcf, 'Color', 'None')
set(gcf,'Renderer','Zbuffer')

subplot(3,1,1)
contour(lon,lat,qi_free_plot','LevelStep',max_ql_step)
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([0 max_qi_plot])
colorbar
title('Specific humidity (kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Longitude','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Latitude','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(3,1,2)
contour(lon,lat,qi_tanh_plot','LevelStep',max_ql_step)
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([0 max_qi_plot])
colorbar
title('Specific humidity (kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Longitude','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Latitude','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(3,1,3)
contour(lon,lat,qi_free_plot'-qi_tanh_plot','LevelStep',max_ql_step)
hold on
plot(coast_lon,coast_lat,'Color',[grey grey grey])
caxis([-max_qi_plot max_qi_plot])
colorbar
title('Specific humidity (kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Longitude','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Latitude','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')