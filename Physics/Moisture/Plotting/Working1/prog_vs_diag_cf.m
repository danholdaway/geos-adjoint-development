close all
clc
clear

%Choose model level to plot.
plot_level1 = 68;
plot_level2 = 60;
plot_level3 = 50;
plot_level4 = 40;

%Choose whether to match the scale.
match_scale = 1;

cd /discover/nobackup/drholdaw/tmp.22292/prog/prog_free/
file = 'x0011dh_a.prog.eta.20130116_00z.nc4';

lon = ncread(file,'lon');
lat = ncread(file,'lat');

cfls_a = ncread(file,'cfls');
cfcn_a = ncread(file,'cfcn');

cd /discover/nobackup/drholdaw/tmp.22293/
file = 'x0011dh_a.prog.eta.20130116_00z.nc4';

cfls_b = ncread(file,'cfls');
cfcn_b = ncread(file,'cfcn');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

cfls_a_lev1 = cfls_a(:,:,plot_level1);
cfls_b_lev1 = cfls_b(:,:,plot_level1);
cfls_a_lev2 = cfls_a(:,:,plot_level2);
cfls_b_lev2 = cfls_b(:,:,plot_level2);
cfls_a_lev3 = cfls_a(:,:,plot_level3);
cfls_b_lev3 = cfls_b(:,:,plot_level3);
cfls_a_lev4 = cfls_a(:,:,plot_level4);
cfls_b_lev4 = cfls_b(:,:,plot_level4);

cfcn_a_lev1 = cfcn_a(:,:,plot_level1);
cfcn_b_lev1 = cfcn_b(:,:,plot_level1);
cfcn_a_lev2 = cfcn_a(:,:,plot_level2);
cfcn_b_lev2 = cfcn_b(:,:,plot_level2);
cfcn_a_lev3 = cfcn_a(:,:,plot_level3);
cfcn_b_lev3 = cfcn_b(:,:,plot_level3);
cfcn_a_lev4 = cfcn_a(:,:,plot_level4);
cfcn_b_lev4 = cfcn_b(:,:,plot_level4);

mincfls1 = min([cfls_a_lev1(:); cfls_b_lev1(:)]);
maxcfls1 = max([cfls_a_lev1(:); cfls_b_lev1(:)]);
mincfls2 = min([cfls_a_lev2(:); cfls_b_lev2(:)]);
maxcfls2 = max([cfls_a_lev2(:); cfls_b_lev2(:)]);
mincfls3 = min([cfls_a_lev3(:); cfls_b_lev3(:)]);
maxcfls3 = max([cfls_a_lev3(:); cfls_b_lev3(:)]);
mincfls4 = min([cfls_a_lev4(:); cfls_b_lev4(:)]);
maxcfls4 = max([cfls_a_lev4(:); cfls_b_lev4(:)]);

mincfcn1 = min([cfcn_a_lev1(:); cfcn_b_lev1(:)]);
maxcfcn1 = max([cfcn_a_lev1(:); cfcn_b_lev1(:)]);
mincfcn2 = min([cfcn_a_lev2(:); cfcn_b_lev2(:)]);
maxcfcn2 = max([cfcn_a_lev2(:); cfcn_b_lev2(:)]);
mincfcn3 = min([cfcn_a_lev3(:); cfcn_b_lev3(:)]);
maxcfcn3 = max([cfcn_a_lev3(:); cfcn_b_lev3(:)]);
mincfcn4 = min([cfcn_a_lev4(:); cfcn_b_lev4(:)]);
maxcfcn4 = max([cfcn_a_lev4(:); cfcn_b_lev4(:)]);

figure
set(gcf,'position',[349 30 930 889])
    
subplot(4,2,1)
contour(lon,lat,((cfls_a_lev1))','LevelStep',(maxcfls1-mincfls1)/5)
caxis([mincfls1 maxcfls1])
title(['cfls Level ' num2str(plot_level1)])
colorbar

subplot(4,2,2)
contour(lon,lat,((cfls_b_lev1))','LevelStep',(maxcfls1-mincfls1)/5)
caxis([mincfls1 maxcfls1])
title(['cfls Level ' num2str(plot_level1)])
colorbar
    
subplot(4,2,3)
contour(lon,lat,((cfls_a_lev2))','LevelStep',(maxcfls2-mincfls2)/5)
caxis([mincfls2 maxcfls2])
title(['cfls Level ' num2str(plot_level2)])
colorbar

subplot(4,2,4)
contour(lon,lat,((cfls_b_lev2))','LevelStep',(maxcfls2-mincfls2)/5)
caxis([mincfls2 maxcfls2])
title(['cfls Level ' num2str(plot_level2)])
colorbar
    
subplot(4,2,5)
contour(lon,lat,((cfls_a_lev3))','LevelStep',(maxcfls3-mincfls3)/5)
caxis([mincfls3 maxcfls3])
title(['cfls Level ' num2str(plot_level3)])
colorbar

subplot(4,2,6)
contour(lon,lat,((cfls_b_lev3))','LevelStep',(maxcfls3-mincfls3)/5)
caxis([mincfls3 maxcfls3])
title(['cfls Level ' num2str(plot_level3)])
colorbar
    
subplot(4,2,7)
contour(lon,lat,((cfls_a_lev4))','LevelStep',(maxcfls4-mincfls4)/5)
caxis([mincfls4 maxcfls4])
title(['cfls Level ' num2str(plot_level4)])
colorbar

subplot(4,2,8)
contour(lon,lat,((cfls_b_lev4))','LevelStep',(maxcfls4-mincfls4)/5)
caxis([mincfls4 maxcfls4])
title(['cfls Level ' num2str(plot_level4)])
colorbar





figure
set(gcf,'position',[1283 55 930 889])

subplot(4,2,1)
contour(lon,lat,((cfcn_a_lev1))','LevelStep',(maxcfcn1-mincfcn1)/5)
caxis([mincfcn1 maxcfcn1])
title(['cfcn Level ' num2str(plot_level1)])
colorbar

subplot(4,2,2)
contour(lon,lat,((cfcn_b_lev1))')%,'LevelStep',(maxcfcn1-mincfcn1)/5)
% caxis([mincfcn1 maxcfcn1])
title(['cfcn Level ' num2str(plot_level1)])
colorbar
    
subplot(4,2,3)
contour(lon,lat,((cfcn_a_lev2))','LevelStep',(maxcfcn2-mincfcn2)/5)
caxis([mincfcn2 maxcfcn2])
title(['cfcn Level ' num2str(plot_level2)])
colorbar

subplot(4,2,4)
contour(lon,lat,((cfcn_b_lev2))')%,'LevelStep',(maxcfcn2-mincfcn2)/5)
% caxis([mincfcn2 maxcfcn2])
title(['cfcn Level ' num2str(plot_level2)])
colorbar
    
subplot(4,2,5)
contour(lon,lat,((cfcn_a_lev3))','LevelStep',(maxcfcn3-mincfcn3)/5)
caxis([mincfcn3 maxcfcn3])
title(['cfcn Level ' num2str(plot_level3)])
colorbar

subplot(4,2,6)
contour(lon,lat,((cfcn_b_lev3))')%,'LevelStep',(maxcfcn3-mincfcn3)/5)
% caxis([mincfcn3 maxcfcn3])
title(['cfcn Level ' num2str(plot_level3)])
colorbar
    
subplot(4,2,7)
contour(lon,lat,((cfcn_a_lev4))','LevelStep',(maxcfcn4-mincfcn4)/5)
caxis([mincfcn4 maxcfcn4])
title(['cfcn Level ' num2str(plot_level4)])
colorbar

subplot(4,2,8)
contour(lon,lat,((cfcn_b_lev4))')%,'LevelStep',(maxcfcn4-mincfcn4)/5)
% caxis([mincfcn4 maxcfcn4])
title(['cfcn Level ' num2str(plot_level4)])
colorbar


