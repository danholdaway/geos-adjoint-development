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

qi_a = ncread(file,'qitot');
ql_a = ncread(file,'qltot');

cd /discover/nobackup/drholdaw/tmp.22293/
file = 'x0011dh_a.prog.eta.20130116_00z.nc4';

qi_b = ncread(file,'qitot');
ql_b = ncread(file,'qltot');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

qi_a_lev1 = qi_a(:,:,plot_level1);
qi_b_lev1 = qi_b(:,:,plot_level1);
qi_a_lev2 = qi_a(:,:,plot_level2);
qi_b_lev2 = qi_b(:,:,plot_level2);
qi_a_lev3 = qi_a(:,:,plot_level3);
qi_b_lev3 = qi_b(:,:,plot_level3);
qi_a_lev4 = qi_a(:,:,plot_level4);
qi_b_lev4 = qi_b(:,:,plot_level4);

ql_a_lev1 = ql_a(:,:,plot_level1);
ql_b_lev1 = ql_b(:,:,plot_level1);
ql_a_lev2 = ql_a(:,:,plot_level2);
ql_b_lev2 = ql_b(:,:,plot_level2);
ql_a_lev3 = ql_a(:,:,plot_level3);
ql_b_lev3 = ql_b(:,:,plot_level3);
ql_a_lev4 = ql_a(:,:,plot_level4);
ql_b_lev4 = ql_b(:,:,plot_level4);

minqi1 = min([qi_a_lev1(:); qi_b_lev1(:)]);
maxqi1 = max([qi_a_lev1(:); qi_b_lev1(:)]);
minqi2 = min([qi_a_lev2(:); qi_b_lev2(:)]);
maxqi2 = max([qi_a_lev2(:); qi_b_lev2(:)]);
minqi3 = min([qi_a_lev3(:); qi_b_lev3(:)]);
maxqi3 = max([qi_a_lev3(:); qi_b_lev3(:)]);
minqi4 = min([qi_a_lev4(:); qi_b_lev4(:)]);
maxqi4 = max([qi_a_lev4(:); qi_b_lev4(:)]);

minql1 = min([ql_a_lev1(:); ql_b_lev1(:)]);
maxql1 = max([ql_a_lev1(:); ql_b_lev1(:)]);
minql2 = min([ql_a_lev2(:); ql_b_lev2(:)]);
maxql2 = max([ql_a_lev2(:); ql_b_lev2(:)]);
minql3 = min([ql_a_lev3(:); ql_b_lev3(:)]);
maxql3 = max([ql_a_lev3(:); ql_b_lev3(:)]);
minql4 = min([ql_a_lev4(:); ql_b_lev4(:)]);
maxql4 = max([ql_a_lev4(:); ql_b_lev4(:)]);

figure
set(gcf,'position',[349 30 930 889])
    
subplot(4,2,1)
contour(lon,lat,((qi_a_lev1))','LevelStep',(maxqi1-minqi1)/5)
caxis([minqi1 maxqi1])
title(['qi Level ' num2str(plot_level1)])
colorbar

subplot(4,2,2)
contour(lon,lat,((qi_b_lev1))','LevelStep',(maxqi1-minqi1)/5)
caxis([minqi1 maxqi1])
title(['qi Level ' num2str(plot_level1)])
colorbar
    
subplot(4,2,3)
contour(lon,lat,((qi_a_lev2))','LevelStep',(maxqi2-minqi2)/5)
caxis([minqi2 maxqi2])
title(['qi Level ' num2str(plot_level2)])
colorbar

subplot(4,2,4)
contour(lon,lat,((qi_b_lev2))','LevelStep',(maxqi2-minqi2)/5)
caxis([minqi2 maxqi2])
title(['qi Level ' num2str(plot_level2)])
colorbar
    
subplot(4,2,5)
contour(lon,lat,((qi_a_lev3))','LevelStep',(maxqi3-minqi3)/5)
caxis([minqi3 maxqi3])
title(['qi Level ' num2str(plot_level3)])
colorbar

subplot(4,2,6)
contour(lon,lat,((qi_b_lev3))','LevelStep',(maxqi3-minqi3)/5)
caxis([minqi3 maxqi3])
title(['qi Level ' num2str(plot_level3)])
colorbar
    
subplot(4,2,7)
contour(lon,lat,((qi_a_lev4))','LevelStep',(maxqi4-minqi4)/5)
caxis([minqi4 maxqi4])
title(['qi Level ' num2str(plot_level4)])
colorbar

subplot(4,2,8)
contour(lon,lat,((qi_b_lev4))','LevelStep',(maxqi4-minqi4)/5)
caxis([minqi4 maxqi4])
title(['qi Level ' num2str(plot_level4)])
colorbar





figure
set(gcf,'position',[1283 55 930 889])

subplot(4,2,1)
contour(lon,lat,((ql_a_lev1))','LevelStep',(maxql1-minql1)/5)
caxis([minql1 maxql1])
title(['ql Level ' num2str(plot_level1)])
colorbar

subplot(4,2,2)
contour(lon,lat,((ql_b_lev1))','LevelStep',(maxql1-minql1)/5)
caxis([minql1 maxql1])
title(['ql Level ' num2str(plot_level1)])
colorbar
    
subplot(4,2,3)
contour(lon,lat,((ql_a_lev2))','LevelStep',(maxql2-minql2)/5)
caxis([minql2 maxql2])
title(['ql Level ' num2str(plot_level2)])
colorbar

subplot(4,2,4)
contour(lon,lat,((ql_b_lev2))','LevelStep',(maxql2-minql2)/5)
caxis([minql2 maxql2])
title(['ql Level ' num2str(plot_level2)])
colorbar
    
subplot(4,2,5)
contour(lon,lat,((ql_a_lev3))','LevelStep',(maxql3-minql3)/5)
caxis([minql3 maxql3])
title(['ql Level ' num2str(plot_level3)])
colorbar

subplot(4,2,6)
contour(lon,lat,((ql_b_lev3))','LevelStep',(maxql3-minql3)/5)
caxis([minql3 maxql3])
title(['ql Level ' num2str(plot_level3)])
colorbar
    
subplot(4,2,7)
contour(lon,lat,((ql_a_lev4))','LevelStep',(maxql4-minql4)/5)
caxis([minql4 maxql4])
title(['ql Level ' num2str(plot_level4)])
colorbar

subplot(4,2,8)
contour(lon,lat,((ql_b_lev4))','LevelStep',(maxql4-minql4)/5)
caxis([minql4 maxql4])
title(['ql Level ' num2str(plot_level4)])
colorbar


