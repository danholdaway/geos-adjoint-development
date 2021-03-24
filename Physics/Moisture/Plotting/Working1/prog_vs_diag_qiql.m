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

qls_a = ncread(file,'qils') + ncread(file,'qlls');
qcn_a = ncread(file,'qicn') + ncread(file,'qlcn');

cd /discover/nobackup/drholdaw/tmp.22293/
file = 'x0011dh_a.prog.eta.20130116_00z_rh102.nc4';

qls_b = ncread(file,'qils1') + ncread(file,'qlls1');
qcn_b = ncread(file,'qicn1') + ncread(file,'qlcn1');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

qls_a_lev1 = qls_a(:,:,plot_level1);
qls_b_lev1 = qls_b(:,:,plot_level1);
qls_a_lev2 = qls_a(:,:,plot_level2);
qls_b_lev2 = qls_b(:,:,plot_level2);
qls_a_lev3 = qls_a(:,:,plot_level3);
qls_b_lev3 = qls_b(:,:,plot_level3);
qls_a_lev4 = qls_a(:,:,plot_level4);
qls_b_lev4 = qls_b(:,:,plot_level4);

qcn_a_lev1 = qcn_a(:,:,plot_level1);
qcn_b_lev1 = qcn_b(:,:,plot_level1);
qcn_a_lev2 = qcn_a(:,:,plot_level2);
qcn_b_lev2 = qcn_b(:,:,plot_level2);
qcn_a_lev3 = qcn_a(:,:,plot_level3);
qcn_b_lev3 = qcn_b(:,:,plot_level3);
qcn_a_lev4 = qcn_a(:,:,plot_level4);
qcn_b_lev4 = qcn_b(:,:,plot_level4);

minqi1 = min([qls_a_lev1(:); qls_b_lev1(:)]);
maxqi1 = max([qls_a_lev1(:); qls_b_lev1(:)]);
minqi2 = min([qls_a_lev2(:); qls_b_lev2(:)]);
maxqi2 = max([qls_a_lev2(:); qls_b_lev2(:)]);
minqi3 = min([qls_a_lev3(:); qls_b_lev3(:)]);
maxqi3 = max([qls_a_lev3(:); qls_b_lev3(:)]);
minqi4 = min([qls_a_lev4(:); qls_b_lev4(:)]);
maxqi4 = max([qls_a_lev4(:); qls_b_lev4(:)]);

minql1 = min([qcn_a_lev1(:); qcn_b_lev1(:)]);
maxql1 = max([qcn_a_lev1(:); qcn_b_lev1(:)]);
minql2 = min([qcn_a_lev2(:); qcn_b_lev2(:)]);
maxql2 = max([qcn_a_lev2(:); qcn_b_lev2(:)]);
minql3 = min([qcn_a_lev3(:); qcn_b_lev3(:)]);
maxql3 = max([qcn_a_lev3(:); qcn_b_lev3(:)]);
minql4 = min([qcn_a_lev4(:); qcn_b_lev4(:)]);
maxql4 = max([qcn_a_lev4(:); qcn_b_lev4(:)]);

figure
set(gcf,'position',[349 30 930 889])
    
subplot(4,2,1)
contour(lon,lat,((qls_a_lev1))','LevelStep',(maxqi1-minqi1)/15)
caxis([minqi1 maxqi1])
title(['qls Level ' num2str(plot_level1)])
colorbar

subplot(4,2,2)
contour(lon,lat,((qls_b_lev1))','LevelStep',(maxqi1-minqi1)/15)
caxis([minqi1 maxqi1])
title(['qls Level ' num2str(plot_level1)])
colorbar
    
subplot(4,2,3)
contour(lon,lat,((qls_a_lev2))','LevelStep',(maxqi2-minqi2)/15)
caxis([minqi2 maxqi2])
title(['qls Level ' num2str(plot_level2)])
colorbar

subplot(4,2,4)
contour(lon,lat,((qls_b_lev2))','LevelStep',(maxqi2-minqi2)/15)
caxis([minqi2 maxqi2])
title(['qls Level ' num2str(plot_level2)])
colorbar
    
subplot(4,2,5)
contour(lon,lat,((qls_a_lev3))','LevelStep',(maxqi3-minqi3)/15)
caxis([minqi3 maxqi3])
title(['qls Level ' num2str(plot_level3)])
colorbar

subplot(4,2,6)
contour(lon,lat,((qls_b_lev3))','LevelStep',(maxqi3-minqi3)/15)
caxis([minqi3 maxqi3])
title(['qls Level ' num2str(plot_level3)])
colorbar
    
subplot(4,2,7)
contour(lon,lat,((qls_a_lev4))','LevelStep',(maxqi4-minqi4)/15)
caxis([minqi4 maxqi4])
title(['qls Level ' num2str(plot_level4)])
colorbar

subplot(4,2,8)
contour(lon,lat,((qls_b_lev4))','LevelStep',(maxqi4-minqi4)/15)
caxis([minqi4 maxqi4])
title(['qls Level ' num2str(plot_level4)])
colorbar





figure
set(gcf,'position',[1283 55 930 889])

subplot(4,2,1)
contour(lon,lat,((qcn_a_lev1))','LevelStep',(maxql1-minql1)/15)
caxis([minql1 maxql1])
title(['qcn Level ' num2str(plot_level1)])
colorbar

subplot(4,2,2)
contour(lon,lat,((qcn_b_lev1))','LevelStep',(maxql1-minql1)/15)
caxis([minql1 maxql1])
title(['qcn Level ' num2str(plot_level1)])
colorbar
    
subplot(4,2,3)
contour(lon,lat,((qcn_a_lev2))','LevelStep',(maxql2-minql2)/15)
caxis([minql2 maxql2])
title(['qcn Level ' num2str(plot_level2)])
colorbar

subplot(4,2,4)
contour(lon,lat,((qcn_b_lev2))','LevelStep',(maxql2-minql2)/15)
caxis([minql2 maxql2])
title(['qcn Level ' num2str(plot_level2)])
colorbar
    
subplot(4,2,5)
contour(lon,lat,((qcn_a_lev3))','LevelStep',(maxql3-minql3)/15)
caxis([minql3 maxql3])
title(['qcn Level ' num2str(plot_level3)])
colorbar

subplot(4,2,6)
contour(lon,lat,((qcn_b_lev3))','LevelStep',(maxql3-minql3)/15)
caxis([minql3 maxql3])
title(['qcn Level ' num2str(plot_level3)])
colorbar
    
subplot(4,2,7)
contour(lon,lat,((qcn_a_lev4))','LevelStep',(maxql4-minql4)/15)
caxis([minql4 maxql4])
title(['qcn Level ' num2str(plot_level4)])
colorbar

subplot(4,2,8)
contour(lon,lat,((qcn_b_lev4))','LevelStep',(maxql4-minql4)/15)
caxis([minql4 maxql4])
title(['qcn Level ' num2str(plot_level4)])
colorbar


