clc
clear
close all

cd /archive/u/drholdaw/x0011dh_a/prog/Y2012/M12/D26/H15/

u_dh = ncread('x0011dh_a.fsens_twe.eta.20121226_15z+20121228_00z-20121227_00z.nc4', 'u');


cd /archive/u/rtodling/g100p1/First/prog/Y2012/M12/D26/H15/

u_rt = ncread('g100p1.fsens_twe.eta.20121226_15z+20121228_00z-20121227_00z.nc4', 'u');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/

plot_level = 60;

u_dh_p = u_dh(:,:,plot_level);
u_rt_p = u_rt(:,:,plot_level);

u = [u_dh_p; u_rt_p];

maxval = max(u(:));
minval = min(u(:));

contourf(u_dh(:,:,plot_level)')
caxis([minval maxval])
colorbar

figure
contourf(u_rt(:,:,plot_level)')
caxis([minval maxval])
colorbar