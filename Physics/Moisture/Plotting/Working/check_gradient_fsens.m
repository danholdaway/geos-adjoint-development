close all
clear
clc

%ENTER DAY AND STARTING HOUR
d2 = '20130102';
H = '21';

%03
%

daten = datenum(d2, 'yyyymmdd');

d1 = datestr(daten-1, 'yyyymmdd');
d3 = datestr(daten+1, 'yyyymmdd');

cd /discover/nobackup/drholdaw/x0011dh_a/grad_fil5_rd5_adjdef/

for i = 1:2

    file = ['x0011dh_a.fsens_twe.eta.',d1,'_',H,'z+',d3,'_00z-',d2,'_00z.nc4']
    
    lon = ncread(file,'lon');
    lat = ncread(file,'lat');
    u = ncread(file,'u');
    v = ncread(file,'v');
    t = ncread(file,'tv');
    q = ncread(file,'sphu');
    p = ncread(file,'delp');
    qi = ncread(file,'qitot');
    ql = ncread(file,'qltot');
    o3 = ncread(file,'ozone');

    pl1 = 72;
    pl2 = 62;
    pl3 = 50;
    pl4 = 30;

    u1 = u(:,:,pl1)';
    u2 = u(:,:,pl2)';
    u3 = u(:,:,pl3)';
    u4 = u(:,:,pl4)';

    v1 = v(:,:,pl1)';
    v2 = v(:,:,pl2)';
    v3 = v(:,:,pl3)';
    v4 = v(:,:,pl4)';

    t1 = t(:,:,pl1)';
    t2 = t(:,:,pl2)';
    t3 = t(:,:,pl3)';
    t4 = t(:,:,pl4)';

    q1 = q(:,:,pl1)';
    q2 = q(:,:,pl2)';
    q3 = q(:,:,pl3)';
    q4 = q(:,:,pl4)';

    p1 = p(:,:,pl1)';
    p2 = p(:,:,pl2)';
    p3 = p(:,:,pl3)';
    p4 = p(:,:,pl4)';

    figure
    if i == 1 
        set(gcf,'position',[ 3    34 2554 885])
    elseif i == 2
        set(gcf,'position',[ 1283 34 2554 885])
    end
    subplot(5,4,1)
    contourf(lon,lat,u1,'LineStyle','none')
    xlabel('u')
    caxis([-max(abs(u1(:))) max(abs(u1(:)))])
    colorbar

    subplot(5,4,2)
    contourf(lon,lat,u2,'LineStyle','none')
    colorbar
    caxis([-max(abs(u2(:))) max(abs(u2(:)))])
    subplot(5,4,3)
    contourf(lon,lat,u3,'LineStyle','none')
    colorbar
    caxis([-max(abs(u3(:))) max(abs(u3(:)))])
    subplot(5,4,4)
    contourf(lon,lat,u4,'LineStyle','none')
    colorbar
    caxis([-max(abs(u4(:))) max(abs(u4(:)))])

    subplot(5,4,5)
    contourf(lon,lat,v1,'LineStyle','none')
    xlabel('v')
    colorbar
    caxis([-max(abs(v1(:))) max(abs(v1(:)))])

    subplot(5,4,6)
    contourf(lon,lat,v2,'LineStyle','none')
    colorbar
    caxis([-max(abs(v2(:))) max(abs(v2(:)))])
    subplot(5,4,7)
    contourf(lon,lat,v3,'LineStyle','none')
    colorbar
    caxis([-max(abs(v3(:))) max(abs(v3(:)))])
    subplot(5,4,8)
    contourf(lon,lat,v4,'LineStyle','none')
    colorbar
    caxis([-max(abs(v4(:))) max(abs(v4(:)))])

    subplot(5,4,9)
    contourf(lon,lat,t1,'LineStyle','none')
    xlabel('t')
    colorbar
    caxis([-max(abs(t1(:))) max(abs(t1(:)))])

    subplot(5,4,10)
    contourf(lon,lat,t2,'LineStyle','none')
    colorbar
    caxis([-max(abs(t2(:))) max(abs(t2(:)))])
    subplot(5,4,11)
    contourf(lon,lat,t3,'LineStyle','none')
    colorbar
    caxis([-max(abs(t3(:))) max(abs(t3(:)))])
    subplot(5,4,12)
    contourf(lon,lat,t4,'LineStyle','none')
    colorbar
    caxis([-max(abs(t4(:))) max(abs(t4(:)))])

    subplot(5,4,13)
    contourf(lon,lat,q1,'LineStyle','none')
    xlabel('q')
    colorbar
    caxis([-max(abs(q1(:))) max(abs(q1(:)))])

    subplot(5,4,14)
    contourf(lon,lat,q2,'LineStyle','none')
    colorbar
    caxis([-max(abs(q2(:))) max(abs(q2(:)))])
    subplot(5,4,15)
    contourf(lon,lat,q3,'LineStyle','none')
    colorbar
    caxis([-max(abs(q3(:))) max(abs(q3(:)))])
    subplot(5,4,16)
    contourf(lon,lat,q4,'LineStyle','none')
    colorbar
    caxis([-max(abs(q4(:))) max(abs(q4(:)))])


    subplot(5,4,17)
    contourf(lon,lat,p1,'LineStyle','none')
    xlabel('p')
    colorbar
    caxis([-max(abs(p1(:))) max(abs(p1(:)))])

    subplot(5,4,18)
    contourf(lon,lat,p2,'LineStyle','none')
    colorbar
    caxis([-max(abs(p2(:))) max(abs(p2(:)))])
    subplot(5,4,19)
    contourf(lon,lat,p3,'LineStyle','none')
    colorbar
    caxis([-max(abs(p3(:))) max(abs(p3(:)))])
    subplot(5,4,20)
    contourf(lon,lat,p4,'LineStyle','none')
    colorbar
    caxis([-max(abs(p4(:))) max(abs(p4(:)))])

end

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/
