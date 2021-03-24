close all
clear
clc

cd /discover/nobackup/drholdaw/tmp.22292/sens.20130117.000000/

jmax = 4;
time = zeros(1,jmax);
increase = zeros(8,jmax);
maxival = zeros(8,jmax+1);

rmax = zeros(8,jmax+1);
cmax = zeros(8,jmax+1);
lmax = zeros(8,jmax+1);

time(5) = 1900;
time(4) = time(5) + 15;
time(3) = time(5) + 30;
time(2) = time(5) + 45;
time(1) = time(5) + 100;

for j = 1:jmax

    if j == 1
        file1 = ['fsens.eta.20130116_',num2str(time(1)),'z_sbac.nc4'];
        file2 = ['fsens.eta.20130116_',num2str(time(2)),'z_sbac.nc4'];
    elseif j == 2
        file1 = ['fsens.eta.20130116_',num2str(time(2)),'z_sbac.nc4'];
        file2 = ['fsens.eta.20130116_',num2str(time(3)),'z_sbac.nc4'];
    elseif j == 3
        file1 = ['fsens.eta.20130116_',num2str(time(3)),'z_sbac.nc4'];
        file2 = ['fsens.eta.20130116_',num2str(time(4)),'z_sbac.nc4'];
    elseif j == 4
        file1 = ['fsens.eta.20130116_',num2str(time(4)),'z_sbac.nc4'];
        file2 = ['fsens.eta.20130116_',num2str(time(5)),'z_sbac.nc4'];
    end

    u_a = ncread(file1,'u');
    v_a = ncread(file1,'v');
    t_a = ncread(file1,'tv');
    q_a = ncread(file1,'sphu');
    p_a = ncread(file1,'delp');
    qi_a = ncread(file1,'qitot');
    ql_a = ncread(file1,'qltot');
    o3_a = ncread(file1,'ozone');
    
    u_b = ncread(file2,'u');
    v_b = ncread(file2,'v');
    t_b = ncread(file2,'tv');
    q_b = ncread(file2,'sphu');
    p_b = ncread(file2,'delp');
    qi_b = ncread(file2,'qitot');
    ql_b = ncread(file2,'qltot');
    o3_b = ncread(file2,'ozone');
    
    increase(1,j) = max(abs(u_b(:)))/max(abs(u_a(:)));
    increase(2,j) = max(abs(v_b(:)))/max(abs(v_a(:)));
    increase(3,j) = max(abs(t_b(:)))/max(abs(t_a(:)));
    increase(4,j) = max(abs(q_b(:)))/max(abs(q_a(:)));
    increase(5,j) = max(abs(p_b(:)))/max(abs(p_a(:)));
    increase(6,j) = max(abs(qi_b(:)))/max(abs(qi_a(:)));
    increase(7,j) = max(abs(ql_b(:)))/max(abs(ql_a(:)));
    increase(8,j) = max(abs(o3_b(:)))/max(abs(o3_a(:)));
    
    [maxival(1,j) location] = max(abs(u_a(:)));
    [rmax(1,j),cmax(1,j),lmax(1,j)] = ind2sub(size(u_a),location);
    [maxival(2,j) location] = max(abs(v_a(:)));
    [rmax(2,j),cmax(2,j),lmax(2,j)] = ind2sub(size(v_a),location);
    [maxival(3,j) location] = max(abs(t_a(:)));
    [rmax(3,j),cmax(3,j),lmax(3,j)] = ind2sub(size(t_a),location);
    [maxival(4,j) location] = max(abs(q_a(:)));
    [rmax(4,j),cmax(4,j),lmax(4,j)] = ind2sub(size(q_a),location);
    [maxival(5,j) location] = max(abs(p_a(:)));
    [rmax(5,j),cmax(5,j),lmax(5,j)] = ind2sub(size(p_a),location);
    [maxival(6,j) location] = max(abs(qi_a(:)));
    [rmax(6,j),cmax(6,j),lmax(6,j)] = ind2sub(size(qi_a),location);
    [maxival(7,j) location] = max(abs(ql_a(:)));
    [rmax(7,j),cmax(7,j),lmax(7,j)] = ind2sub(size(ql_a),location);
    [maxival(8,j) location] = max(abs(o3_a(:)));
    [rmax(8,j),cmax(8,j),lmax(8,j)] = ind2sub(size(o3_a),location);
    
    [maxival(1,j+1) location] = max(abs(u_b(:)));
    [rmax(1,j+1),cmax(1,j+1),lmax(1,j+1)] = ind2sub(size(u_a),location);
    [maxival(2,j+1) location] = max(abs(v_b(:)));
    [rmax(2,j+1),cmax(2,j+1),lmax(2,j+1)] = ind2sub(size(v_a),location);
    [maxival(3,j+1) location] = max(abs(t_b(:)));
    [rmax(3,j+1),cmax(3,j+1),lmax(3,j+1)] = ind2sub(size(t_a),location);
    [maxival(4,j+1) location] = max(abs(q_b(:)));
    [rmax(4,j+1),cmax(4,j+1),lmax(4,j+1)] = ind2sub(size(q_a),location);
    [maxival(5,j+1) location] = max(abs(p_b(:)));
    [rmax(5,j+1),cmax(5,j+1),lmax(5,j+1)] = ind2sub(size(p_a),location);
    [maxival(6,j+1) location] = max(abs(qi_b(:)));
    [rmax(6,j+1),cmax(6,j+1),lmax(6,j+1)] = ind2sub(size(qi_a),location);
    [maxival(7,j+1) location] = max(abs(ql_b(:)));
    [rmax(7,j+1),cmax(7,j+1),lmax(7,j+1)] = ind2sub(size(ql_a),location);
    [maxival(8,j+1) location] = max(abs(o3_b(:)));
    [rmax(8,j+1),cmax(8,j+1),lmax(8,j+1)] = ind2sub(size(o3_a),location);
    
end



figure
set(gcf,'position',[349 30 930 889])

subplot(4,2,1)
plot(1:jmax,increase(1,:))
hold on
plot(1:jmax,ones(1,jmax),'k--')
title('u')

subplot(4,2,2)
plot(1:jmax,increase(2,:))
hold on
plot(1:jmax,ones(1,jmax),'k--')
title('v')

subplot(4,2,3)
plot(1:jmax,increase(3,:))
hold on
plot(1:jmax,ones(1,jmax),'k--')
title('T')

subplot(4,2,4)
plot(1:jmax,increase(4,:))
hold on
plot(1:jmax,ones(1,jmax),'k--')
title('q')

subplot(4,2,5)
plot(1:jmax,increase(5,:))
hold on
plot(1:jmax,ones(1,jmax),'k--')
title('p')

subplot(4,2,6)
plot(1:jmax,increase(6,:))
hold on
plot(1:jmax,ones(1,jmax),'k--')
title('Qi')

subplot(4,2,7)
plot(1:jmax,increase(7,:))
hold on
plot(1:jmax,ones(1,jmax),'k--')
title('Ql')

subplot(4,2,8)
plot(1:jmax,increase(8,:))
hold on
plot(1:jmax,ones(1,jmax),'k--')
title('O_3')






figure
set(gcf,'position',[349 30 930 889])

subplot(4,2,1)
plot(1:jmax+1,maxival(1,:))
title('u')

subplot(4,2,2)
plot(1:jmax+1,maxival(2,:))
title('v')

subplot(4,2,3)
plot(1:jmax+1,maxival(3,:))
title('T')

subplot(4,2,4)
plot(1:jmax+1,maxival(4,:))
title('q')

subplot(4,2,5)
plot(1:jmax+1,maxival(5,:))
title('p')

subplot(4,2,6)
plot(1:jmax+1,maxival(6,:))
title('Qi')

subplot(4,2,7)
plot(1:jmax+1,maxival(7,:))
title('Ql')

subplot(4,2,8)
plot(1:jmax+1,maxival(8,:))
title('O_3')
