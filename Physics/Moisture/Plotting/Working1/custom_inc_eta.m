close all
clear
clc

cd /discover/nobackup/drholdaw/tmp.22292/
file1 = 'inc.eta.cloud.nc4';

lon = ncread(file1,'lon');
lat = ncread(file1,'lat');

u_a = ncread(file1,'u');
v_a = ncread(file1,'v');
t_a = ncread(file1,'tv');
q_a = ncread(file1,'sphu');
p_a = ncread(file1,'delp');
qi_a = ncread(file1,'qitot');
ql_a = ncread(file1,'qltot');
o3_a = ncread(file1,'ozone');

cd /discover/nobackup/drholdaw/tmp.22292/nlm_runs/free/
file1 = 'x0011dh_a.prog.eta.20130116_00z.nc4';

u_b = ncread(file1,'u');
v_b = ncread(file1,'v');
t_b = ncread(file1,'tv');
q_b = ncread(file1,'sphu');
p_b = ncread(file1,'delp');
qi_b = ncread(file1,'qitot');
ql_b = ncread(file1,'qltot');
o3_b = ncread(file1,'ozone');

cd /discover/nobackup/drholdaw/tmp.22292/nlm_runs/free/
file1 = 'x0011dh_a.prog.eta.20130116_06z.nc4';

u_c = ncread(file1,'u');
v_c = ncread(file1,'v');
t_c = ncread(file1,'tv');
q_c = ncread(file1,'sphu');
p_c = ncread(file1,'delp');
qi_c = ncread(file1,'qitot');
ql_c = ncread(file1,'qltot');
o3_c = ncread(file1,'ozone');

cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/

u_b = u_c - u_b;
v_b = v_c - v_b;
t_b = t_c - t_b;
q_b = q_c - q_b;
p_b = p_c - p_b;
qi_b = qi_c - qi_b;
ql_b = ql_c - ql_b;
o3_b = o3_c - o3_b;

% cd /discover/nobackup/drholdaw/tmp.22292/
% ncwrite('inc.eta.cloud.nc4','u',u_b);
% ncwrite('inc.eta.cloud.nc4','v',v_b);
% ncwrite('inc.eta.cloud.nc4','tv',t_b);
% ncwrite('inc.eta.cloud.nc4','sphu',q_b);
% ncwrite('inc.eta.cloud.nc4','delp',p_b);
% ncwrite('inc.eta.cloud.nc4','qitot',qi_b);
% ncwrite('inc.eta.cloud.nc4','qltot',ql_b);
% ncwrite('inc.eta.cloud.nc4','ozone',o3_b);
% cd /home/drholdaw/Lin_Moist_Physics/DAS_plotting/Cloud' Fraction'/



%Choose model level to plot.
plot_level = 54;

match_scale = 0;

% fprintf('\nu:   %g \n\n', max(abs(u_b(:)))/max(abs(u_a(:))))
% fprintf('\nv:   %g \n\n', max(abs(v_b(:)))/max(abs(v_a(:))))
% fprintf('\nt:   %g \n\n', max(abs(t_b(:)))/max(abs(t_a(:))))
% fprintf('\nq:   %g \n\n', max(abs(q_b(:)))/max(abs(q_a(:))))
% fprintf('\np:   %g \n\n', max(abs(p_b(:)))/max(abs(p_a(:))))
% fprintf('\nqi:   %g \n\n', max(abs(qi_b(:)))/max(abs(qi_a(:))))
% fprintf('\nql:   %g \n\n', max(abs(ql_b(:)))/max(abs(ql_a(:))))
% fprintf('\no3:   %g \n\n', max(abs(o3_b(:)))/max(abs(o3_a(:))))

[maxival location] = max(abs(u_a(:)-u_b(:)));
[rmax,cmax,lmax] = ind2sub(size(u_a),location);
fprintf('Max u difference is %f \n',maxival)
fprintf('Located at (%d,%d,%d) \n',rmax,cmax,lmax)
fprintf('u_a = %f, u_b = %f \n\n',u_a(rmax,cmax,lmax),u_b(rmax,cmax,lmax))

[maxival location] = max(abs(v_a(:)-v_b(:)));
[rmax,cmax,lmax] = ind2sub(size(v_a),location);
fprintf('Max v difference is %f \n',maxival)
fprintf('Located at (%d,%d,%d) \n',rmax,cmax,lmax)
fprintf('v_a = %f, v_b = %f \n\n',v_a(rmax,cmax,lmax),v_b(rmax,cmax,lmax))

[maxival location] = max(abs(t_a(:)-t_b(:)));
[rmax,cmax,lmax] = ind2sub(size(t_a),location);
fprintf('Max t difference is %f \n',maxival)
fprintf('Located at (%d,%d,%d) \n',rmax,cmax,lmax)
fprintf('t_a = %f, t_b = %f \n\n',t_a(rmax,cmax,lmax),t_b(rmax,cmax,lmax))

[maxival location] = max(abs(q_a(:)-q_b(:)));
[rmax,cmax,lmax] = ind2sub(size(q_a),location);
fprintf('Max q difference is %f \n',maxival)
fprintf('Located at (%d,%d,%d) \n',rmax,cmax,lmax)
fprintf('q_a = %f, q_b = %f \n\n',q_a(rmax,cmax,lmax),q_b(rmax,cmax,lmax))

[maxival location] = max(abs(p_a(:)-p_b(:)));
[rmax,cmax,lmax] = ind2sub(size(p_a),location);
fprintf('Max p difference is %f \n',maxival)
fprintf('Located at (%d,%d,%d) \n',rmax,cmax,lmax)
fprintf('p_a = %f, p_b = %f \n\n',p_a(rmax,cmax,lmax),p_b(rmax,cmax,lmax))

[maxival location] = max(abs(qi_a(:)-qi_b(:)));
[rmax,cmax,lmax] = ind2sub(size(qi_a),location);
fprintf('Max qi difference is %f \n',maxival)
fprintf('Located at (%d,%d,%d) \n',rmax,cmax,lmax)
fprintf('qi_a = %f, qi_b = %f \n\n',qi_a(rmax,cmax,lmax),qi_b(rmax,cmax,lmax))

[maxival location] = max(abs(ql_a(:)-ql_b(:)));
[rmax,cmax,lmax] = ind2sub(size(ql_a),location);
fprintf('Max ql difference is %f \n',maxival)
fprintf('Located at (%d,%d,%d) \n',rmax,cmax,lmax)
fprintf('ql_a = %f, ql_b = %f \n\n',ql_a(rmax,cmax,lmax),ql_b(rmax,cmax,lmax))

u_a_lev = u_a(:,:,plot_level);
v_a_lev = v_a(:,:,plot_level);
t_a_lev = t_a(:,:,plot_level);
q_a_lev = q_a(:,:,plot_level);
p_a_lev = p_a(:,:,plot_level);
qi_a_lev = qi_a(:,:,plot_level);
ql_a_lev = ql_a(:,:,plot_level);
o3_a_lev = o3_a(:,:,plot_level);

u_b_lev = u_b(:,:,plot_level);
v_b_lev = v_b(:,:,plot_level);
t_b_lev = t_b(:,:,plot_level);
q_b_lev = q_b(:,:,plot_level);
p_b_lev = p_b(:,:,plot_level);
qi_b_lev = qi_b(:,:,plot_level);
ql_b_lev = ql_b(:,:,plot_level);
o3_b_lev = o3_b(:,:,plot_level);

%Zerodiff
all_a = [u_a v_a t_a q_a p_a qi_a ql_a o3_a];
all_b = [u_b v_b t_b q_b p_b qi_b ql_b o3_b];
% all_a = [qi_a];
% all_b = [qi_b];

all_diff = all_a - all_b;
zd = max(abs(all_diff(:)));

if zd == 0 
    fprintf('Zero Diff Passed \n\n');
else
    fprintf('Not Zero Diff, biggest difference is %f \n\n',zd);
end

%Get Maximums for Levels
u_lev = [u_a_lev u_b_lev];
maxu = max(abs(u_lev(:)));
if min(u_lev(:)) >= 0
    minu = 0;
else
    minu = -maxu;
end
v_lev = [v_a_lev v_b_lev];
maxv = max(abs(v_lev(:)));
if min(v_lev(:)) >= 0
    minv = min(v_lev(:));
else
    minv = -maxv;
end
t_lev = [t_a_lev t_b_lev];
maxt = max(abs(t_lev(:)));
if min(t_lev(:)) >= 0
    mint = min(t_lev(:));
else
    mint = -maxt;
end
q_lev = [q_a_lev q_b_lev];
maxq = max(abs(q_lev(:)));
if min(q_lev(:)) >= 0
    minq = min(q_lev(:));
else
    minq = -maxq;
end
p_lev = [p_a_lev p_b_lev];
maxp = max(abs(p_lev(:)));
if min(p_lev(:)) >= 0
    minp = min(p_lev(:));
else
    minp = -maxp;
end
qi_lev = [qi_a_lev qi_b_lev];
maxqi = max(abs(qi_lev(:)));
if min(qi_lev(:)) >= 0
    minqi = min(qi_lev(:));
else
    minqi = -maxqi;
end
ql_lev = [ql_a_lev ql_b_lev];
maxql = max(abs(ql_lev(:)));
if min(ql_lev(:)) >= 0
    minql = min(ql_lev(:));
else
    minql = -maxql;
end
o3_lev = [o3_a_lev o3_b_lev];
maxo3 = max(abs(o3_lev(:)));
if min(o3_lev(:)) >= 0
    mino3 = min(o3_lev(:));
else
    mino3 = -maxo3;
end

if match_scale == 0
    maxu = max(abs(u_a_lev(:)));
    if min(u_a_lev(:)) >= 0
        minu = min(u_lev(:));
    else
        minu = -maxu;
    end
    maxv = max(abs(v_a_lev(:)));
    if min(v_a_lev(:)) >= 0
        minv = min(v_lev(:));
    else
        minv = -maxv;
    end
    maxt = max(abs(t_a_lev(:)));
    if min(t_a_lev(:)) >= 0
        mint = min(t_lev(:));
    else
        mint = -maxt;
    end
    maxq = max(abs(q_a_lev(:)));
    if min(q_a_lev(:)) >= 0
        minq = min(q_lev(:));
    else
        minq = -maxq;
    end
    maxp = max(abs(p_a_lev(:)));
    if min(p_a_lev(:)) >= 0
        minp = min(p_lev(:));
    else
        minp = -maxp;
    end
    maxqi = max(abs(qi_a_lev(:)));
    if min(qi_a_lev(:)) >= 0
        minqi = min(qi_lev(:));
    else
        minqi = -maxqi;
    end
    maxql = max(abs(ql_a_lev(:)));
    if min(ql_a_lev(:)) >= 0
        minql = min(ql_lev(:));
    else
        minql = -maxql;
    end
    maxo3 = max(abs(o3_a_lev(:)));
    if min(o3_a_lev(:)) >= 0
        mino3 = min(o3_lev(:));
    else
        mino3 = -maxo3;
    end
end
     

figure
set(gcf,'position',[349 30 930 889])

subplot(4,2,1)
contourf(lon,lat,((u_a_lev))','LineStyle','none')
caxis([minu maxu])
title('u')
colorbar

subplot(4,2,2)
contourf(lon,lat,((v_a_lev))','LineStyle','none')
caxis([minv maxv])
title('v')
colorbar

subplot(4,2,3)
contourf(lon,lat,((t_a_lev))','LineStyle','none')
caxis([mint maxt])
title('\theta')
colorbar

subplot(4,2,4)
contourf(lon,lat,((q_a_lev))','LineStyle','none')
caxis([minq maxq])
title('q')
colorbar

subplot(4,2,5)
contourf(lon,lat,((p_a_lev))','LineStyle','none')
caxis([minp maxp])
title('dp')
colorbar

subplot(4,2,6)
contourf(lon,lat,((qi_a_lev))','LineStyle','none')
caxis([minqi maxqi])
title('Qi')
colorbar

subplot(4,2,7)
contourf(lon,lat,((ql_a_lev))','LineStyle','none')
caxis([minql maxql])
title('Ql')
colorbar

subplot(4,2,8)
contourf(lon,lat,((o3_a_lev))','LineStyle','none')
caxis([mino3 maxo3])
title('O3')
colorbar


if match_scale == 0
    maxu = max(abs(u_b_lev(:)));
    if min(u_b_lev(:)) >= 0
        minu = min(u_lev(:));
    else
        minu = -maxu;
    end
    maxv = max(abs(v_b_lev(:)));
    if min(v_b_lev(:)) >= 0
        minv = min(v_lev(:));
    else
        minv = -maxv;
    end
    maxt = max(abs(t_b_lev(:)));
    if min(t_b_lev(:)) >= 0
        mint = min(t_lev(:));
    else
        mint = -maxt;
    end
    maxq = max(abs(q_b_lev(:)));
    if min(q_b_lev(:)) >= 0
        minq = min(q_lev(:));
    else
        minq = -maxq;
    end
    maxp = max(abs(p_b_lev(:)));
    if min(p_b_lev(:)) >= 0
        minp = min(p_lev(:));
    else
        minp = -maxp;
    end
    maxqi = max(abs(qi_b_lev(:)));
    if min(qi_b_lev(:)) >= 0
        minqi = min(qi_lev(:));
    else
        minqi = -maxqi;
    end
    maxql = max(abs(ql_b_lev(:)));
    if min(ql_b_lev(:)) >= 0
        minql = min(ql_lev(:));
    else
        minql = -maxql;
    end
    maxo3 = max(abs(o3_b_lev(:)));
    if min(o3_b_lev(:)) >= 0
        mino3 = min(o3_lev(:));
    else
        mino3 = -maxo3;
    end
end


figure
set(gcf,'position',[349 30 930 889])

subplot(4,2,1)
contourf(lon,lat,((u_b_lev))','LineStyle','none')
caxis([minu maxu])
title('u')
colorbar

subplot(4,2,2)
contourf(lon,lat,((v_b_lev))','LineStyle','none')
caxis([minv maxv])
title('v')
colorbar

subplot(4,2,3)
contourf(lon,lat,((t_b_lev))','LineStyle','none')
caxis([mint maxt])
title('\theta')
colorbar

subplot(4,2,4)
contourf(lon,lat,((q_b_lev))','LineStyle','none')
caxis([minq maxq])
title('q')
colorbar

subplot(4,2,5)
contourf(lon,lat,((p_b_lev))','LineStyle','none')
caxis([minp maxp])
title('dp')
colorbar

subplot(4,2,6)
contourf(lon,lat,((qi_b_lev))','LineStyle','none')
caxis([minqi maxqi])
title('Qi')
colorbar

subplot(4,2,7)
contourf(lon,lat,((ql_b_lev))','LineStyle','none')
caxis([minql maxql])
title('Ql')
colorbar

subplot(4,2,8)
contourf(lon,lat,((o3_b_lev))','LineStyle','none')
caxis([mino3 maxo3])
title('O3')
colorbar















