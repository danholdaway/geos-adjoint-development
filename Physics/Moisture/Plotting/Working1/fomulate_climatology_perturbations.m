close all
clear
clc
mydir = pwd;

%Load MATLAB topography
load coast
coast_lat = lat; clear lat
coast_lon = long; clear long

dir1 = '/discover/nobackup/drholdaw/btmp.25956/prog/prog_replay';
dir2 = '/discover/nobackup/drholdaw/btmp.25956/prog/prog_free';

file1 = 'v000_C180.prog.eta.20140201_03z.nc4';
file2 = 'v000_C180.prog.eta.20140201_03z.nc4';

cd(dir1)
file = file1;

lon = ncread(file,'lon'); im = length(lon);
lat = ncread(file,'lat'); jm = length(lat);
lev = ncread(file,'lev'); lm = length(lev);

t_1 = ncread(file,'tv');
q_1 = ncread(file,'sphu');
qi_1 = ncread(file,'qitot');
ql_1 = ncread(file,'qltot');
cfls_1 = ncread(file,'cfls');
cfcn_1 = ncread(file,'cfcn');

cd(dir2)
file = file2;

t_2 = ncread(file,'tv');
q_2 = ncread(file,'sphu');
qi_2 = ncread(file,'qitot');
ql_2 = ncread(file,'qltot');
cfls_2 = ncread(file,'cfls');
cfcn_2 = ncread(file,'cfcn');

cd(mydir)



t_nlm = t_2 - t_1;
q_nlm = q_2 - q_1;
qi_nlm = qi_2 - qi_1;
ql_nlm = ql_2 - ql_1;
cfls_nlm = cfls_2 - cfls_1;
cfcn_nlm = cfcn_2 - cfcn_1;


mean_t = zeros(1,lm);
mean_q = zeros(1,lm);
mean_qi = zeros(1,lm);
mean_ql = zeros(1,lm);
mean_cfls = zeros(1,lm);
mean_cfcn = zeros(1,lm);

for k = 1:lm
    
    mean_t(k) = mean(mean(t_nlm(:,:,k)));
    mean_q(k) = mean(mean(q_nlm(:,:,k)));
    mean_qi(k) = mean(mean(qi_nlm(:,:,k)));
    mean_ql(k) = mean(mean(ql_nlm(:,:,k)));
    mean_cfls(k) = mean(mean(cfls_nlm(:,:,k)));
    mean_cfcn(k) = mean(mean(cfcn_nlm(:,:,k)));
    
end

max_t = zeros(1,lm);
max_q = zeros(1,lm);
max_qi = zeros(1,lm);
max_ql = zeros(1,lm);
max_cfls = zeros(1,lm);
max_cfcn = zeros(1,lm);

for k = 1:lm
    
    max_t(k) = max(max(t_nlm(:,:,k)));
    max_q(k) = max(max(q_nlm(:,:,k)));
    max_qi(k) = max(max(qi_nlm(:,:,k)));
    max_ql(k) = max(max(ql_nlm(:,:,k)));
    max_cfls(k) = max(max(cfls_nlm(:,:,k)));
    max_cfcn(k) = max(max(cfcn_nlm(:,:,k)));
    
end

std_t = zeros(1,lm);
std_q = zeros(1,lm);
std_qi = zeros(1,lm);
std_ql = zeros(1,lm);
std_cfls = zeros(1,lm);
std_cfcn = zeros(1,lm);

for k = 1:lm
    
    std_t(k) = std(std(t_nlm(:,:,k)));
    std_q(k) = std(std(q_nlm(:,:,k)));
    std_qi(k) = std(std(qi_nlm(:,:,k)));
    std_ql(k) = std(std(ql_nlm(:,:,k)));
    std_cfls(k) = std(std(cfls_nlm(:,:,k)));
    std_cfcn(k) = std(std(cfcn_nlm(:,:,k)));
    
end

for k = 1:lm
    
    if mean_t(k) < 1e-8
        mean_t(k) = 0.0;
    end
    if mean_q(k) < 1e-8
        mean_q(k) = 0.0;
    end
    if mean_qi(k) < 1e-8
        mean_qi(k) = 0.0;
    end
    if mean_ql(k) < 1e-8
        mean_ql(k) = 0.0;
    end
    if mean_cfls(k) < 1e-8
        mean_cfls(k) = 0.0;
    end
    if mean_cfcn(k) < 1e-8
        mean_cfcn(k) = 0.0;
    end

end

for k = 1:lm
    
    if max_t(k) < 1e-8
        max_t(k) = 0.0;
    end
    if max_q(k) < 1e-8
        max_q(k) = 0.0;
    end
    if max_qi(k) < 1e-8
        max_qi(k) = 0.0;
    end
    if max_ql(k) < 1e-8
        max_ql(k) = 0.0;
    end
    if max_cfls(k) < 1e-8
        max_cfls(k) = 0.0;
    end
    if max_cfcn(k) < 1e-8
        max_cfcn(k) = 0.0;
    end

end


for k = 1:lm
    
    if std_t(k) < 1e-8
        std_t(k) = 0.0;
    end
    if std_q(k) < 1e-8
        std_q(k) = 0.0;
    end
    if std_qi(k) < 1e-8
        std_qi(k) = 0.0;
    end
    if std_ql(k) < 1e-8
        std_ql(k) = 0.0;
    end
    if std_cfls(k) < 1e-8
        std_cfls(k) = 0.0;
    end
    if std_cfcn(k) < 1e-8
        std_cfcn(k) = 0.0;
    end

end
    
fprintf('mean_T_pert = (/ ')
for k = 1:lm-1
    fprintf('%10.7e, ', abs(mean_t(k)))
end
fprintf('%g', abs(mean_t(lm)))
fprintf(' /) \n')

fprintf('max_T_pert = (/ ')
for k = 1:lm-1
    fprintf('%10.7e, ', abs(max_t(k)))
end
fprintf('%g', abs(max_t(lm)))
fprintf(' /) \n')

fprintf('std_T_pert = (/ ')
for k = 1:lm-1
    fprintf('%10.7e, ', std_t(k))
end
fprintf('%g', std_t(lm))
fprintf(' /) \n')



fprintf('mean_q_pert = (/ ')
for k = 1:lm-1
    fprintf('%10.7e, ', abs(mean_q(k)))
end
fprintf('%g', abs(mean_q(lm)))
fprintf(' /) \n')

fprintf('max_q_pert = (/ ')
for k = 1:lm-1
    fprintf('%10.7e, ', abs(max_q(k)))
end
fprintf('%g', abs(max_q(lm)))
fprintf(' /) \n')

fprintf('std_q_pert = (/ ')
for k = 1:lm-1
    fprintf('%10.7e, ', std_q(k))
end
fprintf('%g', std_q(lm))
fprintf(' /) \n')


fprintf('mean_qi_pert = (/ ')
for k = 1:lm-1
    fprintf('%10.7e, ', abs(mean_qi(k)))
end
fprintf('%g', abs(mean_qi(lm)))
fprintf(' /) \n')

fprintf('max_qi_pert = (/ ')
for k = 1:lm-1
    fprintf('%10.7e, ', abs(max_qi(k)))
end
fprintf('%g', abs(max_qi(lm)))
fprintf(' /) \n')

fprintf('std_qi_pert = (/ ')
for k = 1:lm-1
    fprintf('%10.7e, ', std_qi(k))
end
fprintf('%g', std_qi(lm))
fprintf(' /) \n')


fprintf('mean_ql_pert = (/ ')
for k = 1:lm-1
    fprintf('%10.7e, ', abs(mean_ql(k)))
end
fprintf('%g', abs(mean_ql(lm)))
fprintf(' /) \n')

fprintf('max_ql_pert = (/ ')
for k = 1:lm-1
    fprintf('%10.7e, ', abs(max_ql(k)))
end
fprintf('%g', abs(max_ql(lm)))
fprintf(' /) \n')

fprintf('std_ql_pert = (/ ')
for k = 1:lm-1
    fprintf('%10.7e, ', std_ql(k))
end
fprintf('%g', std_ql(lm))
fprintf(' /) \n')


fprintf('mean_cf_pert = (/ ')
for k = 1:lm-1
    fprintf('%10.7e, ', abs(mean_cfls(k)))
end
fprintf('%g', abs(mean_cfls(lm)))
fprintf(' /) \n')

fprintf('max_cf_pert = (/ ')
for k = 1:lm-1
    fprintf('%10.7e, ', abs(max_cfls(k)))
end
fprintf('%g', abs(max_cfls(lm)))
fprintf(' /) \n')

fprintf('std_cf_pert = (/ ')
for k = 1:lm-1
    fprintf('%10.7e, ', std_cfls(k))
end
fprintf('%g', std_cfls(lm))
fprintf(' /) \n')


fprintf('mean_af_pert = (/ ')
for k = 1:lm-1
    fprintf('%10.7e, ', abs(mean_cfcn(k)))
end
fprintf('%g', abs(mean_cfcn(lm)))
fprintf(' /) \n')

fprintf('max_af_pert = (/ ')
for k = 1:lm-1
    fprintf('%10.7e, ', abs(max_cfcn(k)))
end
fprintf('%g', abs(max_cfcn(lm)))
fprintf(' /) \n')

fprintf('std_af_pert = (/ ')
for k = 1:lm-1
    fprintf('%10.7e, ', std_cfcn(k))
end
fprintf('%g', std_cfcn(lm))
fprintf(' /) \n')












y = 1:lm;





figure
subplot(1,2,1)
plot(abs(mean_t),y)
hold on
plot(abs(mean_t)+std_t,y,'r')
plot(abs(mean_t)-std_t,y,'g')
set(gca,'YDir','reverse')
ylim([1 72])
subplot(1,2,2)
plot(abs(max_t),y)
set(gca,'YDir','reverse')
ylim([1 72])

figure
subplot(1,2,1)
plot(abs(mean_q),y)
set(gca,'YDir','reverse')
ylim([1 72])
subplot(1,2,2)
plot(abs(max_q),y)
set(gca,'YDir','reverse')
ylim([1 72])

figure
subplot(1,2,1)
plot(abs(mean_qi),y)
set(gca,'YDir','reverse')
ylim([1 72])
subplot(1,2,2)
plot(abs(max_qi),y)
set(gca,'YDir','reverse')
ylim([1 72])

figure
subplot(1,2,1)
plot(abs(mean_ql),y)
set(gca,'YDir','reverse')
ylim([1 72])
subplot(1,2,2)
plot(abs(max_ql),y)
set(gca,'YDir','reverse')
ylim([1 72])

figure
subplot(1,2,1)
plot(abs(mean_cfls),y)
set(gca,'YDir','reverse')
ylim([1 72])
subplot(1,2,2)
plot(abs(max_cfls),y)
set(gca,'YDir','reverse')
ylim([1 72])

figure
subplot(1,2,1)
plot(abs(mean_cfcn),y)
set(gca,'YDir','reverse')
ylim([1 72])
subplot(1,2,2)
plot(abs(max_cfcn),y)
set(gca,'YDir','reverse')
ylim([1 72])
