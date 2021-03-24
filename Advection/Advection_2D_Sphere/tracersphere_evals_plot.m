close all
clear
clc
mydir = pwd;

cd /discover/nobackup/drholdaw/Models/Advection_2D_Sphere/Exp/

maxevallw

figure

subplot(2,2,1)
plot(maxevals)
title('Lax-Wendroff')
ylabel('Max Eigenvalue')

maxevallwlim

subplot(2,2,2)
plot(maxevals)
title('Limited Lax-Wendroff')

maxevalut

subplot(2,2,3)
plot(maxevals)
title('Third order')
ylabel('Max Eigenvalue')
xlabel('N')

maxevalutlim

subplot(2,2,4)
plot(maxevals)
title('Limited third order')
xlabel('N')



evals
figure
scatter(e_r,e_i)

clear
revector

cd(mydir); pwd



