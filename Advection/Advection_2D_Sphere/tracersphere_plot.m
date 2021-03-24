close all
clear
clc
mydir = pwd;

cd /discover/nobackup/drholdaw/Models/Advection_2D_Sphere/Exp/

tracersphere;

u_nl = qreplut-qfreeut;


figure
contourf(qfreeut','LineStyle','none')
colorbar

figure
contourf(qreplut','LineStyle','none')
colorbar

figure
contourf(u_nl','LineStyle','none')
colorbar

figure
contourf(qtlmut','LineStyle','none')
colorbar


figure
contourf(qtlmut'-u_nl','LineStyle','none')
colorbar

asd

figure
subplot(3,1,1)
contourf(qtlmlw','LineStyle','none')
colorbar
subplot(3,1,2)
contourf(qrepllw'-qfreelw','LineStyle','none')
colorbar
subplot(3,1,3)
contourf(qtlmlw'-(qrepllw'-qfreelw'),'LineStyle','none')
colorbar

figure
subplot(3,1,1)
contourf(qtlmlwlim','LineStyle','none')
colorbar
subplot(3,1,2)
contourf(qrepllwlim'-qfreelwlim','LineStyle','none')
colorbar
subplot(3,1,3)
contourf(qtlmlwlim'-(qrepllwlim'-qfreelwlim'),'LineStyle','none')
colorbar



clear

tracersphere1;

figure
contourf(qfree','LineStyle','none')
colorbar

figure
contourf(qfreelw','LineStyle','none')
colorbar

figure
contourf(qfreelwlim','LineStyle','none')
colorbar

figure
subplot(3,1,1)
contourf(qtlmlw','LineStyle','none')
colorbar
subplot(3,1,2)
contourf(qrepllw'-qfreelw','LineStyle','none')
colorbar
subplot(3,1,3)
contourf(qtlmlw'-(qrepllw'-qfreelw'),'LineStyle','none')
colorbar

figure
subplot(3,1,1)
contourf(qtlmlwlim','LineStyle','none')
colorbar
subplot(3,1,2)
contourf(qrepllwlim'-qfreelwlim','LineStyle','none')
colorbar
subplot(3,1,3)
contourf(qtlmlwlim'-(qrepllwlim'-qfreelwlim'),'LineStyle','none')
colorbar


cd(mydir); clear mydir