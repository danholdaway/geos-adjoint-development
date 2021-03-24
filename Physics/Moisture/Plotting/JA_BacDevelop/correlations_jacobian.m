close all
clear
clc
mydir = pwd;

cd /discover/nobackup/drholdaw/ExperimentData/BacmeisterPaper/Correlations/

jac_lscloud_t1
var1 = td;
J11 = J1;
J12 = J2;
J13 = J5;

clear td J1 J2 J3 J4 J5 J6 J7 J8

jac_lscloud_q1
var2 = qd;
J21 = J1;
J22 = J2;
J23 = J5;

clear qd J1 J2 J3 J4 J5 J6 J7 J8

jac_lscloud_qlls1
var3 = ql_lsd;
J31 = J1;
J32 = J2;
J33 = J5;

clear ql_lsd J1 J2 J3 J4 J5 J6 J7 J8

jac_lscloud_cfls1
var4 = cf_lsd;
J41 = J1;
J42 = J2;
J43 = J5;

clear cf_lsd J1 J2 J3 J4 J5 J6 J7 J8

cd(mydir)

fontsize = 11;


figure
subplot(4,3,1)
scatter(var1,J11,'kx')
box on
xlabel('\partial T\prime (K)','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('\partial F/ \partial T','FontSize',fontsize,'FontName','TimesNewRoman')
title('(a)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([-0.9 0.3])
ylim([0.65 1.0])

subplot(4,3,2)
scatter(var1,J12,'kx')
box on
xlabel('\partial T\prime (K)','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('\partial F/ \partial q','FontSize',fontsize,'FontName','TimesNewRoman')
title('(b)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([-0.9 0.3])
ylim([0 3000])


subplot(4,3,3)
scatter(var1,J13,'kx')
box on
xlabel('\partial T\prime (K)','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('\partial F/ \partial q_{l,LS}','FontSize',fontsize,'FontName','TimesNewRoman')
title('(c)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([-0.9 0.3])
ylim([-3000 1])




subplot(4,3,4)
scatter(var2,J21,'kx')
box on
xlabel('\partial q\prime (kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('\partial G/ \partial T','FontSize',fontsize,'FontName','TimesNewRoman')
title('(d)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([-1e-4 3.5e-4])
ylim([0 1.5e-4])

subplot(4,3,5)
scatter(var2,J22,'kx')
box on
xlabel('\partial q\prime (kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('\partial G/ \partial q','FontSize',fontsize,'FontName','TimesNewRoman')
title('(e)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([-1e-4 3.5e-4])
ylim([0 1])

subplot(4,3,6)
scatter(var2,J23,'kx')
box on
xlabel('\partial q\prime (kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('\partial G/ \partial q_{l,LS}','FontSize',fontsize,'FontName','TimesNewRoman')
title('(f)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([-1e-4 3.5e-4])
ylim([0 1])




subplot(4,3,7)
scatter(var3,J31,'kx')
box on
xlabel('\partial q_{l,LS}\prime (kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('\partial H/ \partial T','FontSize',fontsize,'FontName','TimesNewRoman')
title('(g)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([-3.5e-4 1e-4])
ylim([-1.5e-4 0])

subplot(4,3,8)
scatter(var3,J32,'kx')
box on
xlabel('\partial q_{l,LS}\prime (kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('\partial H/ \partial q','FontSize',fontsize,'FontName','TimesNewRoman')
title('(h)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([-3.5e-4 1e-4])
ylim([0 0.85])

subplot(4,3,9)
scatter(var3,J33,'kx')
box on
xlabel('\partial q_{l,LS}\prime (kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('\partial H/ \partial q_{l,LS}','FontSize',fontsize,'FontName','TimesNewRoman')
title('(i)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([-3.5e-4 1e-4])
ylim([0 1])




subplot(4,3,10)
scatter(var4,J41,'kx')
box on
xlabel('\partial C_{LS}\prime','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('\partial J/ \partial T','FontSize',fontsize,'FontName','TimesNewRoman')
title('(j)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([-5.5 1.8])
ylim([-5 0])

subplot(4,3,11)
scatter(var4,J42,'kx')
box on
xlabel('\partial C_{LS}\prime','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('\partial J/ \partial q','FontSize',fontsize,'FontName','TimesNewRoman')
title('(k)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([-5.5 1.8])
ylim([0 7e5])

subplot(4,3,12)
scatter(var4,J43,'kx')
box on
xlabel('\partial C_{LS}\prime','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('\partial J/ \partial q_{l,LS}','FontSize',fontsize,'FontName','TimesNewRoman')
title('(l)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([-5.5 1.8])
ylim([0 7e5])

