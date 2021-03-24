close all
clear
clc
mydir = pwd;

cd /home/drholdaw/Advection/Advection_1D_Offline/
jac

x = 1:64;

J1 = J1;
J2 = J3;

[U1,eig1] = eig(J1);
eig1 = diag(eig1);
[U2,eig2] = eig(J2);
eig2 = diag(eig2);

%Sort eigenvalues
[reig1, ind] = sort(real(eig1),'descend');
ieig1 = imag(eig1(ind));
U1 = U1(:,ind);

[reig2, ind] = sort(real(eig2),'descend');
ieig2 = imag(eig2(ind));
U2 = U2(:,ind);


xx = zeros(1,1001);
xxx = 0:(1/63):1;

i = 1;

fontsize = 14;
lin_wid = 1.5;



%PPM SCHEME
figure
set(gcf,'position',[228 553 1051 366])

subplot(2,2,[1,3])
scatter(reig1,ieig1,'LineWidth',lin_wid)
hold on
scatter(reig1(i),ieig1(i),150,'k.','LineWidth',lin_wid)
scatter(reig1(i+1),ieig1(i+1),150,'k.','LineWidth',lin_wid)
plot(xx,-500:500,'k','LineWidth',1.2)
hold off

box on
xlabel('Re(\lambda)','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Im(\lambda)','FontSize',fontsize,'FontName','TimesNewRoman')
title('Eigenvalues','FontSize',fontsize,'FontName','TimesNewRoman')
ylim([-150 150])
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(2,2,2)
plot(xxx,real(U1(:,i)),'b','LineWidth',lin_wid)
hold on
plot(xxx,imag(U1(:,i)),'r','LineWidth',lin_wid)
hold off

box on
xlabel('x','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('q\prime','FontSize',fontsize,'FontName','TimesNewRoman')
title('Eigenvector 1','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

subplot(2,2,4)
plot(xxx,real(U1(:,i+1)),'b','LineWidth',lin_wid)
hold on
plot(xxx,imag(U1(:,i+1)),'r','LineWidth',lin_wid)
hold off

box on
xlabel('x','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('q\prime','FontSize',fontsize,'FontName','TimesNewRoman')
title('Eigenvector 2','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
legend('Real','Imaginary')



















cd(mydir)