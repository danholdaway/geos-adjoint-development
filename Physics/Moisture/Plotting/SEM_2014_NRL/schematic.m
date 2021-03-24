close all
clear
clc

%Choose model level to plot.
plot_level = 63;

fontsize = 13;

line_wid_cont = 1.0;
line_wid_det = 0.6;

grey = 0.75;

Ri = -33:0.001:0;
fm = zeros(1,length(Ri));
fh = zeros(1,length(Ri));

fm_new = zeros(1,length(Ri));
fh_new = zeros(1,length(Ri));

dfmdRi = zeros(1,length(Ri)-1);
dfhdRi = zeros(1,length(Ri)-1);
dfm_newdRi = zeros(1,length(Ri)-1);
dfh_newdRi = zeros(1,length(Ri)-1);

b = 5;
c = 5;
d = 5;

z = 100;
z0 = 0.01;
kappa = 0.41;
gamma = 0.2;
L = 4;
k = 0.1;

lambdam = 150;
lambdah = lambdam*sqrt(1.5*d);

lm = kappa*(z+z0)/(1+k*(z+z0)/lambdam)*(gamma + (1-gamma)/(1 + (z+z0)^2/L^2));
lh = kappa*(z+z0)/(1+k*(z+z0)/lambdah)*(gamma + (1-gamma)/(1 + (z+z0)^2/L^2));

for i = 1:length(Ri)
    
    if Ri(i) < 0 
    
        psi = Ri(i) / (1 + 3 * b * c * sqrt(-Ri(i)));
                
        fm(i) =     (1 - 2*b*psi);
        
        if Ri(i) < -6 
            fh(i) = 0.05+(1 - 3*b*psi);
        else
            fh(i) = (1 - 2*b*psi);
        end
        
    end
    
      
    if Ri(i) == 0 
        
        psi = sqrt(1+d*Ri(i));
        
        fm(i) = 1;
        fh(i) = 1;
        
    end
    
end

Rih = 0.5*(Ri(1:end-1) + Ri(2:end));

figure
set(gcf,'position',[719   553   560   366])

plot(Ri(1:27000),fh(1:27000),'r','LineWidth',1.5)
% asd
hold on
plot(Ri,fm,'b','LineWidth',1.5)

% plot(Ri,fm_new,'b--')
hold on
% plot(Ri,fh_new,'r--')
legend('Perturbed Run','Free Run','location','NorthWest')
legend boxoff
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlabel('Time (UTC)','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Nonlinear trajectory','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'XDir','reverse')
xlim([-33 -3])
set(gca,'Ytick',[])
set(gca,'Xtick',[-33 -30 -27 -24 -21 -18 -15 -12 -9 -6 -3])
set(gca,'XtickLabel',{'0300','0000','2100','1800','1500','1200','0900','0600','0300','0000','2100'})

Hl = 6.0;
Hw = 6.0;

annotation('doublearrow',[0.190 0.190],[0.21 0.38],'LineWidth',1.5...
           ,'Head1Style','cback1','Head2Style','cback1'...
           ,'Head1Length',Hl,'Head1Width',Hw,'Head2Length',Hl,'Head2Width',Hw)
annotation('doublearrow',[0.285 0.285],[0.27 0.47],'LineWidth',1.5...
           ,'Head1Style','cback1','Head2Style','cback1'...
           ,'Head1Length',Hl,'Head1Width',Hw,'Head2Length',Hl,'Head2Width',Hw)
           
annotation('textbox',[0.13 0.23 0.1 0.1],'String','\delta x','LineStyle','none'...
           ,'FontSize',fontsize,'FontName','TimesNewRoman')
       
annotation('textbox',[0.28 0.31 0.1 0.1],'String','TLM restart','LineStyle','none'...
           ,'FontSize',fontsize,'FontName','TimesNewRoman')


%FOR AMS MEETING SLIDES
       
% figure
% set(gcf,'position',[719   553   560   366])
% 
% plot(Ri,fh,'r','LineWidth',2.0)
% hold on
% plot(Ri,fm,'b','LineWidth',2.0)
% 
% % plot(Ri,fm_new,'b--')
% hold on
% % plot(Ri,fh_new,'r--')
% legend('Perturbed Run','Free Run','location','NorthWest')
% set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
% xlabel('Time (hours)','FontSize',fontsize,'FontName','TimesNewRoman')
% ylabel('Nonlinear trajectory','FontSize',fontsize,'FontName','TimesNewRoman')
% set(gca,'XDir','reverse')
% xlim([-33 -3])
% set(gca,'Ytick',[])
% set(gca,'Xtick',[-33 -30 -27 -24 -21 -18 -15 -12 -9 -6 -3])
% set(gca,'XtickLabel',[03 00 21 18 15 12 09 06 03 00 21])
% 
% Hl = 6.0;
% Hw = 6.0;
% 
% annotation('doublearrow',[0.190 0.190],[0.21 0.38],'LineWidth',1.75...
%            ,'Head1Style','cback1','Head2Style','cback1'...
%            ,'Head1Length',Hl,'Head1Width',Hw,'Head2Length',Hl,'Head2Width',Hw)
% annotation('doublearrow',[0.285 0.285],[0.27 0.47],'LineWidth',1.75...
%            ,'Head1Style','cback1','Head2Style','cback1'...
%            ,'Head1Length',Hl,'Head1Width',Hw,'Head2Length',Hl,'Head2Width',Hw)
%            
% annotation('textbox',[0.13 0.23 0.1 0.1],'String','\delta x','LineStyle','none'...
%            ,'FontSize',fontsize,'FontName','TimesNewRoman')
%        
% annotation('textbox',[0.28 0.31 0.1 0.1],'String','TLM restart','LineStyle','none'...
%            ,'FontSize',fontsize,'FontName','TimesNewRoman')

       
       