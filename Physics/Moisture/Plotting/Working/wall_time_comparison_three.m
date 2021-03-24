close all
clear
clc

tlm_6_dry = 243.8051 - 69.7098;
tlm_6_moi_nofilt = 243.8051 - 16.5844;
tlm_6_moi = 243.8051;

tlm_percent = 100*((tlm_6_moi_nofilt/tlm_6_dry)-1)


%TLM: DRY, MOIST_NOFILT, MOIST
walltime(1,1) =   tlm_6_dry; 
walltime(1,2) =   tlm_6_moi_nofilt; 
walltime(1,3) =   tlm_6_moi; 



times = [6 12 18];


bar(walltime./60)
colormap summer

% xlim([0 3.6])

fontsize = 14;

title('24 Hour Adjoint Computational Time','FontSize',fontsize,'FontName','TimesNewRoman')

% xlabel('Window Length','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('Time (Minutes)','FontSize',fontsize,'FontName','TimesNewRoman')

% legend('Dry','Moist','With Filtering','Location','NorthWest')

set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')

saveas(gcf,'processor_dist.eps','psc2')

set(gca,'XtickLabel',['Dry';'Moi';'Moi'])
