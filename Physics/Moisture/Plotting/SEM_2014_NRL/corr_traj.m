close all
clear
clc
mydir = pwd;
fontsize = 16;


cd /discover/nobackup/drholdaw/ExperimentData/BacmeisterPaper/Correlations/

cfn_vs_all

qs_cfn = qs; clear qs
qt_cfn = qt; clear qt
rh_cfn = rh; clear rh
q_cfn = q; clear q
clear i

qcn_vs_all

qs_qcn = qs; clear qs
qt_qcn = qt; clear qt
rh_qcn = rh; clear rh
q_qcn = q; clear q
clear i

cd(mydir)


xmin = 0.8;
xmax = 1.2;
RH1 = xmin:0.0001:xmax;

qoff = 0.0665;
q_min = 1-qoff;
q_max = 1+qoff;
grad = 1/(q_max - q_min);
q1 = 3*grad;

CF_tanh = 0.5 + 0.5*tanh( q1*(RH1 - 1.0)); 

a = q_max/q_min;
b = (0.5*(a+1)*q_min ) ;

CF_lin = zeros(1,length(RH1));
for i = 1:length(RH1)
    
    if RH1(i) < q_min 
        CF_lin(i) = 0;
    elseif q_min <= RH1(i) && RH1(i) < q_max
        CF_lin(i) = RH1(i) * 1/((a-1)*q_min) - 1/(a-1);
    elseif RH1(i) >= q_max
        CF_lin(i) = 1;
    end
        
end

qmq = -2e-4:3e-4/100:1e-4;
qcn_lin = zeros(1,length(qmq));
for i = 1:length(qmq)
    
    if qmq(i) < 0 
        qcn_lin(i) = 0;
    else
        qcn_lin(i) = qmq(i) ;
    end
        
end

qcn_nlin = zeros(1,length(qmq));
for i = 1:length(qmq)
    
    if qmq(i) < -0.25e-4
        qcn_nlin(i) = 0;
    elseif -0.25e-4 <= qmq(i) && qmq(i) < 0.25e-4
        qcn_nlin(i) = 10000*((qmq(i) + 0.25e-4).^2 );
    elseif qmq(i) >= 0.25e-4
        qcn_nlin(i) = qmq(i) ;
    end
        
end



grey1 = 0.75;
grey2 = 0.5;
lin_wid = 2.0;

figure
set(gcf,'position',[494   498   785   840])

subplot(2,1,1)
scatter(rh_cfn,cfn,0.025,'k')
hold on
plot(RH1,CF_lin ,'r','LineWidth',lin_wid)
plot(RH1,CF_tanh,'g','LineWidth',lin_wid)
box on
xlabel('RH','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('C_{LS}','FontSize',fontsize,'FontName','TimesNewRoman')
title('(a)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([1-0.10 1+0.10])
ylim([-0.01 1.01])
legend('Model Output','Linear','Nonlinear','Location','NorthWest');
legend boxoff

grey1 = 0.5;
grey2 = 0.75;

subplot(2,1,2)
scatter(qt_qcn-qs_qcn,qcn,0.025,'k')
hold on
plot(qmq,qcn_lin ,'r','LineWidth',lin_wid)
plot(qmq,qcn_nlin,'g','LineWidth',lin_wid)
box on
xlabel('q_T-q_s (kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
ylabel('q_{LS} (kgkg^{-1})','FontSize',fontsize,'FontName','TimesNewRoman')
title('(b)','FontSize',fontsize,'FontName','TimesNewRoman')
set(gca,'FontSize',fontsize,'FontName','TimesNewRoman')
xlim([-2e-4 2e-4])
ylim([-0.01e-4 1.01e-4])
legend('Model Output','Linear','Nonlinear','Location','NorthWest');
legend boxoff
