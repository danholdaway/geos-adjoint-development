close all
clear
clc
mydir = pwd;


cd /discover/nobackup/drholdaw/btmp.25956/sens.20140202.000000/

% cfnd_vs_sig1
% cfnd_vs_all
% cfnd1 = cfnd; clear cfnd
% cfn1 = cfn; clear cfn
% qcnd_vs_sigma
% tend_vs_af
% qaxd_vs_af1
% cfnd_vs_all1
% qcnd_vs_bigjac
% cf_vs_eval4
% cf_vs_ymax
% cfd_vs_cf
% cflsd_vs_cfls
% qllsd_vs_eval_187
% all_evals
% cfn_vs_all
% qcnd_vs_rh

% cfnd_vs_all_new
% qcnd_vs_all_new
% dqcnd_vs_all_new

% cfnd_vs_lots7
% qcnd_vs_lots2
% dqcnd_vs_lots2

corr_jacob2

% qcn_vs_qtqs

cd(mydir)

% term1 = ((qt+sigma-qs)./(2*sigma));
% term2 = -sigmad.*(qt+sigma-qs)./((2*sigma).^2);

figure
scatter(dcflsd,J71,'kx')

asd

% asd

% figure
% scatter(cfnd,term2,'kx')
% a = get(gca,'xLim');
% b = get(gca,'yLim');


% asd

D1 = 1./sigma;
qcnd1 = dqcnd;
alpha1 = sigma./qs;
sigma1 = sigma;

count = 0;
for i = length(D1):-1:1
    
    if D1(i) < 7500

        D1(i) = [];
        qcnd1(i) = [];
        alpha1(i) = [];
        sigma1(i) = [];
        count = count + 1;
        
    end
    
end

figure
scatter(qcnd1,D1,'kx')


figure
scatter(qcnd1,alpha1,'kx')

figure
scatter(qcnd1,sigma1,'kx')
% xlim([a(1) a(2)])
% ylim([b(1) b(2)])


% set(gca,'YScale','log')
% ylim([10^0 10^0.4])
      
asd

% all_evals = [eval1 eval2 eval3 eval4 eval5 eval6 eval7 eval8];
% 
% eval = sort(all_evals,'ascend');
% 
% plot(eval)
% 
% asd



eval11 = eval1;
eval21 = eval2;
eval31 = eval3;
eval41 = eval4;
eval51 = eval5;
eval61 = eval6;
eval71 = eval7;
eval81 = eval8;

for i = 1:length(eval1)
   
    if eval11(i) == 1
        eval11(i) = 0;   
    end
    if eval21(i) == 1
        eval21(i) = 0;   
    end
    if eval31(i) == 1
        eval31(i) = 0;   
    end
    if eval41(i) == 1
        eval41(i) = 0;   
    end
    if eval51(i) == 1
        eval51(i) = 0;   
    end
    if eval61(i) == 1
        eval61(i) = 0;   
    end
    if eval71(i) == 1
        eval71(i) = 0;   
    end
    if eval81(i) == 1
        eval81(i) = 0;   
    end
    
end

maxeval = zeros(1,length(eval1));
for i = 1:length(eval1)
    maxeval(i) = max(abs([eval1(i) eval2(i) eval3(i) eval4(i) eval5(i) eval6(i) eval7(i) eval8(i)]));
end
maxeval1 = zeros(1,length(eval1));
for i = 1:length(eval1)
    maxeval1(i) = max(abs([eval11(i) eval21(i) eval31(i) eval41(i) eval51(i) eval61(i) eval71(i) eval81(i)]));
end

figure
scatter(qllsd,maxeval,'kx')


% 
% count1 = 0;
% count2 = 0;
% for i = 1:length(denom)
%     
%     if denom(i) > 1.15
%         count1 = count1 + 1;
%     else
%         count2 = count2 + 1;
%     end
% end
% 
% 100*count1/(count1+count2)

% scatter(abs(qt-qstar),(qt+sigma-qstar)./(2*sigma),0.2,'.')
% scatter(cfnd,sigma./qstar,0.2,'.')
% set(gca,'YScale','log')

% x = -6:0.1:2;
% y = ((1)*x + (-0.0));
% hold on
% % plot(x,y,'r')
% 
% cfnd1 = cfnd;
% 
% counter = 0;
% for i = 1:length(cfnd)
%     if abs(cfnd(i) - -((-0.15)*x1(i) + (-0.7))) > 1.0
%         counter = counter + 1;
%     end
% end
% disp(100*counter/i)

% asd

% figure
% scatter(qcnd,(J2),0.1,'.')
% 
% figure
% scatter(qcnd,(J3),0.1,'.')
% 
% figure
% scatter(qcnd,(J4),0.1,'.')
% 
% figure
% scatter(qcnd,(J5),0.1,'.')
% 
% figure
% scatter(qcnd,(J1)*1e2 + (J2)*1e-5 + (J3)*1e-5 + (J4)*1e-5 + (J5)*1e0 ,0.1,'.')

% figure
% scatter(qcnd,(J2),0.1,'.')
% 
% figure
% scatter(qcnd,abs(J1)+abs(J2)+abs(J3),0.1,'.')


% scatter(sig1,cfd,0.1,'.')
% 
% x = 0:0.1:1;
% y = -((-0.15)*x + (-0.7));
% 
% scatter(sigma,sigmad,0.1,'.')
% 
% figure
% scatter(qstar,qstard,0.1,'.')
% 
% figure
% scatter(qt,qtd,0.1,'.')
% 
% figure
% scatter(cfn1,cfnd1,0.1,'.')


% 
% figure
% scatter(qcnd,qt,0.1,'.')
% 
% figure
% scatter(qcnd,qstar,0.1,'.')
% 
% figure
% scatter(qcnd,(qt+sigma-qstar)./(2*sigma),0.1,'.')
% 
% figure
% scatter((qt+sigma-qstar)./(2*sigma),abs(cfnd),0.1,'.')
% hold on
% plot(x,y,'r')
% 
% x1 = (qt+sigma-qstar)./(2*sigma);
% 
% counter = 0;
% for i = 1:length(cfnd)
%     if abs(cfnd(i) - -((-0.15)*x1(i) + (-0.7))) > 1.0
%         counter = counter + 1;
%     end
% end
% disp(100*counter/i)


% scatter((4*sigmaqt1).^2,qcnd,0.5,'.')


% scatter(af,temp,0.1,'.')
% 
% figure
% scatter(af,tempd,0.1,'.')

% scatter(af,qaxd,2,'.')
% 
% qaxd1 = qaxd;
% 
% counter = 0;
% for i = 1:length(qaxd)
%     if af(i) < 5e-4
%         counter = counter + 1;
%         qaxd1(i) = 0;
%     end
% end
% disp(100*counter/i)
%         
% figure
% scatter(af,qaxd1,2,'.')
% 
% figure
% scatter(af,qax,1,'.')