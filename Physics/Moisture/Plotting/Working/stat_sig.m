close all
clear
clc

% generate data:
n1 = 100; m1 = 10.0; s1 = 0.25;
n2 = 100; m2 = 12.0; s2 = 0.50;
data1 = m1 + s1 * randn(n1,1);
data2 = m2 + s2 * randn(n2,1);
%data2 = data1; % test with same data

% get the mean and variance of the sample
mdata1 = mean(data1); vdata1 = var(data1);
mdata2 = mean(data2); vdata2 = var(data2);
% data2 = data1;

scatter(1:n1,data1)
hold on
scatter(1:n2,data2,'r')
plot(1:100,mdata1*ones(1,n1),'-')
plot(1:100,mdata2*ones(1,n2),'-r')

% compute T value
dmdata = mdata2 - mdata1;
sdata = sqrt( vdata1 / n1 + vdata2 / n2 );
tval = dmdata / sdata;

% choose confidence interval
conf_int = 95.0;

ttest2(data1,data2)

% get tcrit for N (pooled sample size) and confidence interval ( 2-tailed test )
% tcrit = stats.t.ppf(1.0 - (1.0-conf_int/100.0)/2.0, (n1 + n2 - 2))

% Test for statistical significance
% Null hypothesis : The two means are essentially the same
% ssig = np.sign(tval-tcrit) != np.sign(tval+tcrit)
% if ( ssig ):
%     print 'the two samples are STATISTICALLY THE SAME at %d%% confidence interval' % int(conf_int)
% else
%     print 'the two samples are STATISTICALLY DIFFERENT at %d%% confidence interval' % int(conf_int)
% end
