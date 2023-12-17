%% Harjoitus 1
%% Tehtävä 1
clearvars;close all;clc

% a)
n_1 = 20;
n_2 = 18;
e1 = 7.2; % expected values
e2 = 8.5;
s1 = 4.2; % standard deviations
s2 = 6.1;

% creating samples
sample1 = s1*randn(n_1,1)+e1;
sample2 = s2*randn(n_2,1)+e2;

% b)
% Calculating means, difference, variances and quotient
m1 = mean(sample1)
m2 = mean(sample2)
D = m1-m2
var1 = sum((sample1-m1).^2/(n_1-1))
var2 = sum((sample2-m2).^2/(n_1-1))
Q = var1/var2

% c) Let's do calculations multiple times using a for loop. I think that the 
% difference is near to the value e1-e2 = -1.3 and that I could see it clearer 
% with bigger samples. Quotient should be near s1^2/s2^2 = 0.47.

n = 10000;
n2 = 50;
n3 = 1000;

D = zeros(n,1);
Q = zeros(n,1);
for i = 1:n
    Sample1 = s1*randn(20,1)+e1;
    Sample2 = s2*randn(18,1)+e2;
    M1 = mean(Sample1);
    M2 = mean(Sample2);
    D(i) = M1-M2;
    var1 = sum((Sample1-M1).^2/(n_1-1));
    var2 = sum((Sample2-M2).^2/(n_1-1));
    Q(i) = var1/var2;
end

% Let's create histograms with different n values. 

histogram(D(1:n2))
title('Difference, n = 50')
figure
histogram(D(1:n3))
title('Difference, n = 1000')
figure
histogram(D)
title('Difference, n = 10000')
figure
histogram(Q(1:n2))
title('Quotient, n = 50')
figure
histogram(Q(1:n3))
title('Quotient, n = 1000')
figure
histogram(Q)
title('Quotient, n = 10000')

% We can see that the histograms look pretty much like what was expected
% :). Results are also clearer when samples are bigger.

%% Tehtävä 2
clearvars;close all;clc

load examgrades;
n = 5; % 5 columns
h = zeros(1,n); % test results
p = zeros(1,n); % p-values
res = zeros(n,2); % matrix to store results
format long % to see p-values better
for i = 1:n
    x = grades(:,i); % one column at a time
    test_cdf = [x,cdf('tlocationscale',x,75,10,1)];
    [res(i,1), res(i,2)] = kstest(x,'CDF',test_cdf);
end
res

% From res-matrix we can see that the result is always 1 (=false) so the 
% grades don't come from Student's t distribution. We can also see that 
% every p-value is under 0.001 which confirms the results. 

%% Tehtävä 3
clearvars;close all;clc

% a)
load 'RandD.mat';
parexp = expfit(RandD); % parameters for distributions
parwbl = wblfit(RandD);

% b)
H = histogram(RandD,'Normalization','pdf');
x = linspace(0,0.12,length(H.Data));
pdwbl = makedist('Weibull',parwbl(1),parwbl(2)); % creating distribution
ywbl = pdf(pdwbl,x); % and prob density function
hold on
plot(x,ywbl,'r','LineWidth',2)

pdexp = makedist('Exponential', parexp);
yexp = pdf(pdexp,x);
hold on
plot(x,yexp,'g','LineWidth',2)

legend('Data','Theoretical Weibull', 'Theoretical exponential')

% c) Let's check test results (and p-values for fun). 
[he, p_exp] = kstest(RandD, 'CDF', pdexp)
[hw, p_wbl] = kstest(RandD, 'CDF', pdwbl)

% Both h values are true, so we could say that data can come from both
% distributions.

% d)
figure
cdfplot(RandD); 
hold on
cdfexp = cdf(pdexp,x); 
cdfwbl = cdf(pdwbl,x);
plot(x,cdfexp,'r')
hold on
plot(x,cdfwbl,'g')
legend('Empirical CDF', 'Hypothesized exponential', 'Hypothesized Weibull','location', 'best')

