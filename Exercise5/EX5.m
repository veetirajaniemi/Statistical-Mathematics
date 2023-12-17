%% Exercise 5
%% Task 1
clearvars;close all;clc

n1 = 40;
n2 = 36;
mean1 = 12.5;
mean2 = 18.3;
sdev1 = 2.2;
sdev2 = 1.8;

% a) When comparing variances of two samples, we can use F-statistic.
% H0: var1 = var2, H_s: var1 != var 2
var1 = sdev1^2;
var2 = sdev2^2;

F = var1/var2
pvalF = 2*(1-fcdf(F,n1-1,n2-1)) % two-sided, F-statistic over 1

% The p-value is 0.2315, so I can clearly accept the null hypothesis. The
% variances of the samples are equal. 

% b) Now we have theoretically unkown but equal variances. Let's calculate
% the t-statistic with this information. 

% H0: mean2-mean1 = 5, H_s: mean2-mean1 != 5

sp2 = ((n2-1)*sdev2^2+(n1-1)*sdev1^2)/(n2+n1-2); 
sp = sqrt(sp2);
se = sp*sqrt(1/n2+1/n1); % standard error for t statistic

t = (mean2-mean1-5)/se
df = n2+n1-2;
pvalt = 2*(1-tcdf(t,df))

% The p-value of 0.089 > 0.05 so I would accept the null hypothesis.
% Shielding has improved the average time by 5 minutes. 


%% Task 2
clearvars;close all;clc

x = 22; % good throws
n = 42; % attempted throws

% H0: p <= 0.32, H_A: p > 0.32

p_exp = x/n;
p0 = 0.32;

% Let's calculate z-statistic. 
z = (p_exp-p0-0.5/n)/sqrt((p0*(1-p0))/n)
pval = 1-normcdf(z)

% P-value 0.0038 < 0.01 so we can clearly reject our null hypothesis. The 
% player has improved her free throw accuracy. 

%% Task 3
clearvars;close all;clc

A = [7 9 5 8 6 8 6 10 7 4]';
B = [4 3 6 2 7 5 5 4 1 3]';
C = [6 1 3 5 3 4 6 5 7 3]';
n = length(A) % same for all groups

jbtest(A); % 0
jbtest(B); % 0
jbtest(C); % 0

% Test scores are normally distributed so we can use anova testing! 
% H0: mean(A) = mean(B) = mean(C), H_A: mean(i) != mean(j) for some i,j in
% A,B,C

[pval, table] = anova1([A B C])
F = table{2,5}

% The p-value 0.0017 < 0.01 so it's easy to say that we can reject our null
% hypothesis, which tell us that some of the expected values is different
% from others so music affects students. It seems like the students listening 
% music have the best results (even though we can't say that after this testing). 


%% Task 4
clearvars;close all;clc

snif = [14 21 18 19 12 14 23 24 16]';
lab = [55 63 67 59 22 31 80 68 54]';
n = length(snif);
% Let's see if there is any correlation between the results with both
% methods. H0: rho = 0, H_A: rho != 0

S2 = sort(lab);
for i = 1:n % calculating ranks
    r2(i) = find(S2 == lab(i));
end
r1 = [2.5 7 5 6 1 2.5 8 9 4]; % easy to do manually
% duplicate values' ranks are mean values of all of their ranks

diff = r1-r2;

rho = 1-(6*sum(diff.^2))/(n*(n^2-1))
t = rho*sqrt((n-2)/(1-rho^2)) % positive
pval = 2*(1-tcdf(t,n-2)) % two-sided

% P-value is very small, so we can easily reject the null hypothesis. That
% means that there is clear correlation between the results, which means
% that both methods produce same kind of results. 


%% Task 5
clearvars;close all;clc

imp = importdata('dataset1.xlsx');
data = imp.data;

% Let's see if there is any correlation between the age of an antique clock
% and the number of bidders or it's selling price. 

[R,P] = corrcoef(data)

% It seems like the number of bidders doesn't correlate very well with the age of the
% clock. The only clearly variable which correlates is the price of the
% clock (index 3). Let's form the matrix X with ones (constant val in model) and 
% the other variables. Y is the dependent variable (so the age of the
% clock).

X = [ones(size(data,1),1) data(:,3)];
Y = data(:,1);
b = X\Y % estimated parameters to linear regression

p = length(b); % number of parameters
n = size(data,1); % number of rows

Yfit = X*b;
res = Y-Yfit; % residuals

RSS = sum(res.^2); % residual sum of squares
var = RSS./(n-p); % variance

% Now let's calculate the t-value of the coefficient:
cov = inv(X'*X)*var; % covariances
std = sqrt(diag(cov)); % standard deviations
t = b./std
pval = 2*(1-tcdf(abs(t),n-p)) % p-values to each parameter

results = [b t pval]

TSS = sum((Y-mean(Y)).^2); % total sum of squares
R2 = 1-RSS/TSS
% 0.5332 is not very good

SSR = sum((Yfit-mean(Y)).^2);
F = (SSR/(p-1))/(RSS/(n-p));
pval_f = 1-fcdf(F,p-1,n-p)

stres = res./sqrt(var); % no outliers!

figure 
plot(data(:,3),stres,'r*')

% Results from MATLAB functions look same.
table = readtable('dataset1.xlsx');
model = fitlm(table, 'Age~Price')