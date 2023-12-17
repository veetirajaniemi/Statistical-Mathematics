%% Exercise 2
%% Task 1
clearvars;close all;clc
            
old = [33.4697;28.5377;35.7013;36.3507;		
30.6141;39.9637;39.9567;35.0495;		
36.5092]

new = [30.5335;29.0158;32.8483;27.3291;
38.9694;29.2271;30.2785;34.2804;30.0490;
29.3983;26.3041;31.0365]

% Let's test for normality first.
h1 = adtest(old)
h2 = adtest(new)

% We get 0 (false) so datas are normally distributed (failed to reject
% the null hypothesis). If the normality test failed, we would probably
% have some measure errors with samples. 

% Null hypothesis: H_0 : mean(new)-mean(old) = 5
% H_A = mean(old)-mean(new) != 5
% (A = old, B = new in lecture material)
% Variances eq. but unknown -> t-procedure
n_old = length(old);
n_new = length(new);
meanold = mean(old);
meannew = mean(new);
d = 5;

sdev_old = sqrt(sum((old-meanold).^2)/(n_old-1));
sdev_new = sqrt(sum((new-meannew).^2)/(n_new-1));

sp2 = ((n_new-1)*sdev_new^2+(n_old-1)*sdev_old^2)/(n_new+n_old-2);
sp = sqrt(sp2);
se = sp*sqrt(1/n_new+1/n_old);

% Let's calculate t-value
t = (meanold-meannew-d)/se
df = n_new+n_old-2
pval = 2*tcdf(t,df)

% P-value 0.3430 is over 0.2 so we accept our null hypothesis. Producer's
% claim is true. 


%% Task 2
clearvars;close all;clc

n_old = 40;
n_new = 54;
svar_old = 17.8;
svar_new = 26.2;

% We don't know the theoretical values of variances. 
% We have the sample variances of the data so let's use F-test.
% Null hypothesis H_0: svar_old = svar_new, H_A = svar_old != svar_new

%  n_old - 1 and n_new - 1 degrees of freedom
df = [n_old-1 n_new-1];
F = svar_old/svar_new % under 1

pval = 2*fcdf(F,df(1),df(2)) % two-sided test

% P-value is just over 0.2 so we accept our null hypothesis.
% We can assume equality of variances.


%% Task 3
clearvars;close all;clc

% The expected weight of a sugar packet is 1 kg. Significance level is 1
% -0.99 = 0.01. 
n = 16;
exp = 1;
mean = 1.053;
s = 0.058;
a = 0.01;

% A mean variable T is t-distributed because variable isn't known. Degree of
% freedom is n-1=15.
% T = (mean-exp)/(s/sqrt(n)) ~ t(15)

% Now the confidence interval is mean +- t_(1-a/2)*s/sqrt(n) (Beta 18.3.2)

t_crit =  tinv(1-a/2,n-1) % critical value
m_lower = mean-t_crit*s/sqrt(n)
m_upper = mean+t_crit*s/sqrt(n)


%% Task 4
clearvars;close all;clc

% H0: p >= 0.1, H_A : p < 0.1, p is probability of false-positive results

% 23 cases out of 324 positive results turned to be false-positive
 
n = 324; % amount of cases
x = 23; % positive cases
p_exp = x/n; % expected probability
p0 = 0.1; % prob value we are comparing to 
a = 1-0.99;

z = (x-n*p0+0.5)/sqrt(n*p0*(1-p0)); % test variable
pval = normcdf(z) % p-value

% P-value is 0.0497. That means we reject null hypothesis --> probability
% of false-positive results is under 0.1 so screening test is acceptable.

% Then upper confidence bound. We can use normal distribution approximation
% because of big enough sample size (n*p0*(1-p0) > 9)

zcrit = norminv(1-a);
p_upper = p_exp+zcrit*sqrt(p_exp*(1-p_exp)/n) % (Beta 18.3.2 applied for one-sided)

% 99% upper confidence bound is 0.1042


%% Task 5
clearvars;close all;clc

acc = 0.02;
a = 1-0.99;
p = 0.5; 
p2 = 0.4;

% When we want to find proportion with accuracy +-2%, we know that the
% confidence interval length L is 2*acc = 2*2% = 2*0.02 = 0.04
L = 2*acc;

% L = 2*z_(a/2)*sqrt((p*1-p)/n) ---> good sample size would be
n1 = 4*norminv(a/2)^2*p*(1-p)/L^2 % with proportion p = 0.5
n2 = 4*norminv(a/2)^2*p2*(1-p2)/L^2 % with proportion p = 0.4


%% Task 6
clearvars;close all;clc

% x is count of dropped packages, n is number of trials
x_a = 35;
n_a = 44;
x_b = 36;
n_b = 52;
p_a = x_a/n_a;
p_b = x_b/n_b;
a = 1-0.99;

% A two-sided confidence interval for p_a-p_b is p_a-p_b +- x, when 
x = norminv(a/2)*sqrt(p_a*(1-p_a)/n_a+p_b*(1-p_b)/n_b);
% We notice that x is negative --->
pdiff_lower = p_a-p_b+x
pdiff_upper = p_a-p_b-x

% Then the hypothesis testing:
% H0: p_a = p_b, H_a: p_a != p_b

% Statistic value z ~ N(0,1)
p_exp = (x_a+x_b)/(n_a+n_b);
z = (p_a-p_b)/sqrt(p_exp*(1-p_exp)*(1/n_a+1/n_b));

% Two-sided test ->
pval = 2*(1-normcdf(z)) % (z is a positive value)
% P-value 0.2512 > 0.2  --> we accept H0 so the radar systems are equally
% effective.

