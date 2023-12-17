%% Exercise 4
%% Task 2 - dataset 1
clearvars;close all;clc

imp = importdata('data1.xlsx');
data = imp.data;

% Let's see is there any correlation between distance travelled in feet and
% other values from the data.

[R,P] = corrcoef(data); % correlation coefficients
i = find(P(:,1)>0.05); % let's decide that p values under 0.05 are significant 
% We can exclude these data values. It seems like all punters were
% right-footed.

inc = setdiff(setdiff(1:size(data,2),i),[2 4 5])
% data in these indexes is included!

[R1,P1] = corrcoef(data(:,inc));
% correlation coefficients and p-values for included data

% The data in column 4 is deleted (smallest correlation to distance variable.
% I tried to delete other data columns as well but couldn't really get to
% any better situation so decided to continue with these values. 

X = [data(:,inc)];
Y = data(:,1);
b = X\Y; % estimated parameters to linear regression 

p = length(b); % number of parameters
n = size(data,1); % number of rows

Yfit = X*b;
res = Y-Yfit; % residuals

RSS = sum(res.^2); % residual sum of squares
var = RSS./(n-p); % variance

% Now let's calculate so called signal to noise ratio or the t-value of the
% coefficient:
cov = inv(X'*X)*var; % covariances
std = sqrt(diag(cov)); % standard deviations
t = b./std;
pval = 2*(1-tcdf(abs(t),n-p)); % p-values to each parameter

results = [b t pval]
% I removed parameters 2 and 5 because of the biggest p-values. I also removed
% ones from X so the constant from the model. 

% Let's calculate the R2 and F value to see the actual goodness of the model

TSS = sum((Y-mean(Y)).^2); % total sum of squares
R2 = 1-RSS/TSS
% 0.7789 is pretty good value

SSR = sum((Yfit-mean(Y)).^2);
F = (SSR/(p-1))/(RSS/(n-p));
pval_f = 1-fcdf(F,p-1,n-p)

% The p-value of our F statistic is really small so our model is very good!

stres = res./sqrt(var); % easy to check that there are no outliners, absolute value
% of every value is under 3

figure 
for i = 1:length(inc)
    subplot(length(inc),1,i)
    plot(data(:,inc(i)),stres,'r*')
end
% Plots look good, there is not anything worrying.  


%% Task 2 - dataset 2
clearvars;close all;clc

imp = importdata('data2.xlsx');
data = imp.data;

% Let's see is there any correlation between number of cigarettes smoked
% and other values from the data. 

[R,P] = corrcoef(data); % correlation coefficients
i = find(P(:,1)>0.05); % let's decide that p values under 0.05 are significant 
% We can exclude these data values. It also feels correct that leukemia and
% number of cigarettes smoked don't have any correlation between them. 

inc = setdiff(1:size(data,2),i)
% data in these indexes is included!

[R1,P1] = corrcoef(data(:,inc));
% correlation coefficients and p-values for included data

% The correlations look fine so I don't delete any columns. 

X = [data(:,inc)];
Y = data(:,1);
b = X\Y; % estimated parameters to linear regression 

p = length(b); % number of parameters
n = size(data,1); % number of rows

Yfit = X*b;
res = Y-Yfit; % residuals

RSS = sum(res.^2); % residual sum of squares
var = RSS./(n-p); % variance

% Now let's calculate so called signal to noise ratio or the t-value of the
% coefficient:
cov = inv(X'*X)*var; % covariances
std = sqrt(diag(cov)); % standard deviations
t = b./std;
pval = 2*(1-tcdf(abs(t),n-p)); % p-values to each parameter

results = [b t pval]
% I removed ones from X (so the constant from the model) because the p-value
% 0.7246 is really big.

% Let's calculate the R2 and F value to see the actual goodness of the model

TSS = sum((Y-mean(Y)).^2); % total sum of squares
R2 = 1-RSS/TSS
% 0.6443 is pretty good value so the model is alright. I could also remove
% second parameter but there is too much difference since R2 would be 0.57
% without that value. 

SSR = sum((Yfit-mean(Y)).^2);
F = (SSR/(p-1))/(RSS/(n-p));
pval_f = 1-fcdf(F,p-1,n-p)

% The p-value of our F statistic is really small so our model is again very good!

stres = res./sqrt(var); % we see there's one outlier! it's only one value which is
% just over 3 so it shouldn't make a huge impact on this model
% I also tried to remove other parametres but when doing that there were
% worse outliers.

figure 
for i = 1:length(inc)
    subplot(length(inc),1,i)
    plot(data(:,inc(i)),stres,'r*')
end

% We can see the outlier but otherwise plots look good and there is not anything worrying.  

