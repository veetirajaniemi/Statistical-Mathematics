%% Exercise 3
%% Task 1
clearvars;close all;clc

test1 = [18;21;16;22;19;24;17
21;23;18;14;16;16;19;18;20;12;
22;15;17]; % before studying

test2 = [22 25 17 24 16 29 20 23 19 ...
20 15 15 18 26 18 24 18 25 19 16]'; % after studying
n = length(test1);

% We can reduce the problem to a one-sample problem
% Differences:
Z = test1-test2;
mean = mean(Z)
% Confidence intervals, now t procedure (we don't know the variance)
a1 = 1-0.95;
a2 = 1-0.99;
s = std(Z)/sqrt(n);
zlower1 = mean-tinv(1-a1/2,n-1)*s
zupper1 = mean+tinv(1-a1/2,n-1)*s
zlower2 = mean-tinv(1-a2/2,n-1)*s
zupper2 = mean+tinv(1-a2/2,n-1)*s

% Then the testing with hypothesis
% H_0: mean(Z) >= 0, HA: mean(Z) < 0

% Test statistic:
t = ((mean-0)/std(Z))*sqrt(n) % negative
pval = tcdf(t,n-1)

% We can see that the p-value 0.0022 < 0.005 so we reject the null hypothesis.
% It means that teaching has improved the students' scores!

%% Task 2
clearvars;close all;clc

bottom = [0.4300 0.2660 0.5670 0.5310... 
    0.7070 0.7160 0.6510 0.5890 0.4690 0.7230]';

surface = [0.4150 0.2380 0.3900 0.4100 0.6050...
    0.6090 0.6320 0.5230 0.4110 0.6120]';

n = length(bottom) % same amount of samples

% H_0: mean(bottom)-mean(surface) < 0,
% H_A: mean(bottom)-mean(surface) >= 0

Z = bottom-surface;
mean = mean(Z)

% We don't know variance --> t-statistic
% Test statistic:

t = ((mean-0)/std(Z))*sqrt(n) % positive
pval = 1-tcdf(t,n-1)

% P-value 0.00044556 < 0.001 --> we reject the null
% hypothesis. It means that the average zinc concentration
% is positive so true average concentration in bottom water
% exceeds that of the surface water.


%% Task 3
clearvars;close all;clc

A = [39 51 58 61 65 72 86]';
B = [22 38 43 47 49 54 72]';
C = [18 31 41 43 44 54 65]';
meanA = 61.71;
meanB = 46.43;
meanC = 42.29;
varA = 225.24;
varB = 232.95;
varC = 229.24;
n = 7; % same for all groups

% We have to perform the test because it's possible that one single value
% can change the value massively so we can't really know if the higher mean
% is actually because of that or because of the policy type. The group
% sizes are also very small which highlihts the individual's effect on
% values.

hA = logical(jbtest(A)) % 0
hB = logical(jbtest(B)) % 0
hC = logical(jbtest(C)) % 0
% All values are normally distributed. Had to change the return values to
% logical to get correct values on tests when returning the task.

% If one of the tests failed so some values weren't normally distributed, we 
% couldn't use ANOVA and would have to stick with some other testing, for
% example some non-parametric tests.

% However, no we can use ANOVA:
% H0: meanA = meanB = meanC
% HA: mean(i) != mean(j) for some i and j (in A,B or C)

% Let's do the anova tests:
[pval, table] = anova1([A B C])
SSD = table{2,2} % variability between groups
SSE = table{3,2} % variability within groups
SST = table{4,2} % total variability in data
F = table{2,5}

% P-value is 0.0648 < 0.2 so we reject the null hypothesis. That means that
% there is some kind of difference between the groups.


%% Task 4
clearvars;close all;clc

nofer = [3.5 4.1 4.0 5.1 4.8 3.5]';
fer1 = [4.2 3.7 5.5 4.7 3.9 4.2]';
fer2 = [6.2 4.8 4.9 6.0 3.9 5.2]';
n = length(nofer);

% Normality tests
h0 = logical(jbtest(nofer))
h1 = logical(jbtest(fer1))
h2 = logical(jbtest(fer2))

% All datas are normally distributed. We can use ANOVA testing:
% H0: mean(nofer) = mean(fer1) = mean(fer2)
% HA: mean(i) != mean(j) with some i and j (in nofer,fer1,fer2) 
[pval, table] = anova1([nofer fer1 fer2])
SST = table{4,2}
SSD = table{2,2}
SSE = table{3,2}
F = table{2,5}

% P-value is 0.0697 < 0.2 so we reject the null hypothesis. There is
% difference between mean annual growth between the three plant
% populations. 


%% Task 5
clearvars;close all;clc

dist = [0 50 150 200 250 300 350 400 450 500];
depth = [0 10 28 42 59 51 73 85 104 96];

% H0: p (rho) = 0, HA: p (rho) != 0
% Let's rank the values 

S1 = sort(dist); % sorting data
S2 = sort(depth);

for i = 1:length(dist) % calculating ranks
    r1(i) = find(S1 == dist(i));
    r2(i) = find(S2 == depth(i));
end

diff = (r1-r2)';
n = length(diff);

% let's calculate rho and test statistic
rho = 1-6*sum(diff.^2)/(n*(n^2-1))
t = rho*sqrt((n-2)/(1-rho^2)) % positive
pval = 2*(1-tcdf(t,n-2))

% P-value 0.0000014675 is very small (<0.2) so we
% reject our null hypothesis. There is some relationship between
% river depth and the distance from the river bank


%% Task 6
clearvars;close all;clc

students = [40 45 51 58 66 68 74 81]';
youngad = [40 49 38 62 80 65 82 71]';

% H0: rho = 0, HA: rho != 0

% Let's sort and calculate ranks

S1 = sort(students);
S2 = sort(youngad);

for i = 1:length(students)
    r1(i) = find(S1(i)==students);
    r2(i) = find(S2(i)==youngad);
end

diff = (r1-r2)';
n = length(diff);

rho = 1-6*sum(diff.^2)/(n*(n^2-1))
t = rho*sqrt((n-2)/(1-rho^2)) 
pval = 2*(1-tcdf(t,n-2))

% P-value 0.0149 is under 0.2. We reject our null hypothesis so there is
% some correlation between the opinions. 

