%% Exercise 6

%% Dataset 1
clearvars;close all;clc
all = readmatrix('DataA.xlsx');
data_1 = all(9:end-3,2);

data = data_1(212:319); % the data to analyze
plot(data)

grid on
hold on

X = [ones(length(data),1) [1:length(data)]'];
Y = data;
b = X\Y;

trendLine = X*b;
plot(1:length(data),trendLine,'k')
legend('Whole data to be analysed', 'Trendline', 'location', 'best')

data_s = data - trendLine; % detrended

% There is some kind of seasonality in data.
p = 12; % seasonal indexes

for i = 1:p
    S(i) = mean(data_s(i:p:end));
end

S_mul = repmat(S',9,1); % Nine years of data
figure
plot(data_s,'r')
hold on
grid on
plot(S_mul,'b')
legend('Detrended data', 'Seasonal indexes','location','best')

I = data_s-S_mul; % irregular parts
figure
grid on
hold on
plot(I)
legend('Irregular parts')

h = kpsstest(I) % 1, not stationary
I_diff = diff(I);
h_d = kpsstest(I_diff) % 0, now it's stationary :)

figure
autocorr(I_diff,40)
figure
parcorr(I_diff, 40)


models = [];
for i = 0:4 % Let's calculate 5*5 different models 
    for j = 0:4
        model = armax(I_diff, [i j]);
        a = aic(model);
        models = [models; i j a];
    end
end

ind = find(models(:,end) == min(models(:,end)));
min = models(ind,:);

% Minimum aic value vith 2 and 3.
ModelFit = armax(I_diff, [min(1) min(2)]);

I_df = forecast(ModelFit,I_diff,p); % forecast for _differences_!

figure
plot(I_diff), hold on
plot(length(I_diff)+1:length(I_diff)+p,I_df,'r')
legend('Differences', 'Forecasts of differences', 'location', 'best')

Xf = [ones(p,1) [length(data)+1:length(data)+p]'];
Tf = Xf*b; % Trend forecast

If(1) = I(end)+I_df(1); % Forecast of irregular series

for i = 2:p
    If(i) = If(i-1)+I_df(i);
end

forecast = Tf + S' + If'; 

figure % had to do little tricks with indexes because choosed them
% manually... not optimal but it works
plot(1:331,data_1(1:331),'b'), hold on, grid on % original data
plot(320:331,forecast,'r') % forecast 
legend('Original data', 'Forecast','location','best')


%% Dataset 2
clearvars;close all;clc

all = readmatrix('DataB.xlsx');
data_1 = all(9:end-3,2); % all data
data = all(239:334,2); % the data to analyse

plot(data), hold on, grid on
X = [ones(length(data),1) [1:length(data)]'];
Y = data;
b = X\Y;

trendLine = X*b;
plot(1:length(data),trendLine,'k')
legend('Whole data to be analysed', 'Trendline', 'location', 'best')
data_s = data - trendLine; % detrended data

% There's clear seasonality in this data! We know that the data is monthly.

p = 12;
for i = 1:p
    S(i) = mean(data_s(i:p:end)); % Seasonal indexes
end

S_mul = repmat(S',8,1); % Eight years of data
figure
plot(data_s,'r')
hold on
grid on
plot(S_mul,'b')
legend('Detrended data', 'Seasonal indexes', 'location', 'best')

I = data_s-S_mul; % irregular parts
figure
grid on
hold on
plot(I)
legend('Irregular parts', 'location', 'best')

h = kpsstest(I) % 1, not stationary
I_diff = diff(I);
h_d = kpsstest(I_diff) % 0, now it's stationary :)

figure
autocorr(I_diff,40)
figure
parcorr(I_diff, 40)

models = [];
for i = 0:4 % Let's calculate 5*5 different models
    for j = 0:4
        model = armax(I_diff, [i j]);
        a = aic(model);
        models = [models; i j a];
    end
end

ind = find(models(:,end) == min(models(:,end)));
min = models(ind,:);

% Minimum aic value this time with 3 and 3.
ModelFit = armax(I_diff, [min(1) min(2)]);

I_df = forecast(ModelFit,I_diff,p); % forecast for _differences_!

figure
plot(I_diff,'b'), hold on, grid on
plot(length(I_diff)+1:length(I_diff)+p,I_df,'r')

legend('Differences', 'Forecasts of differences', 'location', 'best')

Xf = [ones(p,1) [length(data)+1:length(data)+p]'];
Tf = Xf*b; % Trend forecast

If(1) = I(end)+I_df(1); % Forecast of irregular series

for i = 2:p
    If(i) = If(i-1)+I_df(i);
end

forecast = Tf + S' + If'; 

figure % again, I had to do the same tricks with indexes....
plot(1:346,data_1(1:346),'b'), hold on, grid on % original data
plot(335:346,forecast,'r') % forecast 
legend('Original data', 'Forecast','location','best')
