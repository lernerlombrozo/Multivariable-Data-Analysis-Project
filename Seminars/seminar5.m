%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ajay Ganesh, Applied Multivariate Data Analysis CHE 494/694, Seminar 5
% Multivariate Regression Models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
clc;

filename= "Food_Data.xlsx";
FoodData = xlsread(filename,1,'B2:J26');
[~, XLabels, ~] = xlsread(filename,1,'B1:J1'); %Foods
[~, YLabels, ~] = xlsread(filename,1,'A2:A26'); %Countries

%%%%%%%%
x1columns= 5;
x2columns= 6;
ycolumns=9;
%%%%%%%%%

x1=FoodData(:,x1columns);
x2=FoodData(:,x2columns);
y=FoodData(:,ycolumns);

%load('LiquorMortality.mat');
% x1=Liquorconsumptionpercapita;
% x2=Wineconsumptionpercapita;
% y=Cirrhosisdeathrate;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data Visualization %%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
scatter(x1,x2);
xlabel(XLabels(x1columns)); 
ylabel(XLabels(x2columns));
title('2D-bivariate scatter plot');

figure;
subplot(2,2,1);
plot(1:length(x1),x1);
xlabel('sample'); ylabel(XLabels(x1columns));
subplot(2,2,2);
plot(1:length(x1),x2);
xlabel('sample'); ylabel(XLabels(x2columns));
subplot(2,2,3:4);
plot(1:length(x1),y);
xlabel('sample'); ylabel(XLabels(ycolumns));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Regression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = [ones(size(x1)) x1 x2 x1.*x2];
alpha=0.05; %Default

%regress treats NaN values in X or y as missing values. regress omits 
% observations with missing values from the regression fit.

[b,bint,r,rint,stats] = regress(y,X,alpha);

% mvregress for Multiple Multivariate Regression

figure;
scatter3(x1,x2,y,'filled')
hold on
x1fit = min(x1):0.1:max(x1);
x2fit = min(x2):0.1:max(x2);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT; %Regression Model
mesh(X1FIT,X2FIT,YFIT);
xlabel(XLabels(x1columns));
ylabel(XLabels(x2columns));
zlabel(XLabels(ycolumns));
view(50,10);
hold off

%Confidence intervals of the estimate, default alpha = 0.05
disp('Confidence intervals of the estimate');
disp(bint);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Residual Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(2,2,1);
scatter(r,x1);
xlabel(XLabels(x1columns)); 
ylabel('residue');
subplot(2,2,2);
scatter(r,x2);
xlabel(XLabels(x2columns)); 
ylabel('residue');
subplot(2,2,3);
plot(1:length(r),r);
xlabel('n'); ylabel('residue');
subplot(2,2,4);
scatter(r,y+r);
xlabel(sprintf('%s fitted',string(XLabels(ycolumns)))); 
ylabel('residue');
title('Residual analysis plots');

%Intervals to diagnose outliers, returned as a numeric matrix. rint is a 
%p-by-2 matrix, where p is the number of predictors in X. If the interval 
%rint(i,:) for observation i does not contain zero, the corresponding 
%residual is larger than expected in 100*(1-alpha)% of new observations, 
% suggesting an outlier.
disp('Intervals of the residue estimate');
disp(rint);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Outlier Detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Diagnose outliers by finding the residual intervals rint that do not contain 0.

contain0 = (rint(:,1)<0 & rint(:,2)>0);
idx = find(contain0==false);

figure;
hold on;
scatter(y,r);
scatter(y(idx),r(idx),'b','filled');
xlabel(XLabels(ycolumns));
ylabel('Residuals');
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%% Regression Statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%
%Model statistics, returned as a numeric vector including the R2 statistic, 
%the F-statistic and its p-value, and an estimate of the error variance
disp(stats);