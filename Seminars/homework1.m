clear all;clc;

%% homework 1 question 1 c
%%%%%%%%%
m=[2; -3; 4];
Co = [1 1 1; 1 3 2; 1 2 2];
alpha = 0.05;
%%%%%%%%%

p=length(m);

xn = m(1:2,1:1);
x3=m(3:3,1:1);
syms x

Z11= Co(1:2,1:2);
Z12= Co(1:2,3:3);
Z21= Co(3:3,1:2);
Z22= Co(3:3,3:3);

newMean= xn+ Z12 * Z22^-1 * (x - x3)
newZ = Z11 - Z12*Z22^-1*Z21

%% question 1 d
Chi_squared=chi2inv(1-alpha,p)

%% question 2 a
%% Loading Samples
filename= "Food_Data.xlsx";
FoodData = xlsread(filename,1,'B2:J26');
[~, XLabels, ~] = xlsread(filename,1,'B1:J1'); %Foods
[~, YLabels, ~] = xlsread(filename,1,'A2:A26'); %Countries
[Xn, Yn] = size(FoodData);

%% a)

%%%%%%%%%
smallData=[FoodData(:,2) FoodData(:,4) FoodData(:,6)];
%%%%%%%%%

smallData_mean = mean(smallData)
smallData_cov = cov(smallData)
smallData_corr = corr(smallData)

[n,p]=size(smallData);

%%%%%%%%%
col = [2,4,6]; 
columns = length(col);
%%%%%%%%%

% 2D Scatters
figure
count=0;
for xcol = 1:1:columns
    for ycol = 1:1:columns
        count=count+1;
        subplot(columns,columns,count);
        scatter(FoodData(:,col(ycol)),FoodData(:,col(xcol)))
        xlabel(XLabels(col(ycol)))
        ylabel(XLabels(col(xcol)))
    end
end

%3D Scatter
figure
scatter3(FoodData(:,col(1)),FoodData(:,col(2)),FoodData(:,col(3)))
xlabel(XLabels(col(1)))
ylabel(XLabels(col(2)))
zlabel(XLabels(col(3)))
title('3D-trivariate scatter plot');


%% b)
%determine the axes for the 1-alpha confidence elipse for mu. Determine the
%length of the axes
%Get hotellings t^2 distance
C_squared = (n-1)*p/(n-p)*finv(1-alpha,p,n-p);

[e , L]=eig(smallData_cov);
%length of the axes from means along eigenvector ei
half_length=sqrt(L*C_squared/n);
%full length along eigenvector ei
full_length=2*half_length;

%% c)
%%%%%%%%%
mu=[8 17 31];
alpha=0.05;
%%%%%%%%%

C_squared = (n-1)*p/(n-p)*finv(1-alpha,p,n-p);
Fstatistic=n*((smallData_mean-mu) * inv(smallData_cov) *transpose(smallData_mean-mu))
if Fstatistic <= C_squared
    disp('As Fstatistic <= C_squared, mu=mu0 null hypothesis is accepted'); 
else
    disp('As Fstatistic >= C_squared mu=mu0 null hypothesis is rejected'); 
end

%% d-e)
%individual confidence intervals
individual_left_interval=transpose(smallData_mean)+tinv(alpha/2,n-1)*sqrt(diag(smallData_cov)/n);
individual_right_interval=transpose(smallData_mean)-tinv(alpha/2,n-1)*sqrt(diag(smallData_cov)/n);
individual_intervals=[individual_left_interval individual_right_interval];
% checking if mu fits in interval
(transpose(mu)> individual_left_interval) & (transpose(mu) < individual_right_interval)

%symultaneous confidence intervals
symultaneous_left_interval=transpose(smallData_mean)-sqrt(C_squared*diag(smallData_cov)/n);
symultaneous_right_interval=transpose(smallData_mean)+sqrt(C_squared*diag(smallData_cov)/n);
symultaneous_intervals=[symultaneous_left_interval symultaneous_right_interval];
% checking if mu fits in interval
(transpose(mu)> symultaneous_left_interval) & (transpose(mu) < symultaneous_right_interval)

%Bonferroni confidence intervals
Bonferroni_left_interval=transpose(smallData_mean)+tinv(alpha/2/p,n-1)*sqrt(diag(smallData_cov)/n);
Bonferroni_right_interval=transpose(smallData_mean)-tinv(alpha/2/p,n-1)*sqrt(diag(smallData_cov)/n);
Bonferroni_intervals=[Bonferroni_left_interval Bonferroni_right_interval];
% checking if mu fits in interval
(transpose(mu)> Bonferroni_left_interval) & (transpose(mu) < Bonferroni_right_interval)
