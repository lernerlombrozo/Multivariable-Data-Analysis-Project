clear all;clc;

%% seminar 4 question 1 b
%%%Hotelings T^2 statistic

%%%%%%%% can be changed
X=[7 7;1 9;7 6;5 6];% 1b
mu=[7 33];% 1b

X=[5.526 59; 
    10.401 75; 
    9.213 69; 
    8.95 67.5; 
    7.063 62; 
    6.61 62; 
    11.273 74; 
    2.447 47; 
    5.149 59.5; 
    9.158 68; 
    9.004 69; 
    12.132 75; 
    8.199 70.5; 
    6.978 66.5; 
    6.601 64.5; 
    6.89 63; 
    7.622 67.5; 
    10.06 73; 
    10.09 73; 
    10.88 77; 
    7.6 61.5; 
    7.7 66.5; 
    12.01 79.5;
    10.04 74; 
    7 75];% 2a-e
mu=[5 60];% 2c
mu=[6 71];% 2d

alpha=0.05;
%%%%%%%%

[n,p]=size(X);
X_mean = mean(X);
X_cov= cov(X);

%Get hotellings t^2 distance
C_squared = (n-1)*p/(n-p)*finv(1-alpha,p,n-p);

%Get hotellings T^2 statistic to compare
T_squared = n*(X_mean-mu)*X_cov^-1*transpose(X_mean-mu);

if C_squared > T_squared
    disp('can be')
else
    disp('can NOT be')
end

%% seminar 4 question 1 c or 2 a
%determine the axes for the 1-alpha confidence elipse for mu. Determine the
%length of the axes
[e , L]=eig(X_cov);

%length of the axes from means along eigenvector ei
half_length=sqrt(L*C_squared/n);
%full length along eigenvector ei
full_length=2*half_length;

%% seminar 4 question 1 d or 2 b-e
%symultaneous confidence intervals
symultaneous_left_interval=transpose(X_mean)-sqrt(C_squared*diag(X_cov)/n);
symultaneous_right_interval=transpose(X_mean)+sqrt(C_squared*diag(X_cov)/n);
symultaneous_intervals=[symultaneous_left_interval symultaneous_right_interval];

% checking if mu fits in interval
(transpose(mu)> symultaneous_left_interval) & (transpose(mu) < symultaneous_right_interval)

%% seminar 4 question 2 e
%individual
%individual confidence intervals
individual_left_interval=transpose(X_mean)+tinv(alpha/2,n-1)*sqrt(diag(X_cov)/n);
individual_right_interval=transpose(X_mean)-tinv(alpha/2,n-1)*sqrt(diag(X_cov)/n);
individual_intervals=[individual_left_interval individual_right_interval];
% checking if mu fits in interval
(transpose(mu)> individual_left_interval) & (transpose(mu) < individual_right_interval)

%bonferroni
Bonferroni_left_interval=transpose(X_mean)+tinv(alpha/2/p,n-1)*sqrt(diag(X_cov)/n);
Bonferroni_right_interval=transpose(X_mean)-tinv(alpha/2/p,n-1)*sqrt(diag(X_cov)/n);
Bonferroni_intervals=[Bonferroni_left_interval Bonferroni_right_interval];
% checking if mu fits in interval
(transpose(mu)> Bonferroni_left_interval) & (transpose(mu) < Bonferroni_right_interval)

