%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Ajay Ganesh, Applied Multivariate Data Analysis CHE 494/694
% Assignment 2 Solution 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% getting data
close all;
clear all; 
clc;
load('DataSet.mat'); 
x1=DataSet(:,1); 
x2=DataSet(:,2); 
x3=DataSet(:,3); 
y1=DataSet(:,4); 
y2=DataSet(:,5);

%% Data Visualization
figure;
scatter3(x1,x2,x3);
xlabel('x1'); 
ylabel('x2'); 
zlabel('x3'); 
title('3D scatter plot');
figure; 
Varnum=size(DataSet,2); 
for i=1:Varnum
    for j=1:Varnum
        subplot(Varnum,Varnum,((i-1)*Varnum)+j); 
        scatter(DataSet(:,i),DataSet(:,j)); 
        xlabel(sprintf('Column%d',i)); 
        ylabel(sprintf('Column%d',j));
    end
end

%% regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Regression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = [ones(size(x1)) x1 x2 x3 x1.*x2.*x3];
alpha=0.05; %Default

%regress treats NaN values in X or y as missing values. regress omits % observations with missing values from the regression fit.
[b1,bint1,r1,rint1,stats1] = regress(y1,X,alpha); 
[b2,bint2,r2,rint2,stats2] = regress(y2,X,alpha);

% mvregress for Multiple Multivariate Regression
%Confidence intervals of the estimate, default alpha = 0.05
disp('Confidence intervals of the estimate for b1');
disp(bint1);
disp('Confidence intervals of the estimate for b2');
disp(bint2);

%% Residual Analysis
figure;
subplot(2,2,1); scatter(r1,x1);
xlabel('x1'); ylabel('residue'); 
subplot(2,2,2); scatter(r1,x2);
xlabel('x2'); ylabel('residue'); 
subplot(2,2,3); scatter(r1,x3);
xlabel('x3'); ylabel('residue');
figure;
subplot(2,1,1);
plot(1:length(r1),r1);
xlabel('n'); ylabel('residue'); 
subplot(2,1,2);
scatter(r1,y1+r1);
xlabel('y1 fitted'); 
ylabel('residue'); 
title('Residual analysis plots for y2');

figure;
subplot(2,2,1); scatter(r2,x1);
xlabel('x1'); ylabel('residue'); subplot(2,2,2); scatter(r2,x2);
xlabel('x2'); ylabel('residue'); subplot(2,2,3); scatter(r2,x3);
xlabel('x3'); ylabel('residue');
figure;
subplot(2,1,1);
plot(1:length(r2),r2);
xlabel('n'); 
ylabel('residue'); 
subplot(2,1,2);
scatter(r2,y2+r2);
xlabel('y2 fitted'); 
ylabel('residue'); 
title('Residual analysis plots for y2');

%Intervals to diagnose outliers, returned as a numeric matrix. rint is a %p-by-2 matrix, where p is the number of predictors in X. If the interval %rint(i,:) for observation i does not contain zero, the corresponding %residual is larger than expected in 100*(1-alpha)% of new observations, % suggesting an outlier.
disp('Intervals of the residue estimate for output 1');
disp(rint1);
disp('Intervals of the residue estimate for output 2');
disp(rint2);

%%%%%%%%%%%%%%%%%%%%%%%%% Regression Statistics %%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Model statistics, returned as a numeric vector including the R2 statistic,
%the F-statistic and its p-value, and an estimate of the error variance
disp('Statistics for output 1');
disp(stats1);
disp('Statistics for output 2');
disp(stats2);

%% PCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PCA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
X=[x1 x2 x3];
Y=[y1 y2];
[X_coeff,X_score,X_latent]=pca(X);
[Y_coeff,Y_score,Y_latent]=pca(Y);
figure;
subplot(1,2,1);
plot(X_latent);
title('X latent plot');
subplot(1,2,2);
plot(Y_latent);
title('Y latent plot');
disp('Both X and Y seems to be reducable as we can see a "knee" in both of the plots'); 

%% PCR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PCR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
x1_pcr=X_score(:,1);
x2_pcr=X_score(:,2);
x3_pcr=X_score(:,3);

X_pcr = [ones(size(x1)) x1_pcr x2_pcr x3_pcr x1_pcr.*x2_pcr.*x3_pcr]; 
alpha=0.05; %Default

%regress treats NaN values in X or y as missing values. regress omits 
% observations with missing values from the regression fit.

y1_centr=y1-mean(y1); %Mean centering the data 
y2_centr=y2-mean(y2);

[b1_pcr,bint1_pcr,r1_pcr,rint1_pcr,stats1_pcr] = regress(y1_centr,X_pcr,alpha);

[b2_pcr,bint2_pcr,r2_pcr,rint2_pcr,stats2_pcr] = regress(y2_centr,X_pcr,alpha);

% mvregress for Multiple Multivariate Regression 
%Confidence intervals of the estimate, default alpha = 0.05
disp('Confidence intervals of the estimate for b1');
disp(bint1_pcr);
disp('Confidence intervals of the estimate for b2');
disp(bint2_pcr);

%% PLS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
ncomp=2;
x1_centr=x1-mean(x1); 
x2_centr=x2-mean(x2); 
x3_centr=x3-mean(x3);

X_centr=[x1_centr x2_centr x3_centr]; 
Y_centr=[y1_centr y2_centr];

[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats_pls] = plsregress(X_centr,Y_centr,ncomp);
