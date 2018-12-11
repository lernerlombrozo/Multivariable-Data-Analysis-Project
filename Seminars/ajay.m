%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % Ajay Ganesh, Applied Multivariate Data Analysis CHE 494/694,
% Assignment 1 Solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% close all;
clear all; clc;
filename= "Food_Data.xlsx";
FoodData = xlsread(filename,1,'B2:J26');
X=[FoodData(:,2) FoodData(:,4) FoodData(:,6)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('------------------------------------------------------------------'); 
fprintf('\nPart A solution\n'); disp('------------------------------------------------------------------'); disp('Sample mean');
X_mean=mean(X); X_mean=X_mean';
disp('Sample variance'); X_cov=cov(X)
disp('Sample correlation'); X_corr=corr(X)
figure; scatter3(X(:,1),X(:,2),X(:,3)); xlabel('x1'); ylabel('x2'); zlabel('x3'); title('3D-trivariate scatter plot');
figure;
subplot(3,3,1);
scatter(X(:,1),X(:,1));
xlabel('x1'); ylabel('x1');
subplot(3,3,2);
scatter(X(:,2),X(:,1));
xlabel('x2'); ylabel('x1');
subplot(3,3,3);
scatter(X(:,3),X(:,1));
xlabel('x3'); ylabel('x1');
subplot(3,3,4);
scatter(X(:,1),X(:,2));
xlabel('x1'); ylabel('x2');
subplot(3,3,5);
scatter(X(:,2),X(:,2));
xlabel('x2'); ylabel('x2');
subplot(3,3,6);
scatter(X(:,3),X(:,2));
xlabel('x3'); ylabel('x2');
subplot(3,3,7);
scatter(X(:,1),X(:,3));
xlabel('x1'); ylabel('x3');
subplot(3,3,8);
scatter(X(:,2),X(:,3));
xlabel('x2'); ylabel('x3');
subplot(3,3,9);
scatter(X(:,3),X(:,3));
xlabel('x3'); ylabel('x3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part B %%%%%%%%%%%%%%%%%%%%%%%%%%%%% disp('------------------------------------------------------------------');
fprintf('\nPart B solution\n');
disp('------------------------------------------------------------------');
[evtr,eval]=eig(X_cov);
alpha_b=0.1;
p=3;
n=25;
cchi=chi2inv(1-alpha_b,p);
cf=(n-1)*p/(n-p)*finv(1-alpha_b,p,n-p);
MajorAxisLength=cf*sqrt(eval(3,3))/sqrt(n);
fprintf('MajorAxisLength=%4.2f along [%4.2f %4.2f %4.2f]\n\n',MajorAxisLength,evtr(:,3)); MiddleAxisLength=cf*sqrt(eval(2,2))/sqrt(n);
fprintf('MiddleAxisLength=%4.2f along [%4.2f %4.2f %4.2f]\n\n',MiddleAxisLength,evtr(:,2)); MinorAxisLength=cf*sqrt(eval(1,1))/sqrt(n);
fprintf('MinoAxisLength=%4.2f along [%4.2f %4.2f %4.2f]\n\n',MinorAxisLength,evtr(:,1)); %%%%%%%%%%%%%%%%%%%%%%%%%%%% Part C %%%%%%%%%%%%%%%%%%%%%%%%%%%%% disp('------------------------------------------------------------------');
fprintf('\nPart C solution\n');
disp('------------------------------------------------------------------');
alpha_c=0.2;
Threshold=(n-1)*p/(n-p)*finv(1-alpha_c,p,n-p);
mu0=[8;17;31];
Fstatistic=n*(X_mean-mu0)'*inv(X_cov)*(X_mean-mu0);
fprintf('\nThe Fstatistic is %4.2f and the threshold is %4.2f\n',Fstatistic,Threshold);
if Fstatistic <= Threshold
disp('As Fstatistic <= Threshold, mu=mu0 null hypothesis is accepted'); else
disp('As Fstatistic >= Threshold mu=mu0 null hypothesis is rejected'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('------------------------------------------------------------------'); 
fprintf('\nPart D solution\n'); disp('------------------------------------------------------------------'); alpha_d=0.05;
a=[1;0;0];
x1i=[a'*X_mean+tinv(alpha_d/2,n-1)*sqrt(a'*X_cov*a)/sqrt(n) a'*X_mean-tinv(alpha_d/2,n-1)*sqrt(a'*X_cov*a)/sqrt(n)]; 
fprintf('Individual confidence inteval for x1=[%4.2f %4.2f]\n\n',x1i);
a=[0;1;0];
x2i=[a'*X_mean+tinv(alpha_d/2,n-1)*sqrt(a'*X_cov*a)/sqrt(n) a'*X_mean-tinv(alpha_d/2,n-1)*sqrt(a'*X_cov*a)/sqrt(n)]; fprintf('Individual confidence inteval for x2=[%4.2f %4.2f]\n\n',x2i);
a=[0;0;1];
x3i=[a'*X_mean+tinv(alpha_d/2,n-1)*sqrt(a'*X_cov*a)/sqrt(n) a'*X_mean-tinv(alpha_d/2,n-1)*sqrt(a'*X_cov*a)/sqrt(n)]; fprintf('Individual confidence inteval for x3=[%4.2f %4.2f]\n\n',x3i);
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part E %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('------------------------------------------------------------------'); fprintf('\nPart E solution\n'); disp('------------------------------------------------------------------'); alpha_e=0.05;
a=[1;0;0];
x1s=[a'*X_mean-sqrt((n-1)*p/(n-p)*finv(1-alpha_e,p,n-p)*(a'*X_cov*a)/n) a'*X_mean+sqrt((n-1)*p/(n-p)*finv(1-alpha_e,p,n-p)*(a'*X_cov*a)/n)]; fprintf('Simultaneous confidence inteval for x1=[%4.2f %4.2f]\n\n',x1s);
a=[0;1;0];
x2s=[a'*X_mean-sqrt((n-1)*p/(n-p)*finv(1-alpha_e,p,n-p)*(a'*X_cov*a)/n) a'*X_mean+sqrt((n-1)*p/(n-p)*finv(1-alpha_e,p,n-p)*(a'*X_cov*a)/n)]; fprintf('Simultaneous confidence inteval for x2=[%4.2f %4.2f]\n\n',x2s);
a=[0;0;1];
x3s=[a'*X_mean-sqrt((n-1)*p/(n-p)*finv(1-alpha_e,p,n-p)*(a'*X_cov*a)/n) a'*X_mean+sqrt((n-1)*p/(n-p)*finv(1-alpha_e,p,n-p)*(a'*X_cov*a)/n)]; fprintf('Simultaneous confidence inteval for x3=[%4.2f %4.2f]\n\n',x3s);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=[1;0;0];
x1b=[a'*X_mean+tinv(alpha_e/(2*p),n-1)*sqrt(a'*X_cov*a)/sqrt(n) a'*X_mean-tinv(alpha_e/(2*p),n-1)*sqrt(a'*X_cov*a)/sqrt(n)]; fprintf('Bonferroni confidence inteval for x1=[%4.2f %4.2f]\n\n',x1b);
a=[0;1;0];
x2b=[a'*X_mean+tinv(alpha_e/(2*p),n-1)*sqrt(a'*X_cov*a)/sqrt(n) a'*X_mean-tinv(alpha_e/(2*p),n-1)*sqrt(a'*X_cov*a)/sqrt(n)]; fprintf('Bonferroni confidence inteval for x2=[%4.2f %4.2f]\n\n',x2b);
a=[0;0;1];
x3b=[a'*X_mean+tinv(alpha_e/(2*p),n-1)*sqrt(a'*X_cov*a)/sqrt(n) a'*X_mean-tinv(alpha_e/(2*p),n-1)*sqrt(a'*X_cov*a)/sqrt(n)]; fprintf('Bonferroni confidence inteval for x3=[%4.2f %4.2f]\n\n',x3b);
disp('As expected Bonferroni performs the best than simultaneous'); %%%%%%%%%%%%