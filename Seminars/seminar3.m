clear all;clc;

%% seminar 3 question 1B
sigma11 = 2;
sigma22 = 4;
delta12 = -0.8;

x = chi2inv(0.5,2);

variance = [sigma11 delta12*sqrt(sigma11*sigma22); delta12*sqrt(sigma11*sigma22) sigma22];
[e , L]=eig(variance);

L=sum(L);

halfLengthe1=sqrt(max(L)*x);
halfLengthe2=sqrt(min(L)*x);

%% seminar 3 question 3
%% Loading Samples
filename= "Food_Data.xlsx";
FoodData = xlsread(filename,1,'B2:J26');
[~, XLabels, ~] = xlsread(filename,1,'B1:J1'); %Foods
[~, YLabels, ~] = xlsread(filename,1,'A2:A26'); %Countries
[Xn, Yn] = size(FoodData);

%% a)
X_mean = mean(FoodData);
X_covariance = cov(FoodData);
X_correlation = corr(FoodData);

%% b)

columns = 3; %change to get more columns/rows
%col = ceil(rand(1,columns)*Yn); % random (might repeat)
col = 1:1:columns; % first
% 2D Scatters
figure
count=0;
for xcol = 1:1:columns
    for ycol =1:1:columns
        count=count+1;
        subplot(columns,columns,count);
        scatter(FoodData(:,col(xcol)),FoodData(:,col(ycol)))
        xlabel(XLabels(col(xcol)))
        ylabel(XLabels(col(ycol)))
    end
end

%3D Scatter
figure
scatter3(FoodData(:,col(1)),FoodData(:,col(2)),FoodData(:,col(3)))
xlabel(XLabels(col(1)))
ylabel(XLabels(col(2)))
zlabel(XLabels(col(3)))
title('3D-trivariate scatter plot');

%% c)
%%Bivariate Analaysis
%%Fitting a bivariate gaussian PDF to the data
columns=[2 4]; %only two values
X_pair=[FoodData(:,columns(1)) FoodData(:,columns(2))];
figure; 
scatter(X_pair(:,1),X_pair(:,2)) 
xlabel(XLabels(columns(1))); 
ylabel(XLabels(columns(2))); 
title('Bivariate scatter plot');

X_pair_mean=mean(X_pair) 
X_pair_cov=cov(X_pair)
k=0.5;
x1 = -k*X_pair_cov(1,1):.2:k*X_pair_cov(1,1); 
x2 = -k*X_pair_cov(2,2):.2:k*X_pair_cov(2,2);
x1=x1+X_pair_mean(1); x2=x2+X_pair_mean(2);
[X1,X2] = meshgrid(x1,x2);
F_x1x2 = mvnpdf([X1(:) X2(:)],X_pair_mean, X_pair_cov); 
F_x1x2 = reshape(F_x1x2,length(x2),length(x1));
figure;
surf(x1,x2,F_x1x2);
xlabel(XLabels(columns(1))); ylabel(XLabels(columns(2))); zlabel('Probability Density');
title(sprintf('Bivariate distribution of [%s %s] with mu=[%4.2f %4.2f] and sigma=[%4.2f %4.2f;%4.2f %4.2f]',string(XLabels(columns(1))),string(XLabels(columns(2))),X_pair_mean(1),X_pair_mean(2),X_pair_cov(1,1),X_pair_cov(1,2),X_pair_cov(2,1),X_pair_cov(2,2)));

%% d)
dof=length(columns); %p=2
alpha=[0.05 0.5 0.2]; %95% and 50% confidence.

contour_value=exp(-0.5*chi2inv(1-alpha,dof))/(2*pi*sqrt(det(X_pair_cov)));
figure;
scatter(X_pair(:,1),X_pair(:,2));
xlim([x1(1) x1(end)]);
ylim([x2(1) x2(end)]);
hold on;
contour(x1,x2,F_x1x2,[contour_value contour_value]);
%contour(x1,x2,F_x1x2,[contour_value(2) contour_value(2)],'color','b');
legend('Data Points','\alpha=0.05 0.5 (95 50% confidence)','Location','northoutside','Orientation','vertical')
xlabel('x1');
ylabel('x2');
title(sprintf('Confidence ellipse for alpha=%4.2f and alpha=%4.2f',alpha(1),alpha(2)));
hold off;
