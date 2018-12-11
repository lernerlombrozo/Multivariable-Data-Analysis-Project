clear all; clc; 

%% Loading the Data
filename= "NIR_Raman.xlsx";
%Labels were obtained from 
%https://github.com/cran/FastHCS/blob/master/data/Tablets.txt.gz

labels_file="NIR_Labels.xlsx";

%rows 1-310: tablets
%column 1-404: variables (peaks in the range 7400-10507 cm-1)
%column 405: w/w percentage of active substance in tablet
%column 406:production scale (0:laboratory, 1:pilot, 2:full)
%column 407:tablet type (0:A,1:B,2:C,3:D)
%%%%%%%%%%%%%%%%%
NIR_data = xlsread(filename,1);
NIR_variables=NIR_data(:,1:404);
NIR_labels=xlsread(labels_file,1);

%rows 1-120: tablets
%column 1-3401: variables (peaks in the range 3600-200 cm-1)
%column 3402: w/w percentage of active substance in tablet
%column 3403:tablet type (0:A,1:B,2:C,3:D)
%%%%%%%%%%%%%%%%%
Raman_data = xlsread(filename,2);
Raman_variables=Raman_data(:,1:3401);
Raman_labels = 200:3600;

NIR_mean=mean(NIR_variables);
Raman_mean=mean(Raman_variables);

%% visualizing the data
%%% In order to visualize the raw data, plots of the variables were done
%%% bar plots include an average of all the variables 
figure
plot(NIR_labels,NIR_variables)
xlabel('Wavenumber (cm^-1)');
ylabel('Absorbance');
title("NIR Spectroscopy samples")

figure
bar(NIR_labels,NIR_mean)
xlabel('Wavenumber (cm^-1)');
ylabel('Absorbance');
title("NIR Spectroscopy mean");


figure
plot(Raman_labels,Raman_variables)
xlabel('Raman Shift (cm^-1)');
ylabel('Intensity (arbitrary)');
title("Raman Spectroscopy samples");

figure
bar(Raman_labels,Raman_mean)
xlabel('Raman Shift cm^-1');
ylabel('Intensity (arbitrary)');
title("Raman Spectroscopy mean");

%% PCA
[NIR_coeff,NIR_score,NIR_latent,NIR_tsquared,NIR_explained,NIR_mu]=pca(NIR_data);
[Raman_coeff,Raman_score,Raman_latent,Raman_tsquared,Raman_explained,Raman_mu]=pca(Raman_data);

[e,L]=eig(cov(NIR_data));

%%% scree plot
figure;
subplot(1,2,1);
plot(NIR_latent(1:10));
title('NIR latent plot');
subplot(1,2,2);
plot(Raman_latent(1:10));
title('Raman latent plot');

%%% PCA Plot
dof=2; %p=2
alpha1=0.05; %95% confidence

NIR_pair=NIR_score(:,1:dof);
NIR_remainder=NIR_score(:,dof+1:end)

Raman_pair=Raman_score(:,1:dof);
Raman_remainder=Raman_score(:,dof+1:end)
%y1/L1 + y2/L2 =X^2(alpha)
%y2 =+-sqrt(L2(X^2(alpha)-y1/L1))
%y2 is valid while 
Chi=chi2inv(1-alpha1,dof);

%% NIR PCA Plot
NIR_xi=linspace(-sqrt(Chi*NIR_latent(1)),sqrt(Chi*NIR_latent(1)));
NIR_yp =sqrt(NIR_latent(2)*(Chi-(NIR_xi.^2./NIR_latent(1))));
NIR_yn =-sqrt(NIR_latent(2)*(Chi-(NIR_xi.^2./NIR_latent(1))));

outliers=[0 0];
outliersNumber=0;
outliersIndex=[0];
for v = 1:1:length(NIR_pair)
    if(NIR_pair(v,1)>sqrt(Chi*NIR_latent(1)) || NIR_pair(v,1)<-sqrt(Chi*NIR_latent(1)))
        outliersNumber = outliersNumber+1;
        outliers(outliersNumber,:)=NIR_pair(v,:);
        outliersIndex(outliersNumber,:)=v;
    elseif(NIR_pair(v,2)>sqrt(NIR_latent(2)*(Chi-(NIR_pair(v,1)^2./NIR_latent(1)))))
        outliersNumber = outliersNumber+1;
        outliers(outliersNumber,:)=NIR_pair(v,:);
        outliersIndex(outliersNumber,:)=v;
    elseif(NIR_pair(v,2)< -sqrt(NIR_latent(2)*(Chi-(NIR_pair(v,1)^2./NIR_latent(1)))))
        outliersNumber = outliersNumber+1;
        outliers(outliersNumber,:)=NIR_pair(v,:);
        outliersIndex(outliersNumber,:)=v;
    end
end

figure;
scatter(NIR_pair(:,1),NIR_pair(:,2))
xlabel('1st PC'); ylabel('2nd PC');
title('NIR Principal Components');
hold on
if outliers ~= [0 0]
    scatter(outliers(:,1),outliers(:,2),'MarkerEdgeColor',[0 .5 .5],...
                  'MarkerFaceColor',[0 .7 .7],...
                  'LineWidth',1.5)
end
plot(NIR_xi,NIR_yp,'red')
plot(NIR_xi,NIR_yn,'red')
legend('scores','potential outliers','alpha=0.05 (95% confidence)','alpha=0.05 (negatives)','Location','northwest')
hold off

%% Raman PCA Plot
Raman_xi=linspace(-sqrt(Chi*Raman_latent(1)),sqrt(Chi*Raman_latent(1)));
Raman_yp =sqrt(Raman_latent(2)*(Chi-(Raman_xi.^2./Raman_latent(1))));
Raman_yn =-sqrt(Raman_latent(2)*(Chi-(Raman_xi.^2./Raman_latent(1))));

outliers=[0 0];
outliersNumber=0;
outliersIndex=[0];
for v = 1:1:length(Raman_pair)
    if(Raman_pair(v,1)>sqrt(Chi*Raman_latent(1)) || Raman_pair(v,1)<-sqrt(Chi*Raman_latent(1)))
        outliersNumber = outliersNumber+1;
        outliers(outliersNumber,:)=Raman_pair(v,:);
        outliersIndex(outliersNumber,:)=v;
    elseif(Raman_pair(v,2)>sqrt(Raman_latent(2)*(Chi-(Raman_pair(v,1)^2./Raman_latent(1)))))
        outliersNumber = outliersNumber+1;
        outliers(outliersNumber,:)=Raman_pair(v,:);
        outliersIndex(outliersNumber,:)=v;
    elseif(Raman_pair(v,2)< -sqrt(Raman_latent(2)*(Chi-(Raman_pair(v,1)^2./Raman_latent(1)))))
        outliersNumber = outliersNumber+1;
        outliers(outliersNumber,:)=Raman_pair(v,:);
        outliersIndex(outliersNumber,:)=v;
    end
end

figure;
scatter(Raman_pair(:,1),Raman_pair(:,2))
xlabel('1st PC'); ylabel('2nd PC');
title('Raman Principal Components');
hold on
if outliers ~= [0 0]
    scatter(outliers(:,1),outliers(:,2),'MarkerEdgeColor',[0 .5 .5],...
                  'MarkerFaceColor',[0 .7 .7],...
                  'LineWidth',1.5)
end
plot(Raman_xi,Raman_yp,'red')
plot(Raman_xi,Raman_yn,'red')
legend('scores','potential outliers','alpha=0.05 (95% confidence)','alpha=0.05 (negatives)','Location','northwest')
hold off

%% hold for now
figure;
qqplot(X_pair(:,1),X_pair(:,2))
xlabel('x1'); ylabel('x2');
title('Q-Q plot');

Chi3=chi2inv(1-alpha1,3);
xUCL=linspace(0, length(X_remainder)+5);
yUCL=xUCL.*0.+Chi3;
Tsquaredcols=(X_remainder.^2)./[eva(3,3) eva(2,2) eva(1,1)]
Tsquared=sum(Tsquaredcols,2)

newOutliersNumber=0;
newOutliersIndex=[0];
for v = 1:size(Tsquared)
    if (Tsquared(v)>Chi3)
        newOutliersNumber=newOutliersNumber+1;
        newOutliersIndex(newOutliersNumber,:)=v;
    end    
end


b= bar(Tsquared, 'FaceColor','flat');
for v = 1:size(outliersIndex)
    b.CData(outliersIndex(v),:) = [0 .7 .7];
end
for v = 1:size(newOutliersIndex)
    b.CData(newOutliersIndex(v),:) = [0 .5 .2];
end
hold on
plot(xUCL,yUCL,'red');
xlabel('Period'); ylabel('T^2');
hold off


%%%Todo:
%%%PCA
%%%PLS
%%%PCR? 
%%%PLS?
%%%LDA linear discriminant analysis (from machine learning course)
%%%Multiplicative scatter cor- rection (MSC)
%%%First and second derivatives
%%%standard normal variate (SNV)
%%%interval PLS (iPLS) with optimization, a systematic method for vari- able selection developed by Nørgaard and co-workers.14
