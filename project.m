%clear all; clc; 

%%% Loading the Data
filename= "NIR_Raman.xlsx";
%Labels were obtained from 
%https://github.com/cran/FastHCS/blob/master/data/Tablets.txt.gz

labelsFile="NIR_Labels.xlsx";

%rows 1-310: tablets
%column 1-404: variables (peaks in the range 7400-10507 cm-1)
%column 405: w/w percentage of active substance in tablet
%column 406:production scale (0:laboratory, 1:pilot, 2:full)
%column 407:tablet type (0:A,1:B,2:C,3:D)
%NIRData = xlsread(filename,1);
%NIRVariables=NIRData(:,1:404);
%NIRLabels=xlsread(labelsFile,1);

%rows 1-120: tablets
%column 1-3401: variables (peaks in the range 3600-200 cm-1)
%column 3402: w/w percentage of active substance in tablet
%column 3403:tablet type (0:A,1:B,2:C,3:D)
%RamanData = xlsread(filename,2);
RamanVariables=RamanData(:,1:3401);
RamanLabels = 200:3600;

NIRTotals=mean(NIRVariables);
RamanTotals=mean(RamanVariables);

%%% In order to visualize the raw data, plots of the variables were done
%%% bar plots include an average of all the variables 
figure
plot(NIRLabels,NIRVariables)
figure
bar(NIRLabels,NIRTotals)
figure
plot(RamanLabels,RamanVariables)
figure
bar(RamanLabels,RamanTotals)


