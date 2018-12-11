clear all; clc; 

%%% Loading the Data
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

%%%Todo:
%%%PLS
%%%Multiplicative scatter cor- rection (MSC)
%%%First and second derivatives
%%%standard normal variate (SNV)
%%%interval PLS (iPLS) with optimization, a systematic method for vari- able selection developed by Nørgaard and co-workers.14
