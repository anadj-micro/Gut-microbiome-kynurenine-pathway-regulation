%% Script to plot a heatmap with genus abundances in samples from probiotic
% experiments to show that the microbiome composition is similar between
% probiotic responders and non-responders
close all
clear all

% read genus table
asv = readtable('genus_percentage_dada2.txt', 'Delimiter', '\t');

%% add info about ido1 expression

ido = readtable('ido_tdo_levels.txt');
% join with microbiome data
asv1 = innerjoin(asv, ido);
% rearrange
asv1 = asv1(:, [1 end-7:end 2:end-8]);

%% select only top 100 genera
T1 = 100;

abundanceData = asv1{:, 10:end};
uTaxa = asv1.Properties.VariableNames(10:end);
abundanceMean =  mean(abundanceData, 1);

% sorting taxa by their abundance mean
[abundanceMean1, iSort] = sort(abundanceMean, 'descend');
abundanceMatrix = abundanceData(:, iSort);
uTaxa = uTaxa(iSort);

% consolidate 'Other' if abundant taxa number is lower than uTaxa number. 
if length(uTaxa) > T1 
    abundanceMatrix1 = abundanceMatrix(:, 1:T1);
    uTaxa1 = uTaxa(1:T1);
    rowU = abundanceMatrix( :, (T1+1):end);
    other = sum(rowU,2);
    %test if it fits expected
    sum(abundanceMatrix1,2) + other
    abundanceMatrix1 = [abundanceMatrix1 other];
    uTaxa1 =[uTaxa1 'Other'];
end
asv2 = [asv1(:,1:9) array2table(abundanceMatrix1)];
asv2.Properties.VariableNames(10:end) = uTaxa1;

%% Sort the table based on sample group
asv2 = sortrows(asv2, 'SampleGroup');

% change values with following rules:
%          C1(0)    < 0.001
% 0.001 =< C2(0.001)< 0.01
% 0.01  =< C3(0.01) < 0.1
% 0.1   =< C4(0.1)  < 1
% 1     =< C5(1)    < 10
% 10    =< C6(10)

abundanceData = asv2{:, 10:end};

mask1 = abundanceData >= 10;
abundanceData(mask1) = 10;
mask2 = abundanceData < 10 & abundanceData >= 1;
abundanceData(mask2) = 8;
mask3 = abundanceData < 1 & abundanceData >= 0.1;
abundanceData(mask3) = 6;
mask4 = abundanceData < 0.1 & abundanceData >= 0.01;
abundanceData(mask4) = 4;
mask5 = abundanceData < 0.01 & abundanceData >= 0.001;
abundanceData(mask5) = 2;
mask6 = abundanceData < 0.001;
abundanceData(mask6) = 0;

asv3 = [asv2(:,1:9) array2table(abundanceData)];
asv3.Properties.VariableNames(10:end) = uTaxa1;

% sort taxa based on the abundances in Control group
abundanceData1 = asv3{strcmp(asv3.SampleGroup, "C") == 1,10:end};
abundanceData2 = asv3{:,10:end};
uTaxa2 = asv3.Properties.VariableNames(10:end);
abundanceMean =  mean(abundanceData1, 1);

% sorting taxa by their abundance mean
[abundanceMean1, iSort] = sort(abundanceMean, 'descend');
abundanceMatrix1 = abundanceData2(:, iSort);
uTaxa3 = uTaxa2(iSort);

asv4 = [asv3(:,1:9) array2table(abundanceMatrix1)];
asv4.Properties.VariableNames(10:end) = uTaxa3;

% remove "other" and put probiotics at the end
asv4 = removevars(asv4, 'Other');
asv4 = asv4(:, [1:18 20:27 29:96 99:101 103:end 19 28 97 98 102]);

%% plot

% set colormap
map = [255 255 255;
    213 228 245;
    131 175 225;
    48 122 206;
    29 73 123;
    9 24 41];
map1 = map./255;
cmap = colormap(map1);

xl = asv4.Properties.VariableNames(10:end);

subplot(4,1,1)
ylC = asv4{strcmp(asv4.SampleGroup, "C") == 1, 'Sample_ID'};
h = heatmap(asv4{strcmp(asv4.SampleGroup, "C") == 1,10:end}, 'Colormap', cmap)
set(h, "YLabel", 'Untreated')
h.XDisplayLabels = nan(size(h.XDisplayData));
h.YDisplayLabels = nan(size(h.YDisplayData));

subplot(4,1,2)
ylAVN = asv4{strcmp(asv4.SampleGroup, "AVN") == 1, 'Sample_ID'};
h = heatmap(asv4{strcmp(asv4.SampleGroup, "AVN") == 1,10:end}, 'Colormap', cmap)
set(h, "YLabel", 'AVN')
h.XDisplayLabels = nan(size(h.XDisplayData));
h.YDisplayLabels = nan(size(h.YDisplayData));

subplot(4,1,3)
ylAVNNR = asv4{strcmp(asv4.SampleGroup, "AVN_P_NR") == 1, 'Sample_ID'};
h = heatmap(asv4{strcmp(asv4.SampleGroup, "AVN_P_NR") == 1,10:end}, 'Colormap', cmap)
set(h, "YLabel", 'Non-responders')
h.XDisplayLabels = nan(size(h.XDisplayData));
h.YDisplayLabels = nan(size(h.YDisplayData));

subplot(4,1,4)
ylAVNR = asv4{strcmp(asv4.SampleGroup, "AVN_P_R") == 1, 'Sample_ID'};
h = heatmap(xl, ylAVNR, asv4{strcmp(asv4.SampleGroup, "AVN_P_R") == 1,10:end}, 'Colormap', cmap)
set(h, "YLabel", 'Responders')
h.XDisplayLabels = nan(size(h.XDisplayData));
h.YDisplayLabels = nan(size(h.YDisplayData));

print -painters -dpdf 10_16_24.pdf


