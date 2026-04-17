%% Script for plotting stacked barplots: one bar per sample, different colors
% correspond to different ASVs- phylogenetically close ASVs similar colors.
% first run main.m to assgin colors to ASVs and get tblASVtaxonomy_new.csv
clear all;
close all;

% read abundance table
genusPercentage = readtable('ASV_unstacked_percentage.txt');

%% Reading ido1 table and merging
ido = readtable('ido_tdo_levels', 'Delimiter', '\t');
genusPercentage1 = innerjoin(genusPercentage, ido);
% reorganize columns
genusPercentage1 = genusPercentage1(:,[1 end-7:end 2:end-8]);

% sort samples based on ido1 expression
genusPercentage1 = sortrows(genusPercentage1,'ido1','ascend');

%%
% now select based on the abundance threshold of 0.1%
% setup threshold for prevalence
T1 = 500;

abundanceData = genusPercentage1{:, 10:end};
uTaxa = genusPercentage1.Properties.VariableNames(10:end);
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
genusPercentage2 = genusPercentage1;
genusPercentage2 = [genusPercentage1(:,1:9) array2table(abundanceMatrix1)];
genusPercentage2.Properties.VariableNames(10:end) = uTaxa1;

T3Interest1 = genusPercentage2;

%% Add information about 16S counts
% read 16S table
BLoad = readtable('16S_qPCR', 'Delimiter', '\t');
BLoad.Properties.VariableNames{1} = 'Sample_ID';
% merge
genusPercentage3 = innerjoin(genusPercentage2, BLoad(:,[1 end]));
genusPercentage3 = genusPercentage3(:, [1:9 end 10:end-1]);

%% Add information about bacterial diversity
% read diversity table
diversity = readtable('Diversity.txt', 'Delimiter', '\t');
% merge
genusPercentage3 = innerjoin(genusPercentage3, diversity(:,[1 end]));
genusPercentage3 = genusPercentage3(:, [1:10 end 11:end-1]);
% sort samples based on ido1 expression
genusPercentage3 = sortrows(genusPercentage3,'ido1','ascend');

T3Interest1 = genusPercentage3;

%% Preparing for plotting
T3Interest1 = T3Interest1(:, [1 12:end]);

% get the taxonomy table and subset only used asvs
taxonomy = readtable('tblASVtaxonomy_new.csv');
uTaxa = T3Interest1.Properties.VariableNames(2:end);

subtaxonomy1 = taxonomy(ismember(taxonomy.ASV, uTaxa(1:end-1)) == 1,:);
% sort by taxonomy
subtaxonomy1 = sortrows(subtaxonomy1, 'ColorOrder');
% sort table based on taxonomy
sortedNames = subtaxonomy1{:,1};
T3Interest1 = [T3Interest1(:,1) T3Interest1(:,sortedNames) T3Interest1(:,end)];

%plotting

map = hex2rgb(subtaxonomy1.HexColor);
% % to add color for "others"
map = [map; [0.5 0.5 0.5]];
cmap = colormap(map);
abundanceTable = T3Interest1{:,2:end};

subplot(4,1,1)
h = bar(1:size(T3Interest1,1),  abundanceTable, 'stacked', 'EdgeColor', 'none');
for i = 1:length(h)
    set(h(i), 'FaceColor', cmap(i, :));
end
set(gca, 'Box', 'on', 'XColor', [0 0 0], 'YColor', [0 0 0]);
set(gca, 'XTickLabels', [])
xlim([0 size(T3Interest1,1)+1])
ylim([0 100])
ylabel('Abundance(%)', 'FontSize',12, 'FontWeight','bold', 'FontName', 'Arial');
title('Barplot of samples from experiments 3, 4, 7, 8 and 11', 'FontSize',12, 'FontWeight','bold', 'FontName', 'Arial')
hold on

subplot(4,1,2)
bar(1:size(T3Interest1,1),genusPercentage3.ido1, 'k')
ylabel('{\it ido1} expression', 'FontSize',12, 'FontWeight','bold', 'FontName', 'Arial')

subplot(4,1,3)
bar(1:size(T3Interest1,1),genusPercentage3.Simpson, 'k')
ylabel('Inverse Simpson Index',  'FontSize',12, 'FontWeight','bold', 'FontName', 'Arial')
set(gca, 'XTickLabels', [])
hold on

subplot(4,1,4)
bar(1:size(T3Interest1,1),genusPercentage3.x16S_qPCR_per_g, 'k')
set(gca, 'YScale', 'log')
set(gca, 'YTick', [10^8 10^9 10^10 10^11 10^12])
ylabel('16S copy number/g of feces',  'FontSize',12, 'FontWeight','bold', 'FontName', 'Arial')
set(gca, 'XTickLabels', [])

print -painters -dpdf 10_15_24.pdf


