%% Script to plot a figure with sublots where each subplot is a barplot with:
% a) ASV abundances in samples from mice treated with different antibiotic
% regimes; b) IDO1 expression levels in corresponding colon samples; 
% c) Inverse Simpson in fecal samples; and d) 16S rDNA copy number per gram
% of fecal matter for each of the animals. 

close all
clear all

% read taxonomy table
tbltaxonomy = readtable('tblASVtaxonomy_dada2.csv');

% read count table and unstack
tblcounts = readtable('tblASVcounts.csv', 'Delimiter', ',');
tblcounts = unstack(tblcounts, 'count', 'ASV');

% only keep taxonomic ASVs that appear in the count table
tbltaxonomy = tbltaxonomy(ismember(tbltaxonomy.ASV,tblcounts.Properties.VariableNames(2:end)), :);

% convert counts to relative abundance table
counts_matrix = tblcounts{:, 2:end}; % the first column is "SampleID"
counts_matrix(isnan(counts_matrix)) = 0; % missing count value is filled with 0
relative_abundance_matrix = counts_matrix ./ sum(counts_matrix, 2); % convert to relative abundance

% assign colors to asvs
%
% make sure that ASVs in tbltaxonomy and relative_abundance_matrix are the
% same (sometimes one table has more ASVs than the other)
%
% tbltaxonomy_new appends two columns to tbltaxonomy: HexColor and
% ColorOrder. HexColor is the color for each ASV and ColorOrder specifies
% the position in stacked barplot (the lower color order, the lower
% position is more towards to the bottom part of the barplot)
%
% color_spectrum_dict links the gradients of colors (maximum 21) assigned to each
% taxonomy
[tbltaxonomy_new, color_spectrum_dict] = brad_assign_colors_Ana_colors(tbltaxonomy, relative_abundance_matrix);

% save new taxonomy file
writetable(tbltaxonomy_new, 'tblASVtaxonomy_new.csv', 'Delimiter', ',', 'WriteRowNames',true);
