%% This script analyzes kynurenine-to-tryptophan ratio in plasma and 
% stool of allogeneic hematopoietic stem cell transplant (allo-HCT) patients, 
% and correlates it with neopterin (a marker for immune activation). 
% Allo-HCT patient data has been published previously (PMID 37329880 and 
% PMID 33654104).

%% load the metabolites analyzed for samples
tblMetabolites = readtable('tblMetabolitesUnstackedStoolPlasma.txt');

%% get names of metabolites
uMetabolites = tblMetabolites.Properties.VariableNames(3:end);

%% load the taxUmap coordinates (This is from a PMID 37329880)
tblTaxUmap = readtable('current_taxumap_embedding.csv', 'ReadVariableNames',true);
tblTaxUmap.Properties.VariableNames = {'SampleID' 'taxUmapX' 'taxUmapy'};

% get names of samples missing from taxUmap
disp('Tese 4 samples are missing from Scientific Data/TaxUMAP:')
samplesMissing = setdiff(tblMetabolites.SampleID, tblTaxUmap.SampleID)
% add taxumap coordinates to samples
tblMetabolites = innerjoin(tblMetabolites, tblTaxUmap);

%% get microbiota compositional information 
tblCountsAsvMelt = readtable("scientific_data_tables/tblcounts_asv_melt.csv");
tblDepth = groupsummary(tblCountsAsvMelt, 'SampleID', 'sum', 'Count');
tblRelAbundAsvMelt = innerjoin(tblCountsAsvMelt, tblDepth);
tblRelAbundAsvMelt.relAbund = tblRelAbundAsvMelt.Count ./ tblRelAbundAsvMelt.sum_Count;
tblRelAbundAsvMelt.relAbund2 = tblRelAbundAsvMelt.relAbund .^2;

%% biodiversity  (Inverse Simpson)
tblSamples = groupsummary(tblRelAbundAsvMelt, "SampleID", "sum","relAbund2");
tblSamples.inverseSimpson = 1./tblSamples.sum_relAbund2;
tblSamples.GroupCount = [];
tblSamples.sum_relAbund2 = [];
tblSamples = innerjoin(tblSamples, tblTaxUmap);
tblSamples.Properties.RowNames = tblSamples.SampleID;

%% add diversity to metabolite tables
tblMetabolites.inverseSimpson = tblSamples.inverseSimpson(tblMetabolites.SampleID);

%% find top ASV in sample (for coloring of ther samples based on the dominant taxa)
topAsv = sortrows(tblCountsAsvMelt, {'SampleID' 'Count'}, {'ascend' 'descend'});
[~, idxTop] = unique(topAsv.SampleID);
topAsv = topAsv(idxTop, :);
% get their color
tblASVtaxonomy = readtable("scientific_data_tables/tblASVtaxonomy_silva132_v4v5_filter.csv");
tblASVtaxonomy.Properties.RowNames = tblASVtaxonomy.ASV;
% extract the color
topAsv.color = tblASVtaxonomy.HexColor(topAsv.ASV);
% covert hex to rgb
topAsv.color = hex2rgb(topAsv.color);
topAsv.Properties.RowNames = topAsv.SampleID;
% add to metabolite table
tblMetabolites.color = topAsv.color(tblMetabolites.SampleID, :);

%% add kyn/trp ratios
tblMetabolites.kyn2trpRatio_stool = tblMetabolites.kynurenine_stool./tblMetabolites.tryptophan_stool;
tblMetabolites.kyn2trpRatio_plasma = tblMetabolites.kynurenine_plasma./tblMetabolites.tryptophan_plasma;

%% Make a figure

% A. draw the taxUMAP showing 92 samples selected for paired fecal and blood
% LC-MS/MS

% this dims the colors of the samples shown in the taxUMAP
% to highlight the 92 samples used in metabolite analysis
c = topAsv.color(tblSamples.SampleID, :);
c = c + (1-c) * 0.5;

figure(4)
set(gcf, 'Position', [1000        1015         784         322])
subplot(1, 3, 1)
scatter(tblSamples.taxUmapX, tblSamples.taxUmapy, 2,...
    c, "filled")
colormap gray
hold on
scatter(tblMetabolites.taxUmapX, tblMetabolites.taxUmapy, 10,...
    tblMetabolites.color, 'filled', 'MarkerEdgeColor','k')
hold off
axis equal tight
%grid on
xlabel ('TaxUMAP 1')
ylabel ('TaxUMAP 2')
set(gca, 'XTickLabels', [], 'YTickLabels', [])
axis off

% B: Inflamation correlates positively
% with KYN/TRP systemically...
subplot(1, 3, 2)
xMetabolite = 'neopterin_plasma';
yMetabolite = 'kyn2trpRatio_plasma';
scatter(tblMetabolites{:, xMetabolite}+1,...
    tblMetabolites{:, yMetabolite}+eps,...
    50, tblMetabolites.color, 'filled', 'MarkerEdgeColor','k')
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel(xMetabolite, 'Interpreter','none')
ylabel(yMetabolite, 'Interpreter','none')
grid
[r, p] = corr(tblMetabolites{:, xMetabolite},...
    tblMetabolites{:, yMetabolite},...
    "Type","Spearman", "Rows","pairwise")
title(sprintf('R=%0.2f, p=%1.0e',r, p))

% C: ...but not in gut, suggesting
% microbial influence.
subplot(1, 3, 3)
xMetabolite = 'neopterin_stool';
yMetabolite = 'kyn2trpRatio_stool';
scatter(tblMetabolites{:, xMetabolite}+10,...
    tblMetabolites{:, yMetabolite}+1e-4,...
    50, tblMetabolites.color, 'filled', 'MarkerEdgeColor','k')
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel(xMetabolite, 'Interpreter','none')
ylabel(yMetabolite, 'Interpreter','none')
grid
[r, p] = corr(tblMetabolites{:, xMetabolite},...
    tblMetabolites{:, yMetabolite},...
    "Type","Spearman", "Rows","pairwise")
title(sprintf('R=%0.2f, p=%1.0e',r, p))

% save the figure to .eps file for editing in adobe illustrator
print('trpCat_alloHCT.eps','-depsc')
