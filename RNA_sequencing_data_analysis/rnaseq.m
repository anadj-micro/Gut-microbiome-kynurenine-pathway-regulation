%% This script processes RNA sequencing data from untreated, antibiotic-treated,
% and LPS-treated mice and analyzes the expression of the kynurenine pathway
% enzymes in the intestine

% load the table of sample names
tblSamples = readtable("tblSamples.csv");
tblSamples.Properties.RowNames = tblSamples.sample;
tblSamples.sample = categorical(tblSamples.sample);
tblSamples.treatment = categorical(tblSamples.treatment);
tblSamples.organ = categorical(tblSamples.organ);

%% load the normalized counts of antibiotics vs control (provided by Genewiz)
tblNormCounts1 = readtable('control-vs-antibiotics_deseq2_normalized_counts.csv');
tblNormCounts2 = readtable('control-vs-LPS_deseq2_normalized_counts.csv');

%% Load tables for gene ontology enrichemnt analysis
go_terms = readtable('go_terms.mgi.txt', 'ReadVariableNames', false);
go_terms.Properties.VariableNames = {'category' 'go' 'description'};
go_terms = go_terms(: ,{'category' 'go' 'description'});

mrk_ensembl = readtable('MRK_ENSEMBL.rpt.txt', 'ReadVariableNames', false);
mrk_ensembl.Properties.VariableNames{1} = 'mgi';
mrk_ensembl.Properties.VariableNames{2} = 'symbol';
mrk_ensembl.Properties.VariableNames{3} = 'description';
mrk_ensembl.Properties.VariableNames{6} = 'gene';
mrk_ensembl = mrk_ensembl(:, {'mgi' 'symbol' 'description' 'gene'});

gene_association = readtable('gene_association.mgi', 'ReadVariableNames', false, 'NumHeaderLines', 3, 'FileType', 'text');
gene_association.Properties.VariableNames{2} = 'mgi';
gene_association.Properties.VariableNames{4} = 'role';
gene_association.Properties.VariableNames{5} = 'go';
gene_association = gene_association(:, {'mgi' 'go' 'role'});

gene_association = innerjoin(gene_association, mrk_ensembl(:, {'mgi' 'gene'}));
gene_association = gene_association(:, {'gene' 'go' 'role'});

%% join the two deseq2 tables, remove the duplicates, and change names to more intuitive ones
tblNormCounts = innerjoin(tblNormCounts1, tblNormCounts2, 'Keys', 'Var1');
% find and remove duplicated columns
idxDuplicated = contains(tblNormCounts.Properties.VariableNames, 'tblNormCounts2');
tblNormCounts(:, idxDuplicated) = [];
tblNormCounts.Properties.VariableNames =...
    strrep(tblNormCounts.Properties.VariableNames, '_tblNormCounts1', '');
tblNormCounts.Properties.VariableNames =...
    strrep(tblNormCounts.Properties.VariableNames, 'II', 'L');
tblNormCounts.Properties.VariableNames =...
    strrep(tblNormCounts.Properties.VariableNames, 'I1_', 'A1_');
tblNormCounts.Properties.VariableNames =...
    strrep(tblNormCounts.Properties.VariableNames, 'I2_', 'A2_');
tblNormCounts.Properties.VariableNames =...
    strrep(tblNormCounts.Properties.VariableNames, 'I3_', 'A3_');

%% normalize each sample based on the geometric mean of counts
% estimate pseudo-reference with geometric mean row by row
pseudoRefSample = geomean(tblNormCounts{:, 2:end},2);
nz = pseudoRefSample > 0;
ratios = bsxfun(@rdivide,tblNormCounts{nz, 2:end},pseudoRefSample(nz));
sizeFactors = median(ratios,1);
tblNormCounts{:, 2:end} = bsxfun(@rdivide, tblNormCounts{:, 2:end}, sizeFactors);

%% stack the measurements made per mouse
tblNormCountsStack = stack(tblNormCounts, 2:width(tblNormCounts));
tblNormCountsStack.Properties.VariableNames = {'gene' 'sample' 'normCounts'};
tblNormCountsStack = innerjoin(tblNormCountsStack, tblSamples);
tblNormCountsStack.mouse = cellfun(@(x) x(1:2), cellstr(tblNormCountsStack.sample), 'UniformOutput', false);
tblNormCountsStack.sample = [];
tblNormCountsStack.treatment = [];
tblNormCountsUnstack = unstack(tblNormCountsStack, 'normCounts', 'mouse');

%% compute the means, dispersan and fold cange 
tblNormCountsUnstack.meanC = mean(tblNormCountsUnstack{:, {'C1' 'C2' 'C3'}}, 2);
tblNormCountsUnstack.meanA = mean(tblNormCountsUnstack{:, {'A1' 'A2' 'A3'}}, 2);
tblNormCountsUnstack.meanL = mean(tblNormCountsUnstack{:, {'L1' 'L2' 'L3'}}, 2);
tblNormCountsUnstack.dispC = std(tblNormCountsUnstack{:, {'C1' 'C2' 'C3'}}, 0, 2)./tblNormCountsUnstack.meanC;
tblNormCountsUnstack.dispA = std(tblNormCountsUnstack{:, {'A1' 'A2' 'A3'}}, 0, 2)./tblNormCountsUnstack.meanA;
tblNormCountsUnstack.dispL = std(tblNormCountsUnstack{:, {'L1' 'L2' 'L3'}}, 0, 2)./tblNormCountsUnstack.meanL;

% compute the mean and the log2FC
tblNormCountsUnstack.meanBaseA  = (tblNormCountsUnstack.meanA + tblNormCountsUnstack.meanC) / 2;
tblNormCountsUnstack.foldChangeA = tblNormCountsUnstack.meanA ./ tblNormCountsUnstack.meanC;
tblNormCountsUnstack.log2FCA = log2(tblNormCountsUnstack.foldChangeA);
tblNormCountsUnstack.meanBaseL  = (tblNormCountsUnstack.meanL + tblNormCountsUnstack.meanC) / 2;
tblNormCountsUnstack.foldChangeL = tblNormCountsUnstack.meanL ./ tblNormCountsUnstack.meanC;
tblNormCountsUnstack.log2FCL = log2(tblNormCountsUnstack.foldChangeL);

%% compute the p-value of the differential expression for each tissue
% compute the P-values (takes ~1 min)
organ = {'Intestine' 'Spleen' 'Liver'};
h = waitbar(0,'Computing p-values for diff expressed genes');
c = 0;
for i = 1:length(organ)
    c = c+1; 
    waitbar(c/(length(organ)*2),h);
    idx = tblNormCountsUnstack.organ == organ{i};
    tLocalA = nbintest(tblNormCountsUnstack{idx, {'A1' 'A2' 'A3'}},...
        tblNormCountsUnstack{idx, {'C1' 'C2' 'C3'}},'VarianceLink','LocalRegression');
    c = c+1; 
    waitbar(c/(length(organ)*2),h);
    tblNormCountsUnstack.pValueA(idx) = tLocalA.pValue;
    tLocalL = nbintest(tblNormCountsUnstack{idx, {'L1' 'L2' 'L3'}},...
        tblNormCountsUnstack{idx, {'C1' 'C2' 'C3'}},'VarianceLink','LocalRegression');
    tblNormCountsUnstack.pValueL(idx) = tLocalL.pValue;
end
close(h);

%% do a global multiple hypothesis correction
pValue = [tblNormCountsUnstack.pValueA; tblNormCountsUnstack.pValueL];
pValueAdujusted = mafdr(pValue,'BHFDR',true);

tblNormCountsUnstack.pValueAdjusA =...
    pValueAdujusted(1:height(tblNormCountsUnstack));
tblNormCountsUnstack.pValueAdjusL =...
    pValueAdujusted(height(tblNormCountsUnstack)+1:end);

%% analyze differential expression in intestine, liver and spleed
figure(101)
for i = 1:length(organ)
    idx = tblNormCountsUnstack.organ == organ{i};
    subplot(1, 3, i)
    plot(tblNormCountsUnstack.log2FCA(idx), tblNormCountsUnstack.log2FCL(idx), 'o',...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0.8 0.8 0.8]);
    hold on
    significantDe = tblNormCountsUnstack.pValueAdjusA < 0.05; % the genes statistically significant
    plot(tblNormCountsUnstack.log2FCA(significantDe&idx), tblNormCountsUnstack.log2FCL(significantDe&idx), 'o',...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'g');
    significantDe = tblNormCountsUnstack.pValueAdjusL <0.05;
    plot(tblNormCountsUnstack.log2FCA(significantDe&idx), tblNormCountsUnstack.log2FCL(significantDe&idx), 'o',...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b');
    significantDe = tblNormCountsUnstack.pValueAdjusA < 0.05 & tblNormCountsUnstack.pValueAdjusL <0.05;
    plot(tblNormCountsUnstack.log2FCA(significantDe&idx), tblNormCountsUnstack.log2FCL(significantDe&idx), 'o',...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k');
    hold off
    axis equal;
    grid on;
    xlabel('Antibiotics vs control')
    ylabel('LPS vs control')
    title(organ{i});
    set(gca, 'XLim', [-10 10], 'YLim', [-15 15])
end

%% choose genes to highlight
geneSet = mrk_ensembl(strcmp(mrk_ensembl.symbol, 'Ido1'), :);
% find go terms with this gene
%go = go_terms.go(strcmp(go_terms.description, 'tryptophan catabolic process to kynurenine'), :)
gos = gene_association.go(strcmp(gene_association.gene, geneSet.gene));
go_terms(ismember(go_terms.go, gos), :)

%
% Waterfall plots of diff expression
figure(102)
for i = 1:3
    % Abx
    subplot(3, 2, (i-1)*2 + 1)
    idx = find(tblNormCountsUnstack.organ == organ{i} & tblNormCountsUnstack.pValueAdjusA < 0.05);
    %idx = find(tblNormCountsUnstack.organ == organ{i});
    [~, iSort] = sort(tblNormCountsUnstack.log2FCA(idx)); 
    idx = idx(iSort);
    plot(tblNormCountsUnstack.log2FCA(idx), 'o',...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0.8 0.8 0.8]);
    hold on
    idxGene = find(ismember(tblNormCountsUnstack.gene(idx), geneSet.gene)); 
    plot(idxGene, tblNormCountsUnstack.log2FCA(idx(idxGene)), 'o',...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r');
    text(idxGene+30, tblNormCountsUnstack.log2FCA(idx(idxGene)), geneSet.symbol);    
    hold off
    title(['Abx in ' organ{i}])
    set(gca, 'YLim', [-10 15]); grid on; ylabel('log2(FC)')
    % LPS
    subplot(3, 2, i*2)
    idx = find(tblNormCountsUnstack.organ == organ{i} & tblNormCountsUnstack.pValueAdjusL < 0.05);
    [~, iSort] = sort(tblNormCountsUnstack.log2FCL(idx)); idx = idx(iSort);
    plot(tblNormCountsUnstack.log2FCL(idx), 'o',...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0.8 0.8 0.8]);
    hold on
    idxGene = find(ismember(tblNormCountsUnstack.gene(idx), geneSet.gene)); 
    plot(idxGene, tblNormCountsUnstack.log2FCL(idx(idxGene)), 'o',...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'r');    
    text(idxGene+30, tblNormCountsUnstack.log2FCA(idx(idxGene)), geneSet.symbol);
    hold off
    title(['LPS in ' organ{i}])
    set(gca, 'YLim', [-10 15]); grid on; ylabel('log2(FC)')
end


%% GO terms downreg and upreg in treatments and organs
treatmentForDe = {'Antibiotic' 'LPS'};
directionOfDe = {'down' 'up' };
% focus analysis only on biological processes
tblGoEnrichment = go_terms(strcmp(go_terms.category, 'Biological Process'), {'go' 'description'});
for i = 1:length(organ)
    idx = (tblNormCountsUnstack.organ == organ{i});
    for j = 1:length(treatmentForDe)
        for k = 1:length(directionOfDe)
            if j == 1 && k == 1
                idx2 = idx & tblNormCountsUnstack.pValueAdjusA < 0.05  &...
                    tblNormCountsUnstack.log2FCA < 0;
            elseif j == 1 && k == 2
                idx2 = idx & tblNormCountsUnstack.pValueAdjusA < 0.05  &...
                    tblNormCountsUnstack.log2FCA > 0;
            elseif j == 2 && k == 1
                idx2 = idx & tblNormCountsUnstack.pValueAdjusL < 0.05  &...
                    tblNormCountsUnstack.log2FCL < 0;
            elseif j == 2 && k == 2
                idx2 = idx & tblNormCountsUnstack.pValueAdjusL < 0.05  &...
                    tblNormCountsUnstack.log2FCL > 0;
            end            
            diffRegGenes = tblNormCountsUnstack{idx2, {'gene'}};      
            diffRegTable = gene_association;
            diffRegTable.de = ismember(diffRegTable.gene, diffRegGenes);
            goEnriched = grpstats(diffRegTable, 'go', 'sum', 'DataVars', 'de');

            M = length(unique(diffRegTable.gene)); % total number of genes
            K = sum(idx2); % total number of genes diff exp
            % compute pValue of GO using hyper geometric distribution
            goEnriched.pValue = hygepdf(goEnriched.sum_de, M, K, goEnriched.GroupCount);
            goEnriched = goEnriched(:, {'go' 'pValue'});
            if j == 1 && k == 1
                goEnriched.Properties.VariableNames{2} = ['pValueDownInAbx' organ{i}];
            elseif j == 1 && k == 2
                goEnriched.Properties.VariableNames{2} = ['pValueUpInAbx' organ{i}];
            elseif j == 2 && k == 1
                goEnriched.Properties.VariableNames{2} = ['pValueDownInLPS' organ{i}];
            elseif j == 2 && k == 2
                goEnriched.Properties.VariableNames{2} = ['pValueUpInLPS' organ{i}];
            end            
            tblGoEnrichment = outerjoin(tblGoEnrichment, goEnriched,...
                'Type', 'left', 'MergeKeys', true);
        end
    end
end
%% multiple hypothesis correction of GO association
pValue = tblGoEnrichment{:, 3:end};
pValueAdujusted = mafdr(pValue(:),'BHFDR',true);
tblGoEnrichmentAdjusted = tblGoEnrichment;
tblGoEnrichmentAdjusted{:, 3:end} = reshape(pValueAdujusted, size(pValue));

%% draw PCA of genes differentially expressed
mice = {'C1' 'C2' 'C3' 'A1' 'A2' 'A3' 'L1' 'L2' 'L3'};
treatment = {'Control' 'Control' 'Control'...
    'Antibiotics' 'Antibiotics' 'Antibiotics'...
    'LPS' 'LPS' 'LPS'}';

significantDe = tblNormCountsUnstack.pValueAdjusA < 0.05 |...
    tblNormCountsUnstack.pValueAdjusL < 0.05;


figure(2)

idx = tblNormCountsUnstack.organ == 'Intestine' ;
m = zscore(tblNormCountsUnstack{idx&significantDe,...
    mice}');
[coeff, score, latent, tsquared, explained] = pca(m, 'Centered',true);
subplot(1, 3, 1)
gscatter(score(:, 1), score(:, 2), treatment)
text(score(:, 1), score(:, 2), mice)
xlabel(sprintf('PC1 (%0.2f%% EV)', explained(1)))
ylabel(sprintf('PC2 (%0.2f%% EV)', explained(2)))
title('Intestine')
axis square; grid on;

idx = tblNormCountsUnstack.organ == 'Spleen';
m = zscore(tblNormCountsUnstack{idx&significantDe,...
    mice}');
[coeff, score, latent, tsquared, explained] = pca(m, 'Centered',true);
subplot(1, 3, 2)
gscatter(score(:, 1), score(:, 2), treatment)
text(score(:, 1), score(:, 2), mice)
xlabel(sprintf('PC1 (%0.2f%% EV)', explained(1)))
ylabel(sprintf('PC2 (%0.2f%% EV)', explained(2)))
title('Spleen')
axis square; grid on;

idx = tblNormCountsUnstack.organ == 'Liver';
m = zscore(tblNormCountsUnstack{idx&significantDe,...
    mice}');
[coeff, score, latent, tsquared, explained] = pca(m, 'Centered',true);
subplot(1, 3, 3)
gscatter(score(:, 1), score(:, 2), treatment)
text(score(:, 1), score(:, 2), mice)
xlabel(sprintf('PC1 (%0.2f%% EV)', explained(1)))
ylabel(sprintf('PC2 (%0.2f%% EV)', explained(2)))
title('Liver')
axis square; grid on;

%% IDO1 level
ido1Idx = mrk_ensembl.gene(strcmp(mrk_ensembl.symbol, 'Ido1'));
ido1Idx = strcmp(tblNormCountsUnstack.gene, ido1Idx);

figure(3)
idx = tblNormCountsUnstack.organ == 'Intestine' ;
m = zscore(tblNormCountsUnstack{idx&significantDe, mice}');
[coeff, score, latent, tsquared, explained] = pca(m, 'Centered',true);
scatter(score(:, 1), score(:, 2), 400,...
    tblNormCountsUnstack{idx&ido1Idx, mice}', 'filled')
text(score(:, 1)+1, score(:, 2)+1,mice)
xlabel(sprintf('PC1 (%0.2f%% EV)', explained(1)))
ylabel(sprintf('PC2 (%0.2f%% EV)', explained(2)))
title('Intestine')
axis square; grid on;

%% figures for paper
% AVN impacts genes expressed in host organs.

% get all the genes in the kynurine pathway
kpGo = tblGoEnrichmentAdjusted.go(strcmp(tblGoEnrichmentAdjusted.description, 'tryptophan catabolic process to kynurenine'), :);
kpGenes = mrk_ensembl.gene(ismember(mrk_ensembl.gene, gene_association.gene(strcmp(gene_association.go, kpGo))));


%% AVN
maxLog2Fc = 5;
log2FC = tblNormCountsUnstack.log2FCA;
log2FC(log2FC<-maxLog2Fc) = -maxLog2Fc;
log2FC(log2FC>maxLog2Fc) = maxLog2Fc;
meanCounts = (tblNormCountsUnstack.meanA + tblNormCountsUnstack.meanC)/2;

figure(102)
set(gcf, 'Position', [1000         261         453         977], 'Renderer', 'painters')
for i = 1:length(organ)
    idx = tblNormCountsUnstack.organ == organ{i};
    subplot(3, 1, i)
    scatter(log2FC(idx), meanCounts(idx), 'o',...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerFaceAlpha',0.2);
     hold on
     significantDe = tblNormCountsUnstack.pValueAdjusA < 0.05; % the genes statistically significant
     scatter(log2FC(significantDe&idx),...
         meanCounts(significantDe&idx), 'o',...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerFaceAlpha',0.2);
    kpGeneInOrgan = find(ismember(tblNormCountsUnstack.gene, kpGenes) & idx);
    scatter(log2FC(kpGeneInOrgan),...
        meanCounts(kpGeneInOrgan), 'o',...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerFaceAlpha',0.8);
    % plot gene name next to red dot
    for j = 1:length(kpGeneInOrgan)
        symbol = mrk_ensembl.symbol(strcmp(mrk_ensembl.gene,...
            tblNormCountsUnstack.gene(kpGeneInOrgan(j))));
        text(log2FC(kpGeneInOrgan(j))+0.2,...
        meanCounts(kpGeneInOrgan(j))+j, symbol, 'Color','r')
    end
    hold off
    grid on;

    title(sprintf('%s (%d genes down, %d genes up)',...
        organ{i}, sum(log2FC(significantDe&idx)<0),...
        sum(log2FC(significantDe&idx)>0)));
    xlabel('log_2(FC) AVN vs Untreated')
    %ylabel('Mean normalized abunance counts')
    set(gca, 'XLim', [-maxLog2Fc maxLog2Fc])
    set(gca, 'YScale', 'log', 'YLim', [0.1 1E6])
    %fontsize('scale', 1.5)
    box on
end
print(gcf, 'genesInAVN_03_11_26_alone.png', '-dpng');

%%
% LPS
figure(103)
set(gcf, 'Position', [1000         261         453         977], 'Renderer', 'painters')

log2FC = tblNormCountsUnstack.log2FCL;
log2FC(log2FC<-maxLog2Fc) = -maxLog2Fc;
log2FC(log2FC>maxLog2Fc) = maxLog2Fc;

meanCounts = (tblNormCountsUnstack.meanL + tblNormCountsUnstack.meanC)/2;

for i = 1:length(organ)
    idx = tblNormCountsUnstack.organ == organ{i};
    subplot(3, 1, i)
    scatter(log2FC(idx), meanCounts(idx), 'o',...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerFaceAlpha',0.2);
     hold on
     significantDe = tblNormCountsUnstack.pValueAdjusL < 0.05; % the genes statistically significant
     scatter(log2FC(significantDe&idx),...
         meanCounts(significantDe&idx), 'o',...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerFaceAlpha',0.2);
    kpGeneInOrgan = find(ismember(tblNormCountsUnstack.gene, kpGenes) & idx);
    scatter(log2FC(kpGeneInOrgan),...
        meanCounts(kpGeneInOrgan), 'o',...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerFaceAlpha',0.8);
    % plot gene name next to red dot
    for j = 1:length(kpGeneInOrgan)
        symbol = mrk_ensembl.symbol(strcmp(mrk_ensembl.gene,...
            tblNormCountsUnstack.gene(kpGeneInOrgan(j))));
        text(log2FC(kpGeneInOrgan(j))+0.2,...
        meanCounts(kpGeneInOrgan(j)), symbol, 'Color','r')
    end
    hold off
    grid on;

    title(sprintf('%s (%d genes down, %d genes up)',...
        organ{i}, sum(log2FC(significantDe&idx)<0),...
        sum(log2FC(significantDe&idx)>0)));
    xlabel('log_2(FC) LPS vs Untreated')
    %ylabel('Mean normalized abunance counts')
    set(gca, 'XLim', [-maxLog2Fc maxLog2Fc])
    set(gca, 'YScale', 'log', 'YLim', [0.1 1E6])
    %fontsize('scale', 1.5)
    box on
end
print(gcf, 'genesInLPS_03_11_26_alone.png', '-dpng');

%% heat map of GO enrichment
%t = tblGoEnrichmentAdjusted(cellfun(@(x) ~isempty(x), regexpi(tblGoEnrichmentAdjusted.description, 'tryptophan')), :);
t = tblGoEnrichmentAdjusted(contains(tblGoEnrichmentAdjusted.description, {'tryptophan' 'L-tryptophan'}, 'IgnoreCase', true), :);
t = t(~contains(t.description, 'tryptophan' + lettersPattern, 'IgnoreCase', true), :);
t = t(~contains(t.description, '-L-tryptophan' , 'IgnoreCase', true), :);

columnLabels = t.Properties.VariableNames(3:end);
columnLabels = strrep(columnLabels, 'pValue', '');
columnLabels = strrep(columnLabels, 'Abx', 'AVN');

figure(104)
set(gcf, 'Position', [1000         906        1019         332])
h = heatmap(columnLabels, join([t.go, t.description],' ',2), t{:, 3:end});
h.CellLabelFormat = '%.3f';
h.Colormap = [0 0 0; 1 1 1; 1 1 1];
colorbar('off')

