function run_figureS2_rnaseq()
%RUN_FIGURES2_RNASEQ Reproduce Figure S2C-D using the frozen active GO set.
%
% This dated analysis preserves the July 14 run and repeats differential
% expression and GO Biological Process over-representation testing from the
% archived integer counts.  It replaces the obsolete GO:0019441 display with
% the active GO:0006569 term and a frozen current MGI/GOC annotation snapshot.

p = localPaths();
diaryFile = fullfile(p.logs, 'matlab_run.log');
if isfile(diaryFile)
    delete(diaryFile);
end
diary(diaryFile);
cleanup = onCleanup(@() diary('off'));

fprintf('Figure S2 RNA-seq rerun started: %s\n', string(datetime('now')));
fprintf('MATLAB: %s\n', version);

[allStats, sampleManifest] = recomputeDifferentialExpression(p);
writetable(allStats, fullfile(p.results, 'fig1_rnaseq_gene_statistics_current_go.csv'));
writetable(sampleManifest, fullfile(p.results, 'fig1_rnaseq_sample_manifest_current_go.csv'));

[terms, associations, release] = readCurrentGoInputs(p);
writetable(release, fullfile(p.results, 'go_annotation_release.csv'));

targetGO = "GO:0006569";
targetTerm = terms(terms.GO == targetGO, :);
assert(height(targetTerm) == 1, 'Current GO term GO:0006569 was not found exactly once.');

targetGenes = unique(associations.Gene(associations.GO == targetGO));
assert(~isempty(targetGenes), 'Current MGI annotations contain no genes for GO:0006569.');
targetGeneTable = unique(allStats(ismember(allStats.Gene, targetGenes), {'Gene','Symbol'}), 'rows');
targetGeneTable = sortrows(targetGeneTable, 'Symbol');
writetable(targetGeneTable, fullfile(p.results, 'fig1_rnaseq_ltryptophan_catabolic_process_genes.csv'));

targetStatistics = allStats(ismember(allStats.Gene, targetGenes), :);
writetable(targetStatistics, fullfile(p.results, 'fig1_rnaseq_ltryptophan_catabolic_process_gene_statistics.csv'));

[goResults, testSummary] = computeGoEnrichment(allStats, terms, associations);
writetable(goResults, fullfile(p.results, 'fig1_rnaseq_go_enrichment_current_mgi.csv'));
writetable(testSummary, fullfile(p.results, 'fig1_rnaseq_go_enrichment_test_summary.csv'));
targetResults = goResults(goResults.GO == targetGO, :);
writetable(targetResults, fullfile(p.results, 'fig1_rnaseq_ltryptophan_catabolic_process_enrichment.csv'));

for contrast = ["AVN", "LPS"]
    intestineStats = allStats(allStats.Organ == "Intestine" & ...
        allStats.Contrast == contrast, :);
    makeRnaPanel(intestineStats, targetGenes, contrast, p);
end

avnTarget = targetResults(targetResults.Organ == "Intestine" & ...
    targetResults.Contrast == "AVN" & targetResults.Direction == "Down", :);
lpsTarget = targetResults(targetResults.Organ == "Intestine" & ...
    targetResults.Contrast == "LPS" & targetResults.Direction == "Down", :);
fprintf('Current GO term: %s (%s); mapped gene set = %d genes.\n', ...
    targetGO, targetTerm.Description, height(targetGeneTable));
disp(targetGeneTable);
fprintf('Intestinal AVN-down: overlap %d/%d, P = %.6g, FDR = %.6g.\n', ...
    avnTarget.Overlap, avnTarget.SetSize, avnTarget.PValue, avnTarget.FDR);
fprintf('Intestinal LPS-down: overlap %d/%d, P = %.6g, FDR = %.6g.\n', ...
    lpsTarget.Overlap, lpsTarget.SetSize, lpsTarget.PValue, lpsTarget.FDR);
fprintf('Figure S2 RNA-seq rerun completed: %s\n', string(datetime('now')));
end

function p = localPaths()
codeDir = fileparts(mfilename('fullpath'));
p.run = codeDir;
p.inputs = fullfile(codeDir, 'current_go_inputs');
p.results = fullfile(p.run, 'results');
p.panels = fullfile(p.run, 'panels');
p.logs = fullfile(p.run, 'logs');
for folder = {p.results,p.panels,p.logs}
    if ~isfolder(folder{1})
        mkdir(folder{1});
    end
end
end

function [allStats, manifest] = recomputeDifferentialExpression(p)
avn = readtable(fullfile(p.inputs, 'raw_counts_C_AVN.csv'), ...
    'VariableNamingRule', 'preserve');
lps = readtable(fullfile(p.inputs, 'raw_counts_C_LPS.csv'), ...
    'VariableNamingRule', 'preserve');
avn.Properties.VariableNames{1} = 'Gene';
lps.Properties.VariableNames{1} = 'Gene';
avn.Gene = string(avn.Gene);
lps.Gene = string(lps.Gene);
assert(numel(unique(avn.Gene)) == height(avn) && numel(unique(lps.Gene)) == height(lps), ...
    'RNA-seq gene identifiers must be unique.');

commonGenes = intersect(avn.Gene, lps.Gene, 'stable');
avn = avn(ismember(avn.Gene, commonGenes), :);
lps = lps(ismember(lps.Gene, commonGenes), :);
[present, lpsOrder] = ismember(avn.Gene, lps.Gene);
assert(all(present), 'Not all AVN genes occur in the LPS count matrix.');
lps = lps(lpsOrder, :);
assert(isequal(avn.Gene, lps.Gene), 'Raw-count gene rows do not align.');

controlNames = avn.Properties.VariableNames(startsWith(avn.Properties.VariableNames, 'C'));
assert(isequal(avn{:,controlNames}, lps{:,controlNames}), ...
    'Control counts differ between archived contrast files.');
avnNames = avn.Properties.VariableNames(startsWith(avn.Properties.VariableNames, 'I'));
lpsNames = lps.Properties.VariableNames(startsWith(lps.Properties.VariableNames, 'II'));
counts = [avn{:,controlNames}, avn{:,avnNames}, lps{:,lpsNames}];
sampleNames = [controlNames, avnNames, lpsNames];
assert(all(isfinite(counts), 'all') && all(counts >= 0, 'all') && ...
    all(mod(counts, 1) == 0, 'all'), 'Raw RNA-seq matrix must contain nonnegative integers.');

manifest = makeSampleInfo(sampleNames);
sizeFactors = deseqSizeFactors(counts);
manifest.DisplaySizeFactor = sizeFactors';
normalized = counts ./ sizeFactors;

mapping = readtable(fullfile(p.inputs, 'MRK_ENSEMBL.rpt'), ...
    'FileType', 'text', 'ReadVariableNames', false, 'VariableNamingRule', 'preserve');
mapping.Properties.VariableNames{1} = 'MGI';
mapping.Properties.VariableNames{2} = 'Symbol';
mapping.Properties.VariableNames{6} = 'Gene';
mapping.MGI = string(mapping.MGI);
mapping.Gene = string(mapping.Gene);
mapping.Symbol = string(mapping.Symbol);
mapping = unique(mapping(:, {'Gene','Symbol'}), 'rows', 'stable');
[~, firstSymbolRow] = unique(mapping.Gene, 'stable');
mapping = mapping(firstSymbolRow, :);

allStats = table();
for organ = ["Intestine", "Liver", "Spleen"]
    controlIdx = manifest.Organ == organ & manifest.Treatment == "Control";
    for contrast = ["AVN", "LPS"]
        treatmentIdx = manifest.Organ == organ & manifest.Treatment == contrast;
        assert(sum(controlIdx) == 3 && sum(treatmentIdx) == 3, ...
            'Expected three biological replicates per group and organ.');
        controlCounts = counts(:, controlIdx);
        treatmentCounts = counts(:, treatmentIdx);
        eligible = sum(controlCounts, 2) + sum(treatmentCounts, 2) > 0;
        pValue = nan(height(avn), 1);
        test = nbintest(treatmentCounts(eligible, :), controlCounts(eligible, :), ...
            'VarianceLink', 'LocalRegression');
        pValue(eligible) = test.pValue;
        qValue = bhAdjust(pValue);
        meanControl = mean(normalized(:, controlIdx), 2);
        meanTreatment = mean(normalized(:, treatmentIdx), 2);
        log2FC = log2((meanTreatment + 0.5) ./ (meanControl + 0.5));
        block = table(avn.Gene, repmat(organ, height(avn), 1), ...
            repmat(contrast, height(avn), 1), meanControl, meanTreatment, ...
            log2FC, pValue, qValue, 'VariableNames', ...
            {'Gene','Organ','Contrast','MeanControl','MeanTreatment','Log2FC','PValue','FDR'});
        block = outerjoin(block, mapping, 'Keys', 'Gene', 'MergeKeys', true, ...
            'Type', 'left');
        allStats = [allStats; block]; %#ok<AGROW>
    end
end
end

function info = makeSampleInfo(names)
n = numel(names);
sample = strings(n,1);
treatment = strings(n,1);
organ = strings(n,1);
for i = 1:n
    token = string(names{i});
    parts = split(token, '-');
    mouse = parts(1);
    tissue = parts(2);
    if startsWith(mouse, "II")
        treatment(i) = "LPS";
        mouse = extractAfter(mouse, 1);
    elseif startsWith(mouse, "I")
        treatment(i) = "AVN";
        mouse = "A" + extractAfter(mouse, 1);
    else
        treatment(i) = "Control";
    end
    switch tissue
        case "I1"
            organ(i) = "Intestine";
        case "L1"
            organ(i) = "Liver";
        case "S1"
            organ(i) = "Spleen";
        otherwise
            error('Unexpected tissue code: %s', tissue);
    end
    sample(i) = mouse + "_" + tissue;
end
info = table(sample, string(names(:)), treatment, organ, ...
    'VariableNames', {'Sample','SourceColumn','Treatment','Organ'});
end

function [terms, associations, release] = readCurrentGoInputs(p)
gafGz = fullfile(p.inputs, 'mgi.gaf.gz');
oboPath = fullfile(p.inputs, 'go-basic.obo');
mapPath = fullfile(p.inputs, 'MRK_ENSEMBL.rpt');
assert(isfile(gafGz) && isfile(oboPath) && isfile(mapPath), ...
    'One or more frozen GO/MGI inputs are missing.');

cacheDir = fullfile(p.run, 'cache');
if ~isfolder(cacheDir)
    mkdir(cacheDir);
end
gafPath = fullfile(cacheDir, 'mgi.gaf');
if ~isfile(gafPath)
    gunzip(gafGz, cacheDir);
end

gaf = readtable(gafPath, 'FileType', 'text', 'Delimiter', '\t', ...
    'ReadVariableNames', false, 'CommentStyle', '!', 'VariableNamingRule', 'preserve');
assert(width(gaf) >= 9, 'Current MGI GAF must contain at least nine columns.');
gaf.Properties.VariableNames(1:9) = {'Database','MGI','Symbol','Qualifier','GO', ...
    'Reference','Evidence','WithFrom','Aspect'};
gaf.MGI = string(gaf.MGI);
gaf.GO = string(gaf.GO);
gaf.Qualifier = string(gaf.Qualifier);
gaf.Aspect = string(gaf.Aspect);
gaf = gaf(gaf.Aspect == "P" & ~contains(gaf.Qualifier, "NOT"), ...
    {'MGI','GO'});

mapping = readtable(mapPath, 'ReadVariableNames', false, ...
    'FileType', 'text', 'VariableNamingRule', 'preserve');
mapping.Properties.VariableNames{1} = 'MGI';
mapping.Properties.VariableNames{6} = 'Gene';
mapping.MGI = string(mapping.MGI);
mapping.Gene = string(mapping.Gene);
mapping = unique(mapping(:, {'MGI','Gene'}), 'rows');
associations = innerjoin(gaf, mapping, 'Keys', 'MGI');
associations = unique(associations(:, {'Gene','GO'}), 'rows');

terms = readBiologicalProcessTerms(oboPath);
associations = associations(ismember(associations.GO, terms.GO), :);

gafHeader = readGafHeader(gafPath);
gafDate = firstHeaderValue(gafHeader, "!date-generated:");
gafGoVersion = firstHeaderValue(gafHeader, "!go-version:");
oboHeader = readlines(oboPath, 'EmptyLineRule', 'skip');
oboRelease = firstHeaderValue(oboHeader, "data-version:");
release = table("Gene Ontology Consortium", ...
    "https://current.geneontology.org/annotations/mgi.gaf.gz", gafDate, gafGoVersion, ...
    "GO basic ontology", "https://current.geneontology.org/ontology/go-basic.obo", ...
    oboRelease, "2026-07-16", ...
    'VariableNames', {'AnnotationProvider','AnnotationURL','GAFGenerated', ...
    'GAFGoVersion','OntologyProduct','OntologyURL','OntologyRelease','AccessDate'});
end

function lines = readGafHeader(gafPath)
fid = fopen(gafPath, 'r');
assert(fid >= 0, 'Could not open current MGI GAF: %s', gafPath);
cleanup = onCleanup(@() fclose(fid));
lines = strings(0, 1);
while true
    raw = fgetl(fid);
    if ~ischar(raw) || ~startsWith(string(raw), "!")
        break
    end
    lines(end+1,1) = string(raw); %#ok<AGROW>
end
end

function terms = readBiologicalProcessTerms(oboPath)
lines = readlines(oboPath);
termStart = find(lines == "[Term]");
termEnd = [termStart(2:end)-1; numel(lines)];
go = strings(numel(termStart), 1);
description = strings(numel(termStart), 1);
namespace = strings(numel(termStart), 1);
obsolete = false(numel(termStart), 1);
for i = 1:numel(termStart)
    block = lines(termStart(i):termEnd(i));
    go(i) = firstHeaderValue(block, "id:");
    description(i) = firstHeaderValue(block, "name:");
    namespace(i) = firstHeaderValue(block, "namespace:");
    obsolete(i) = any(block == "is_obsolete: true");
end
terms = table(go, description, 'VariableNames', {'GO','Description'});
terms = terms(namespace == "biological_process" & ~obsolete & go ~= "", :);
assert(any(terms.GO == "GO:0006569"), 'GO:0006569 is not active in the frozen ontology.');
assert(~any(terms.GO == "GO:0019441"), 'Obsolete GO:0019441 was incorrectly retained.');
end

function value = firstHeaderValue(lines, prefix)
hit = lines(startsWith(lines, prefix));
if isempty(hit)
    value = "";
else
    value = strtrim(extractAfter(hit(1), strlength(prefix)));
end
end

function [result, testSummary] = computeGoEnrichment(stats, terms, associations)
% The direct, active annotations in the frozen current MGI GAF are used.
% This deliberately matches the July 14 direct-annotation test while updating
% the term set, annotation snapshot, gene mapping, and term nomenclature.
[isMapped, termIndex] = ismember(associations.GO, terms.GO);
associations = associations(isMapped, :);
termIndex = termIndex(isMapped);
associations.TermIndex = termIndex;
allAnnotatedGenes = unique(associations.Gene);
[~, associations.GeneIndex] = ismember(associations.Gene, allAnnotatedGenes);
associations = sortrows(associations, 'TermIndex');

activeIndex = unique(associations.TermIndex);
termGeneLists = cell(height(terms), 1);
groupStart = [1; find(diff(associations.TermIndex) ~= 0) + 1];
groupEnd = [groupStart(2:end)-1; height(associations)];
for i = 1:numel(groupStart)
    termGeneLists{associations.TermIndex(groupStart(i))} = unique( ...
        associations.GeneIndex(groupStart(i):groupEnd(i)));
end

result = table();
testSummary = table();
for organ = unique(stats.Organ, 'stable')'
    for contrast = unique(stats.Contrast, 'stable')'
        block = stats(stats.Organ == organ & stats.Contrast == contrast, :);
        [inUniverse, universeIndex] = ismember( ...
            unique(block.Gene(isfinite(block.PValue))), allAnnotatedGenes);
        universeIndex = unique(universeIndex(inUniverse));
        universeMask = false(numel(allAnnotatedGenes), 1);
        universeMask(universeIndex) = true;
        M = numel(universeIndex);
        for direction = ["Down", "Up"]
            if direction == "Down"
                selected = block.Gene(block.FDR < 0.05 & block.Log2FC < 0);
            else
                selected = block.Gene(block.FDR < 0.05 & block.Log2FC > 0);
            end
            [isSelected, selectedIndex] = ismember(unique(selected), allAnnotatedGenes);
            selectedIndex = unique(selectedIndex(isSelected));
            selectedIndex = selectedIndex(universeMask(selectedIndex));
            selectedMask = false(numel(allAnnotatedGenes), 1);
            selectedMask(selectedIndex) = true;
            K = numel(selectedIndex);

            nTerms = numel(activeIndex);
            setSize = zeros(nTerms, 1);
            overlap = zeros(nTerms, 1);
            pValue = ones(nTerms, 1);
            for j = 1:nTerms
                genes = termGeneLists{activeIndex(j)};
                setSize(j) = sum(universeMask(genes));
                overlap(j) = sum(selectedMask(genes));
                if K > 0 && setSize(j) > 0 && overlap(j) > 0
                    pValue(j) = hygecdf(overlap(j)-1, M, K, setSize(j), 'upper');
                end
            end
            rows = terms(activeIndex, :);
            rows.Organ = repmat(organ, nTerms, 1);
            rows.Contrast = repmat(contrast, nTerms, 1);
            rows.Direction = repmat(direction, nTerms, 1);
            rows.UniverseSize = repmat(M, nTerms, 1);
            rows.SelectedGenes = repmat(K, nTerms, 1);
            rows.SetSize = setSize;
            rows.Overlap = overlap;
            rows.PValue = pValue;
            rows.FDR = bhAdjust(pValue);
            result = [result; rows]; %#ok<AGROW>
            testSummary = [testSummary; table(organ, contrast, direction, M, K, nTerms, ...
                'VariableNames', {'Organ','Contrast','Direction','UniverseSize', ...
                'SelectedGenes','TestedTerms'})]; %#ok<AGROW>
        end
    end
end
end

function makeRnaPanel(stats, targetGenes, contrast, p)
stats.PlotMean = (stats.MeanControl + stats.MeanTreatment) / 2;
fig = figure('Visible', 'off', 'Units', 'inches', 'Position', [1 1 3.6 3.3]);
significant = stats.FDR < 0.05;
target = ismember(stats.Gene, targetGenes);
background = find(~significant & ~target);
if numel(background) > 5000
    rng(20260716, 'twister');
    background = background(randperm(numel(background), 5000));
end
scatter(stats.Log2FC(background), log10(stats.PlotMean(background) + 1), ...
    10, [0.70 0.70 0.70], 'filled', 'MarkerFaceAlpha', 0.28);
hold on;
scatter(stats.Log2FC(significant), log10(stats.PlotMean(significant) + 1), ...
    13, [0.25 0.25 0.25], 'filled', 'MarkerFaceAlpha', 0.55);
scatter(stats.Log2FC(target), log10(stats.PlotMean(target) + 1), ...
    30, [0.75 0.15 0.15], 'filled', 'MarkerEdgeColor', 'white');

if contrast == "AVN"
    labels = "Ido1";
else
    % Both genes are stated in the main text; labelling them makes the LPS
    % panel directly inspectable without implying that Ido1 is significant.
    labels = ["Ido1", "Tdo2"];
end
for symbol = labels
    marker = strcmpi(stats.Symbol, symbol);
    if any(marker)
        text(stats.Log2FC(marker), log10(stats.PlotMean(marker)+1), "  " + symbol, ...
            'FontSize', 8, 'Color', [0.55 0 0]);
    end
end
xline(0, ':', 'Color', [0.4 0.4 0.4]);
xlabel("log_2 fold change (" + contrast + " / control)");
ylabel('log_{10}(mean normalized count + 1)');
title("Intestine: " + contrast + " versus control", 'FontWeight', 'normal');
exportPanel(fig, "figS2_rnaseq_intestine_" + lower(contrast), p);
close(fig);
end

function exportPanel(fig, stem, p)
set(fig, 'Color', 'white', 'Renderer', 'painters');
axesHandles = findall(fig, 'Type', 'axes');
set(axesHandles, 'FontName', 'Arial', 'FontSize', 9, 'LineWidth', 0.8, ...
    'TickDir', 'out', 'Box', 'off');
pdfFile = fullfile(p.panels, stem + ".pdf");
pngFile = fullfile(p.panels, stem + ".png");
exportgraphics(fig, pdfFile, 'ContentType', 'vector', 'BackgroundColor', 'white');
exportgraphics(fig, pngFile, 'Resolution', 600, 'BackgroundColor', 'white');
exportgraphics(fig,fullfile(p.panels,stem + ".eps"),'ContentType','vector');
end

function sizeFactors = deseqSizeFactors(counts)
assert(all(counts >= 0, 'all') && all(mod(counts, 1) == 0, 'all'), ...
    'Raw counts must be nonnegative integers.');
usable = all(counts > 0, 2);
assert(any(usable), 'No genes are positive in every sample.');
geometricMeans = exp(mean(log(counts(usable, :)), 2));
ratios = counts(usable, :) ./ geometricMeans;
sizeFactors = median(ratios, 1, 'omitnan');
sizeFactors = sizeFactors ./ geomean(sizeFactors);
assert(all(isfinite(sizeFactors) & sizeFactors > 0), 'Invalid size factors.');
end

function q = bhAdjust(p)
q = nan(size(p));
valid = isfinite(p);
pv = p(valid);
if isempty(pv)
    return
end
[sorted, order] = sort(pv(:));
m = numel(sorted);
adjusted = sorted .* m ./ (1:m)';
adjusted = flipud(cummin(flipud(adjusted)));
adjusted = min(adjusted, 1);
unsorted = nan(m, 1);
unsorted(order) = adjusted;
q(valid) = reshape(unsorted, size(pv));
end
