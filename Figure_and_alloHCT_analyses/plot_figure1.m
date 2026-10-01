% Reproduce Figure 1, Figure S1, and their allo-HCT association statistics.
% Run this script from any folder. Requires Statistics and Machine Learning
% Toolbox. Inputs are public original tables plus SILVA 138.1 assignments
% derived from public ASV sequences. No private clinical files are needed.
% Nothing is read from another revision folder or a precomputed result table.

%% 1. Read inputs and select the same 92 samples as the manuscript.
here = fileparts(mfilename('fullpath'));
input = fullfile(here, 'inputs');
output = fullfile(here, 'outputs');
if ~isfolder(output)
    mkdir(output);
end
metabolites = readtable(fullfile(input, 'tblMetabolitesUnstackedStoolPlasma.txt'), 'TextType', 'string');
coordinates = readtable(fullfile(input, 'current_taxumap_embedding.csv'), 'TextType', 'string');
coordinates.Properties.VariableNames = {'SampleID', 'TaxUMAP1', 'TaxUMAP2'};
counts = readtable(fullfile(input, 'tblcounts_asv_melt.csv'), 'TextType', 'string');
originalTaxonomy = readtable(fullfile(input, 'tblASVtaxonomy_silva132_v4v5_filter.csv'), 'TextType', 'string');
taxonomy = readtable(fullfile(input, 'asv_taxonomy_reclassified.csv'), 'TextType', 'string');
% Reclassification must preserve original ASV sequences and plotting colors.
[found, originalRows] = ismember(taxonomy.ASV, originalTaxonomy.ASV);
assert(all(found));
assert(all(taxonomy.Sequence == originalTaxonomy.Sequence(originalRows)));
assert(all(taxonomy.OldHexColor == originalTaxonomy.HexColor(originalRows)));
% Force identifiers to text: FMT patient IDs must not become numeric NaNs.
options = detectImportOptions(fullfile(input, 'tblASVsamples.csv'));
options = setvartype(options, {'SampleID', 'PatientID'}, 'string');
metadata = readtable(fullfile(input, 'tblASVsamples.csv'), options);
assert(numel(unique(metabolites.SampleID)) == height(metabolites));
assert(numel(unique(coordinates.SampleID)) == height(coordinates));
cohort = innerjoin(metabolites, coordinates, 'Keys', 'SampleID');
cohort = sortrows(cohort, 'SampleID');
% The obsolete Category (DIV/STR/VRE) field is not used in any analysis.
cohort.Category = [];
assert(height(cohort) == 92);
[found, rows] = ismember(cohort.SampleID, metadata.SampleID);
assert(all(found));
cohort.PatientID = string(metadata.PatientID(rows));
assert(numel(unique(cohort.PatientID)) == 73);
cohort.KYNTRP_stool = cohort.kynurenine_stool ./ cohort.tryptophan_stool;
cohort.KYNTRP_plasma = cohort.kynurenine_plasma ./ cohort.tryptophan_plasma;

%% 2. Calculate ASV relative abundances and inverse Simpson diversity.
% The denominator is every original read, including unclassified ASVs.
[inCohort, sampleRows] = ismember(counts.SampleID, cohort.SampleID);
c = counts(inCohort, :);
sampleRows = sampleRows(inCohort);
[found, asvColumns] = ismember(c.ASV, taxonomy.ASV);
assert(all(found), 'An observed ASV is missing from the revised taxonomy.');
assert(numel(unique(taxonomy.ASV)) == height(taxonomy));
assert(numel(unique(c.ASV)) == 2106);
asvCounts = accumarray([sampleRows, asvColumns], c.Count, [92, height(taxonomy)]);
readDepth = sum(asvCounts, 2);
asvRelative = asvCounts ./ readDepth;
cohort.InverseSimpson = 1 ./ sum(asvRelative.^2, 2);
assert(all(readDepth > 0));
assert(max(abs(sum(asvRelative, 2) - 1)) < 1e-12);

%% 3. Inherit the original atlas coordinates and dominant-ASV colors.
% Sort exactly as the original plot: ties retain original count-table order.
sortedCounts = sortrows(counts, {'SampleID', 'Count'}, {'ascend', 'descend'});
[~, first] = unique(sortedCounts.SampleID, 'stable');
dominant = sortedCounts(first, {'SampleID', 'ASV'});
[found, rows] = ismember(dominant.ASV, originalTaxonomy.ASV);
assert(all(found));
dominant.RGB = hex2rgb(originalTaxonomy.HexColor(rows));
atlas = innerjoin(coordinates, dominant, 'Keys', 'SampleID');
[found, rows] = ismember(cohort.SampleID, atlas.SampleID);
assert(all(found));
cohort.DominantASV = atlas.ASV(rows);
cohort.RGB = atlas.RGB(rows, :);

%% 4. Aggregate the frozen SILVA 138.1 assignments into families.
% Missing/uncultured calls stay in the denominator, not the named-family screen.
family = taxonomy.NewFamily;
missingFamily = ismissing(family) | ismember(lower(family), ...
    ["", "na", "<not present>", "unclassified", "uncultured", "unknown family"]);
family(missingFamily) = "Unclassified";
[families, ~, familyIndex] = unique(family);
familyPercent = zeros(92, numel(families));
familyRGB = zeros(numel(families), 3);
asvRGB = hex2rgb(taxonomy.OldHexColor);
weights = mean(asvRelative, 1)';
for j = 1:numel(families)
    members = familyIndex == j;
    familyPercent(:, j) = 100 * sum(asvRelative(:, members), 2);
    % Mean-abundance-weighted ASV colors reproduce the established palette.
    familyRGB(j, :) = sum(asvRGB(members, :) .* weights(members), 1) / sum(weights(members));
end
assert(max(abs(sum(familyPercent, 2) - 100)) < 1e-10);
prevalence = mean(familyPercent > 0, 1)';
meanPercent = mean(familyPercent, 1)';
keep = families ~= "Unclassified" & prevalence >= 0.20 & meanPercent >= 0.10;
assert(nnz(keep) == 29);

%% 5. Compute raw and neopterin-adjusted family Spearman correlations.
% All 92 samples are separate observations. There is no patient averaging
% or selection-stratum adjustment. Rank residuals implement partial Spearman.
% The archived family screen uses two-sided t-approximation P values.
x = tiedrank(familyPercent(:, keep));
y = tiedrank(cohort.KYNTRP_stool);
design = [ones(92, 1), tiedrank(cohort.neopterin_stool)];
rho = corr(x, y);
p = 2 * tcdf(-abs(rho) .* sqrt(90 ./ (1 - rho.^2)), 90);
partialRho = corr(x - design * (design \ x), y - design * (design \ y));
partialP = 2 * tcdf(-abs(partialRho) .* sqrt(89 ./ (1 - partialRho.^2)), 89);
familyStats = table(families(keep), prevalence(keep), meanPercent(keep), ...
    rho, p, bh(p), partialRho, partialP, bh(partialP), ...
    'VariableNames', {'Family', 'Prevalence', 'MeanPercent', 'SpearmanRho', ...
    'P', 'BHq', 'PartialSpearmanRho', 'PartialP', 'PartialBHq'});
familyStats.RGB = familyRGB(keep, :);
familyStats = sortrows(familyStats, 'SpearmanRho');
writetable(familyStats, fullfile(output, 'family_statistics.csv'));

%% 6. Recompute Figure 1B-D and Figure S1 statistics on unshifted data.
% Unlike the family-screen t approximation, these panels use MATLAB's
% corr(...,'Type','Spearman') P values, matching their original scripts.
enterobacteria = familyPercent(:, families == "Enterobacteriaceae");
predictors = {cohort.neopterin_plasma, cohort.neopterin_stool, ...
    cohort.InverseSimpson, enterobacteria, enterobacteria};
outcomes = {cohort.KYNTRP_plasma, cohort.KYNTRP_stool, ...
    cohort.KYNTRP_stool, cohort.tryptophan_stool, cohort.kynurenine_stool};
panel = ["1B"; "1C"; "1D"; "S1A"; "S1B"];
n = zeros(5, 1);
r = zeros(5, 1);
p = zeros(5, 1);
for j = 1:5
    complete = isfinite(predictors{j}) & isfinite(outcomes{j});
    n(j) = nnz(complete);
    [r(j), p(j)] = corr(predictors{j}(complete), outcomes{j}(complete), 'Type', 'Spearman');
end
panelStats = table(panel, n, r, p, 'VariableNames', {'Panel', 'N', 'SpearmanRho', 'P'});
panelStats.HolmAdjustedP = nan(5, 1);
[sortedP, order] = sort(p(4:5));
panelStats.HolmAdjustedP(3 + order) = min(1, cummax(sortedP .* [2; 1]));
writetable(panelStats, fullfile(output, 'panel_statistics.csv'));

%% 7. Reproduce only the two diversity sensitivities in the Figure 1 legend.
% These models retain all samples but use a patient random intercept.
% They are separate from the sample-level family tests above. No DIV/STR/VRE.
modelData = table(tiedrank(cohort.KYNTRP_stool), tiedrank(cohort.InverseSimpson), ...
    tiedrank(cohort.neopterin_stool), categorical(cohort.PatientID), ...
    'VariableNames', {'Outcome', 'Diversity', 'Neopterin', 'Patient'});
formulas = ["Outcome ~ Diversity + (1|Patient)"; ...
    "Outcome ~ Diversity + Neopterin + (1|Patient)"];
beta = zeros(2, 1);
modelP = zeros(2, 1);
for j = 1:2
    model = fitlme(modelData, formulas(j), 'FitMethod', 'REML');
    row = strcmp(model.Coefficients.Name, 'Diversity');
    beta(j) = model.Coefficients.Estimate(row);
    modelP(j) = model.Coefficients.pValue(row);
end
mixedStats = table(formulas, beta, modelP, 'VariableNames', {'Model', 'DiversityBeta', 'P'});
writetable(mixedStats, fullfile(output, 'diversity_sensitivity_statistics.csv'));
% Save the values needed to inspect joins, denominators, and sample inclusion.
sampleValues = cohort(:, {'SampleID', 'PatientID', 'DominantASV', ...
    'KYNTRP_stool', 'KYNTRP_plasma', 'neopterin_stool', 'neopterin_plasma', ...
    'tryptophan_stool', 'kynurenine_stool', 'InverseSimpson'});
sampleValues.ReadDepth = readDepth;
sampleValues.UnclassifiedPercent = familyPercent(:, families == "Unclassified");
sampleValues.EnterobacteriaceaePercent = enterobacteria;
writetable(sampleValues, fullfile(output, 'sample_values.csv'));

%% 8. Assemble Figure 1 using the established patient-only layout.
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 12 7.6]);
positions = [0.020 0.615 0.185 0.315; 0.370 0.625 0.170 0.268; ...
    0.575 0.625 0.170 0.268; 0.780 0.625 0.170 0.268; 0.060 0.250 0.910 0.290];
ax = axes(fig, 'Position', positions(1, :));
scatter(ax, atlas.TaxUMAP1, atlas.TaxUMAP2, 4, atlas.RGB + (1-atlas.RGB)*0.60, 'filled');
hold(ax, 'on');
scatter(ax, cohort.TaxUMAP1, cohort.TaxUMAP2, 30, cohort.RGB, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.55);
axis(ax, 'equal');
axis(ax, 'tight');
axis(ax, 'off');
title(ax, 'Allo-HCT microbiome atlas', 'FontSize', 8.5, 'FontWeight', 'normal');

% Keep the legend separate from the atlas, including the orange family.
key = axes(fig, 'Position', [0.207 0.620 0.130 0.300]);
hold(key, 'on');
axis(key, [0 1 0 1]);
axis(key, 'off');
labels = ["Clostridia", "Bacteroidota", "Actinomycetota", "Pseudomonadota", ...
    "Enterococcus", "Staphylococcus", "Streptococcus", "Lactobacillus", ...
    "Erysipelotrichaceae", "Other bacteria"];
colors = hex2rgb(["#BEA89A", "#1CF8EC", "#D0D0D0", "#EE2C2C", ...
    "#0D7E2B", "#F4EE26", "#AFCF3C", "#5875DE", "#FBA22E", "#CA0BE8"]);
text(key, 0, 0.98, 'Dominant ASV', 'FontWeight', 'bold', 'FontSize', 6.8);
for j = 1:numel(labels)
    yPosition = 0.85 - (j-1)*0.83/9;
    scatter(key, 0.09, yPosition, 24, colors(j, :), 's', 'filled');
    text(key, 0.19, yPosition, labels(j), 'FontSize', 5.6);
end

% Offsets are display-only and reproduce the existing Figure 1 panels.
displayX = {cohort.neopterin_plasma+1, cohort.neopterin_stool+10, cohort.InverseSimpson};
displayY = {cohort.KYNTRP_plasma+eps, cohort.KYNTRP_stool+1e-4, cohort.KYNTRP_stool+1e-4};
xlabels = ["Plasma neopterin (nM)", "Stool neopterin (nM)", "Inverse Simpson diversity"];
ylabels = ["Plasma KYN/TRP", "Stool KYN/TRP", "Stool KYN/TRP"];
for j = 1:3
    ax = axes(fig, 'Position', positions(j+1, :));
    scatter(ax, displayX{j}, displayY{j}, 27, cohort.RGB, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.45);
    set(ax, 'XScale', 'log', 'YScale', 'log');
    xlabel(ax, xlabels(j));
    ylabel(ax, ylabels(j));
    text(ax, 0.97, 0.04, sprintf('\\rho = %.2f; P = %.2g; n = %d', r(j), p(j), n(j)), ...
        'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontSize', 7);
    style_axes(ax, 7.2);
    axis(ax, 'square');
    if j == 3
        xlim(ax, [0.9 55]);
    end
end
ax = axes(fig, 'Position', positions(5, :));
b = bar(ax, familyStats.SpearmanRho, 0.72, 'FaceColor', 'flat', 'LineWidth', 0.3);
b.CData = familyStats.RGB;
hold(ax, 'on');
yline(ax, 0, 'k-');
significant = familyStats.BHq < 0.05;
scatter(ax, find(significant), familyStats.SpearmanRho(significant), 42, 'ko', 'LineWidth', 1.1);
set(ax, 'XTick', 1:height(familyStats), 'XTickLabel', familyStats.Family, ...
    'XTickLabelRotation', 90, 'TickLabelInterpreter', 'none');
ylabel(ax, 'Spearman \rho with stool KYN/TRP');
title(ax, 'Family associations with stool KYN/TRP (SILVA 138.1; 29 families)', 'FontWeight', 'normal', 'FontSize', 9);
text(ax, 0.995, 0.97, 'o  BH q < 0.05', 'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontSize', 7);
style_axes(ax, 6.5);
ylim(ax, [-0.33 0.50]);
xlim(ax, [0.25 29.75]);
for j = 1:5
    annotation(fig, 'textbox', [max(0.003, positions(j,1)-0.027), ...
        positions(j,2)+positions(j,4)+0.004, 0.035, 0.030], ...
        'String', char('A'+j-1), 'FontWeight', 'bold', 'FontSize', 12, 'EdgeColor', 'none', 'Margin', 0);
end
save_plot(fig, fullfile(output, 'Figure1'));

%% 9. Draw the Enterobacteriaceae ratio-component supplement.
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 9 4.6]);
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
titles = ["A   Stool tryptophan", "B   Stool kynurenine"];
for j = 1:2
    ax = nexttile(layout);
    y = outcomes{j+3};
    scatter(ax, log10(enterobacteria+0.01), log10(y+1), 35, cohort.RGB, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    xlabel(ax, 'log_{10}[Enterobacteriaceae (%) + 0.01]');
    ylabel(ax, 'log_{10}[concentration (nM) + 1]');
    title(ax, titles(j), 'FontWeight', 'normal');
    text(ax, 0.97, 0.97, sprintf('\\rho = %.3f; P = %.3g; n = %d', r(j+3), p(j+3), n(j+3)), ...
        'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    style_axes(ax, 10);
    xlim(ax, [-2.2 2]);
    ylim(ax, [min(log10(y+1))-0.2 max(log10(y+1))+0.65]);
    axis(ax, 'square');
end
save_plot(fig, fullfile(output, 'FigureS1'));
disp(panelStats);
disp(mixedStats);
disp(familyStats(ismember(familyStats.Family, ...
    ["Enterobacteriaceae", "Ruminococcaceae", "Oscillospiraceae", "Lachnospiraceae"]), :));
fprintf('Finished: %d samples, %d patients, %d screened families. MATLAB %s\n', ...
    height(cohort), numel(unique(cohort.PatientID)), height(familyStats), version);

function q = bh(p)
% Benjamini-Hochberg adjustment within one complete screen.
[sorted, order] = sort(p);
adjusted = sorted .* numel(p) ./ (1:numel(p))';
q = zeros(size(p));
q(order) = min(1, flipud(cummin(flipud(adjusted))));
end

function rgb = hex2rgb(hex)
% Convert the original six-digit taxonomy colors to MATLAB RGB values.
hex = erase(string(hex), '#');
hex = hex(:);
rgb = zeros(numel(hex), 3);
for j = 1:3
    rgb(:, j) = hex2dec(char(extractBetween(hex, 2*j-1, 2*j))) / 255;
end
end

function style_axes(ax, fontSize)
set(ax, 'FontName', 'Arial', 'FontSize', fontSize, 'Box', 'off', ...
    'TickDir', 'out', 'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.12);
end

function save_plot(fig, stem)
% Keep editable vector exports, plus a convenient image preview.
drawnow;
exportgraphics(fig, stem + ".png", 'Resolution', 300);
exportgraphics(fig, stem + ".pdf", 'ContentType', 'vector');
exportgraphics(fig, stem + ".eps", 'ContentType', 'vector');
end
