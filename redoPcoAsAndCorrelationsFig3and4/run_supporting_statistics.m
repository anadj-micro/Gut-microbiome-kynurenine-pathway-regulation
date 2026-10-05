% Supporting mouse analyses using the SAME original inputs as Figures 3–4.
% Run from any directory. No manuscript, private data, or archived-result input.
% Permutation procedures preserve the archived algorithms, seeds, and cohorts.
here = fileparts(mfilename('fullpath'));
addpath(here);
statsDir = fullfile(here,'statistics');
if ~isfolder(statsDir)
    mkdir(statsDir);
end

%% 1. Reuse Figure 3's original-input matching, then save the 54-mouse cohort.
plot_figure3;
statsDir = fullfile(here,'statistics');
writetable(m,fullfile(statsDir,'figure3_54_mouse_manifest.csv'));
[loadR,loadP] = corr(m.Log10Load,m.IDO1,'Type','Pearson');
writetable(table(loadR,loadR^2,loadP,height(m),'VariableNames', ...
    {'PearsonR','R2','P','N'}),fullfile(statsDir,'figS4A_load_ido1.csv'));
fig = figure('Visible','off','Color','w');
scatter(m.Log10Load,m.IDO1,35,'filled');
hold on;
lineX = linspace(min(m.Log10Load),max(m.Log10Load),100);
plot(lineX,polyval(polyfit(m.Log10Load,m.IDO1,1),lineX),'k-');
xlabel('log_{10} fecal 16S copies/g');
ylabel('Colonic Ido1');
title(sprintf('Pearson r = %.3f; R^2 = %.3f',loadR,loadR^2));
export_eps(fig,out,'figS4A_load_ido1',12);

%% 2. Family P and BH q: raw Spearman and load/experiment-adjusted ranks.
% The correct partial-correlation t statistic divides by sqrt(1-rho^2).
% Save the historical denominator error separately; do not perpetuate it.
rawR = corr(d.familyPercent,m.IDO1,'Type','Spearman');
% Use a stable two-sided t-tail approximation for raw Spearman tests too;
% MATLAB's corr can round extremely small Spearman P values to zero.
rawP = 2*tcdf(-abs(rawR).*sqrt((height(m)-2)./max(eps,1-rawR.^2)),height(m)-2);
familyRanks = tiedrank(d.familyPercent);
familyResidual = familyRanks-design*(design\familyRanks);
partialR = corr(familyResidual,yResidual);
df = height(m)-rank(design)-1;
partialP = 2*tcdf(-abs(partialR).*sqrt(df./max(eps,1-partialR.^2)),df);
screen = mean(d.familyPercent>0,1)>=.20 & mean(d.familyPercent,1)>=.10;
rawQ = nan(size(rawP));
partialQ = nan(size(partialP));
rawQ(screen) = bh_fdr(rawP(screen));
partialQ(screen) = bh_fdr(partialP(screen));
historicalP = 2*tcdf(-abs(partialR).*sqrt(df./max(1-eps,1-partialR.^2)),df);
historicalQ = nan(size(partialP));
historicalQ(screen) = bh_fdr(historicalP(screen));
families = table(d.familyNames(:),mean(d.familyPercent,1)', ...
    100*mean(d.familyPercent>0,1)',screen',rawR,rawP,rawQ,partialR, ...
    partialP,partialQ,historicalP,historicalQ,'VariableNames', ...
    {'Family','MeanAbundancePercent','PrevalencePercent','PassesScreen', ...
    'RawRho','RawP','RawQ','PartialRho','PartialP','PartialQ', ...
    'HistoricalIncorrectP','HistoricalIncorrectQ'});
families = sortrows(families,'PartialRho','descend');
writetable(families,fullfile(statsDir,'family_statistics_and_legacy_audit.csv'));
disp(families(ismember(families.Family, ...
    ["Lachnospiraceae","Oscillospiraceae","Enterobacteriaceae"]),:));
assert(sum(screen)==42);

%% 3. Six prespecified resolutions; prevalence >=20%, no mean-abundance filter.
% Adjust for log10 load and experiment; permute residualized IDO1 within
% experiment, using the existing archived permanova_term implementation.
taxonomy = readtable(fullfile(folder,'tblASVtaxonomy_dada2.csv'), ...
    'TextType','string','VariableNamingRule','preserve');
[matched,loc] = ismember(d.asvNames,string(taxonomy.ASV));
assert(all(matched),'Mouse taxonomy must cover every ASV.');
levels = ["phylum","class","order","family","genus","asv"];
resolution = table;
dummy = dummyvar(categorical(m.Experiment));
reduced = [ones(height(m),1),zscore(m.Log10Load),dummy(:,2:end)];
for k = 1:numel(levels)
    if levels(k)=="asv"
        abundance = d.asvPercent;
    else
        labels = taxonomy.(levels(k))(loc);
        labels(ismissing(labels) | labels=="" | labels=="<not present>") = "Unclassified";
        [names,~,index] = unique(labels);
        abundance = zeros(height(m),numel(names));
        for j = 1:numel(names)
            abundance(:,j) = sum(d.asvPercent(:,index==j),2);
        end
    end
    retained = mean(abundance>0,1)>=.20 & sum(abundance,1)>0;
    relative = abundance(:,retained)/100;
    D = squareform(pdist(relative,@bray_curtis_distance));
    test = permanova_term(D,reduced,zscore(m.IDO1),9999,700+k,m.Experiment);
    test = addvars(test,levels(k),sum(retained),height(m),'Before',1, ...
        'NewVariableNames',{'TaxonomicLevel','FeaturesRetained','Samples'});
    resolution = [resolution;test];
end
writetable(resolution,fullfile(statsDir,'figS4B_permanova_resolution.csv'));
disp(resolution(:,{'TaxonomicLevel','PartialR2','PValue'}));
fig = figure('Visible','off','Color','w');
bar(resolution.PartialR2);
set(gca,'XTick',1:6,'XTickLabel',levels);
ylabel('Partial R^2 for Ido1');
for k = 1:6
    text(k,resolution.PartialR2(k),sprintf('P=%.4g',resolution.PValue(k)), ...
        'VerticalAlignment','bottom','HorizontalAlignment','center');
end
ylim([0 .23]);
export_eps(fig,out,'figS4B_permanova_resolution',12);

%% 4. Two recovery contrasts on unfiltered ASV relative abundance.
% Dispersion reproduces the archived positive-PCoA-axis centroid permutation
% test (not a new bias-corrected or negative-eigenvalue PERMDISP method).
recovery = table;
experiments = [7 8];
groupA = ["AVN","C"];
groupB = ["AVN_FMT","AVN_recovery"];
for k = 1:2
    use = double(m.Experiment)==experiments(k) & ismember(m.Group,[groupA(k),groupB(k)]);
    group = double(m.Group(use)==groupA(k));
    D = squareform(pdist(d.asvPercent(use,:),@bray_curtis_distance));
    test = group_test(D,group,820+k,920+k);
    test = addvars(test,experiments(k),groupA(k),sum(group==1),groupB(k),sum(group==0), ...
        'Before',1,'NewVariableNames',{'Experiment','Group1','N1','Group2','N2'});
    recovery = [recovery;test];
end
writetable(recovery,fullfile(statsDir,'recovery_permanova_dispersion.csv'));
disp(recovery);

%% 5. Read Figure 4 independently; primary pooled-mixture versus AVN test.
plot_figure4;
statsDir = fullfile(here,'statistics');
pc = d.samples(d.pcoaRows,:);
% Restore the workbook order used by the archived seeded test. PCoA display
% uses sorted IDs; order does not change R2, but changes sampled permutations.
[matched,order] = ismember(original.Sample_ID,pc.SampleID);
assert(all(matched));
pc = pc(order,:);
asv = d.asvPercent(order,:);
use = ismember(pc.Group,["AVN","AVN_P"]);
assert(sum(use)==29);
group = double(pc.Group(use)=="AVN_P");
D = squareform(pdist(asv(use,:),@bray_curtis_distance));
test = group_test(D,group,461,963);
writetable(test,fullfile(statsDir,'figS5_pooled_gavage_permanova_dispersion.csv'));
disp(test);

%% 6. Sensitivity analyses retain exactly Figure 4's 36 mice and normalization.
% Categorical experiment + treatment group; not responder status. No ranks.
m = d.samples(d.metabolomicsRows,:);
dummyExperiment = dummyvar(categorical(m.Experiment));
dummyGroup = dummyvar(categorical(m.Group));
Z = [dummyExperiment(:,2:end),dummyGroup(:,2:end)];
X = [ones(height(m),1),Z];
assert(rank(X)==size(X,2));
ids = [888;930;948];
[matched,rows] = ismember(ids,d.ions.IonID);
assert(all(matched));
summary = table;
leaveOneOut = table;
for k = 1:3
    peak = d.peakArea(:,rows(k));
    [r,p] = corr(m.IDO1,peak,'Type','Pearson');
    looR = zeros(height(m),1);
    looP = zeros(height(m),1);
    for j = 1:height(m)
        keep = (1:height(m))'~=j;
        [looR(j),looP(j)] = corr(m.IDO1(keep),peak(keep),'Type','Pearson');
    end
    [partialR,partialP] = partialcorr(m.IDO1,peak,Z,'Type','Pearson');
    residualX = m.IDO1-X*(X\m.IDO1);
    residualY = peak-X*(X\peak);
    checkR = corr(residualX,residualY);
    df = height(m)-rank(X)-1;
    checkP = 2*tcdf(-abs(checkR)*sqrt(df/(1-checkR^2)),df);
    assert(abs(partialR-checkR)<1e-12 && abs(partialP-checkP)<1e-10);
    summary = [summary;table(ids(k),r,p,min(looR),max(looR),partialR,partialP,df, ...
        'VariableNames',{'IonID','RawPearsonR','RawP','LOOMinR','LOOMaxR', ...
        'PartialPearsonR','PartialP','PartialDF'})];
    leaveOneOut = [leaveOneOut;table(repmat(ids(k),height(m),1),m.SampleID,looR,looP, ...
        'VariableNames',{'IonID','OmittedMouse','PearsonR','P'})];
end
writetable(summary,fullfile(statsDir,'metabolite_sensitivity.csv'));
writetable(leaveOneOut,fullfile(statsDir,'metabolite_leave_one_out.csv'));
writetable(m,fullfile(statsDir,'metabolite_36_mouse_manifest.csv'));
disp(summary);

function test = group_test(D,group,permanovaSeed,dispersionSeed)
test = permanova_term(D,ones(numel(group),1),group,9999,permanovaSeed,[]);
[coordinates,eigenvalues] = cmdscale(D);
coordinates = coordinates(:,1:min(size(coordinates,2),sum(eigenvalues>0)));
observed = dispersion_statistic(coordinates,group);
rng(dispersionSeed,'twister');
exceedances = 0;
for k = 1:9999
    permuted = group(randperm(numel(group)));
    exceedances = exceedances + (dispersion_statistic(coordinates,permuted)>=observed);
end
test.DispersionP = (1+exceedances)/10000;
test.DispersionSeed = dispersionSeed;
test.PERMANOVASeed = permanovaSeed;
test.N = numel(group);
end

function statistic = dispersion_statistic(coordinates,group)
means = zeros(2,1);
for k = 0:1
    points = coordinates(group==k,:);
    means(k+1) = mean(sqrt(sum((points-mean(points,1)).^2,2)));
end
statistic = abs(diff(means));
end
