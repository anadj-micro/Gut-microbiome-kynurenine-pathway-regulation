function test_revision_workflows
% Execute candidate code on a default path, and verify prespecified outputs.
here = fileparts(mfilename('fullpath'));
oldPath = path;
oldFolder = pwd;
cleanup = onCleanup(@() restore_environment(oldPath,oldFolder));
restoredefaultpath;
candidate = fileparts(here);
rnaFolder = fullfile(candidate,'RNA_sequencing_data_analysis');
mouseFolder = fullfile(candidate,'redoPcoAsAndCorrelationsFig3and4');
tissueFolder = fullfile(candidate,'Supporting_statistics');
addpath(rnaFolder,mouseFolder,tissueFolder);
cd(here);
diary(fullfile(here,'test_revision_workflows.log'));
diaryCleanup = onCleanup(@() diary('off'));

run_figureS2_rnaseq;
run_supporting_statistics;
tissueInput = fullfile(tissueFolder,'inputs','metabolomics2_sampleform.xlsx');
hasTissueInput = isfile(tissueInput);
if hasTissueInput
    run_tissue_cecal_interaction;
else
    fprintf('Tissue interaction skipped: original specimen-weight form must be supplied separately.\n');
end

%% Verify RNA-seq tables already recomputed by the candidate entry point.
results = fullfile(rnaFolder,'results');
genes = readtable(fullfile(results,'fig1_rnaseq_ltryptophan_catabolic_process_genes.csv'));
assert(height(genes)==10);
rna = readtable(fullfile(results,'fig1_rnaseq_gene_statistics_current_go.csv'),'TextType','string');
ido = rna(rna.Symbol=="Ido1" & rna.Organ=="Intestine" & rna.Contrast=="AVN",:);
tdo = rna(rna.Symbol=="Tdo2" & rna.Organ=="Intestine" & rna.Contrast=="LPS",:);
assert(abs(ido.FDR-2.53368847281026e-5)<1e-12);
assert(abs(tdo.FDR-.0118142610120652)<1e-12);

%% Verify mouse statistics, cohorts, archived coefficients and fixed errors.
folder = fullfile(mouseFolder,'statistics');
resolution = readtable(fullfile(folder,'figS4B_permanova_resolution.csv'));
assert(max(abs(resolution.PValue-[.1493;.0021;.0001;.0001;.0001;.0001]))<1e-12);
recovery = readtable(fullfile(folder,'recovery_permanova_dispersion.csv'));
assert(max(abs(recovery.PValue-[.0036;.0951]))<1e-12);
assert(max(abs(recovery.DispersionP-[.0087;.1020]))<1e-12);
load = readtable(fullfile(folder,'figS4A_load_ido1.csv'));
assert(abs(load.PearsonR-.695649717388531)<1e-12 && load.N==54);
family = readtable(fullfile(folder,'family_statistics_and_legacy_audit.csv'),'TextType','string');
assert(sum(family.PassesScreen)==42);
top = family(family.Family=="Lachnospiraceae",:);
assert(abs(top.PartialRho-.72731546972585)<1e-12);
assert(abs(top.HistoricalIncorrectQ-.000148222427772792)<1e-12);
assert(top.PartialQ<1e-7);
beta = readtable(fullfile(folder,'figS5_pooled_gavage_permanova_dispersion.csv'));
assert(beta.N==29 && abs(beta.PartialR2-.024626)<1e-6);
assert(abs(beta.PValue-.5216)<1e-12,'Archived seeded P not reproduced.');
sensitivity = readtable(fullfile(folder,'metabolite_sensitivity.csv'));
assert(max(abs(sensitivity.RawPearsonR-[.715362346975782;.710028071810167;.732609244226734]))<1e-12);
assert(max(abs(sensitivity.PartialP-[.0189804060471696;.611723963263424;.0451015034449827]))<1e-10);
assert(height(readtable(fullfile(folder,'metabolite_leave_one_out.csv')))==108);
if hasTissueInput
    interaction = readtable(fullfile(tissueFolder,'outputs','fig5G_interaction_likelihood_ratio.csv'));
    assert(abs(interaction.pValue(2)-9.4306e-8)<1e-12);
end
fprintf('ALL LOCAL NUMERICAL CHECKS PASSED. MATLAB %s\n',version);
writetable(table("PASS",hasTissueInput,string(version),datetime('now'), ...
    'VariableNames',{'NumericalChecks','TissueInteractionTested','MATLAB','Completed'}), ...
    fullfile(fileparts(mfilename('fullpath')),'test_results.csv'));
end

function restore_environment(oldPath,oldFolder)
path(oldPath);
cd(oldFolder);
end
