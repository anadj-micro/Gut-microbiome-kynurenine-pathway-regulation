% Figure 5G: treatment-by-compartment interaction, mouse random intercept.
% The public ion workbook lacks specimen weights. Original sample form is
% therefore required, and must be added to the deposit for reproducibility.
here = fileparts(mfilename('fullpath'));
folder = fullfile(here,'inputs');
out = fullfile(here,'outputs');
if ~isfolder(out)
    mkdir(out);
end
ions = readtable(fullfile(folder,'DATA_LOESS_NORM_filtered_cecal_content_and_tissue.xlsx'), ...
    'Sheet','ions','VariableNamingRule','preserve','TextType','string');
samples = readtable(fullfile(folder,'metabolomics2_sampleform.xlsx'), ...
    'VariableNamingRule','preserve','TextType','string');
row = find(strcmpi(ions.ionTopName,'Gamma-Tocotrienol'));
assert(numel(row)==1);
ids = string(ions.Properties.VariableNames(8:end))';
intensity = ions{row,8:end}';
samples = samples(startsWith(samples.name,'KP_'),:);
[matched,loc] = ismember(samples.name,ids);
assert(all(matched));
samples.Intensity = intensity(loc);
samples.Weight = samples.('weight (g)');
assert(all(samples.Weight>0) && all(samples.Intensity>0));
samples.Treatment = categorical(samples.treatment,["C","AVN"]);
samples.Compartment = categorical(samples.type,["cecal_content","tissue"]);
samples.AnimalID = categorical(samples.treatment+"_"+string(samples.replicate));
samples.Log2PeakAreaPerGram = log2(samples.Intensity./samples.Weight);
assert(height(samples)==15 && numel(unique(samples.AnimalID))==8);
additive = fitlme(samples,'Log2PeakAreaPerGram ~ Treatment+Compartment+(1|AnimalID)', ...
    'FitMethod','ML');
full = fitlme(samples,'Log2PeakAreaPerGram ~ Treatment*Compartment+(1|AnimalID)', ...
    'FitMethod','ML');
comparison = dataset2table(compare(additive,full));
writetable(samples(:,{'name','AnimalID','Treatment','Compartment','Weight', ...
    'Intensity','Log2PeakAreaPerGram'}),fullfile(out,'fig5G_tissue_cecal_values.csv'));
writetable(comparison,fullfile(out,'fig5G_interaction_likelihood_ratio.csv'));
writetable(dataset2table(full.Coefficients),fullfile(out,'fig5G_mixed_model_coefficients.csv'));
disp(comparison);
