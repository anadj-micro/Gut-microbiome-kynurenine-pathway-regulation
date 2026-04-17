%% Script to run lasso and select taxa associated with IDO1 expression
% run on ASVs, genera and families, standardize abundances with minmax function,
% remove samples with less than 1000 sequences and use CV=5

close all
clear all

% first load ASV table
asv = readtable('ASV_unstacked.txt', 'Delimiter', '\t');
% load the IDO1 table with samples of interest and subset asv table based
% on it
ido = readtable('ido_tdo_levels');
asv1 = asv(ismember(asv.Sample_ID, ido.Sample_ID) == 1,:); 

% check the depth of each sample
s = sum(asv1{:,2:end},2);
%mean(s)
min(s) %10198 sequences in the smalest sample

goodSamples = asv1.Sample_ID;   

%% Make a wide table containing abundances at all the taxonomical levels of
% interest (percentage)

abundTbl = readtable('ASV_unstacked_percentage.txt');
genus = readtable('genus_percentage_dada2.txt');
abundTbl = innerjoin(abundTbl, genus);
family = readtable('family_percentage_dada2.txt');
abundTbl = innerjoin(abundTbl, family);
asv = abundTbl;
asv1 = asv{:,2:end}./100; % divide every column with 100 to get relative abundances
asv{:, 2:end} = asv1;

% Select only those samples that have >1000 sequences
asv = asv(ismember(asv.Sample_ID, goodSamples) == 1,:);

% join IDO1 data with microbiome data
asv1 = innerjoin(asv, ido);
% rearrange
asv1 = asv1(:, [1 end-7:end 2:end-8]);

%% Normalize relative abundances by 16S counts

% read 16S table
BLoad = readtable('16S_qPCR', 'Delimiter', '\t');
BLoad.Properties.VariableNames{1} = 'Sample_ID';
% merge
tbl = innerjoin(asv1, BLoad(:,[1 end]));
tbl = tbl(:, [1:9 end 10:end-1]);

% normalize by multiplying
normMatrix = tbl{:, 11:end} .* tbl.x16S_qPCR_per_g;
tbl1 = tbl;
tbl1{:, 11:end} = normMatrix;
asv2 = tbl1;

%% Standardize each taxa with minmax function so that minimum becomes 0 
% and maximum 1

%asv2 = asv1; %run this if you haven't normalized relative abundances with
%16S counts

asv3 = normalize(asv2(:,11:end), 'range');
tblLasso = asv2;
tblLasso(:,11:end) = asv3;

%% Run lasso
nam = tblLasso.Properties.VariableNames(11:end);
opt = statset('UseParallel',true);
alpha = 1;
CV = 5;
Lambda = 100;
LambdaRatio = 1e-2;
filename = strcat('alpha', sprintf('%0.3f', alpha), '_', 'CV', sprintf('%0.0f', CV), '_', 'lambda', sprintf('%0.0f', Lambda), '_', 'lambdaRatio', sprintf('%0.0e', LambdaRatio));
rng(42) %set the seed
[B, FitInfo] = lasso(tblLasso{:, 11:end},tblLasso.ido1,'NumLambda',Lambda, 'LambdaRatio',LambdaRatio,'CV', CV ,'Options',opt, 'Alpha',alpha, 'MCReps', 100, 'PredictorNames', nam);

lassoPlot(B,FitInfo,'PlotType','Lambda','XScale','log');
filename1 = 'Trace_plot_coefficients_16S_normalized_minmax_scale';
filename2 = strcat(filename1, '_', filename, '_', '10_14_24');
print(figure(1), strcat(filename2, '.pdf'), '-dpdf')

lassoPlot(B,FitInfo,'PlotType','CV');
filename3 = 'Lambda_plot_16S_normalized_minmax_scale';
filename4 = strcat(filename3, '_', filename, '_', '10_14_24');
print(figure(2), strcat(filename4, '.pdf'), '-dpdf')

figure(3)
set(groot,'defaultAxesTickLabelInterpreter','none');
indx = FitInfo.Index1SE;
B0 = B(:,indx);
nonzeros = sum(B0 ~= 0)
predictors = find(B0);
PredNam = FitInfo.PredictorNames(predictors);
B1 = B0(predictors);
[~, isort] = sort(B1);
B1 = B1(isort);
PredNam = PredNam(isort);
set(gca, 'XTick', 1:length(PredNam), 'XTickLabel', PredNam);
hold on
plot([1 size(predictors,1)], [0 0], 'k--', 'LineWidth', 1)
ylabel('\beta', 'FontSize', 16, 'FontWeight', 'bold')
for i=1:length(B1)
    if B1(i) < 0
        plot([i i], [0 B1(i)], 'b-')
        hold on
        plot(i, B1(i),'bo','MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'MarkerSize', 10)
    else
        plot([i i], [0 B1(i)], 'r-')
        hold on
        plot(i, B1(i),'ro','MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'MarkerSize', 10)
    end
end
view(90,90)
box off

tbltaxonomy = readtable('tblASVtaxonomy_dada2.csv');
% replace ASV with taxonomy
PredNam1 = PredNam;
for i = 1:size(PredNam,2)
    g = PredNam{i};
    if contains(g, 'ASV') == 1 && contains (g, 'unclassified') ~= 1
        taxonomy = tbltaxonomy(strcmp(tbltaxonomy.ASV, g) == 1,:);
        newg = strcat(taxonomy.family, ';',taxonomy.genus, ';',taxonomy.ASV);
        PredNam{i} = cell2mat(newg);
    end
end
set(gca, 'XTick', 1:length(PredNam), 'XTickLabel', PredNam);

filename5 = 'Coefficientes_lasso_16S_normalized_minmax_scale';
filename6 = strcat(filename5, '_', filename, '_', '10_14_24');
print(figure(3), strcat(filename6, '.pdf'), '-dpdf')

%% Plot selected taxa
for i = 1:length(PredNam1)
    t = PredNam1{i};
    if contains(t, 'ASV') == 1
        tFull = contains(PredNam, t);
        t1 = PredNam{tFull};
    else
        t1 = t;
    end
    figure(i)
    scatter(tblLasso.ido1, tblLasso{:,t}, 100, 'ob','filled')
    xlabel('IDO1')
    ylabel(sprintf('Scaled abundance %s',t1), 'Interpreter', 'none')
    ylim([0 1])
    xlim([0 1.2])
    h = lsline;
    set(h,'color','g','LineWidth',3)
    ylim([0 1])
    xlim([0 1.2])
    box on
    print('-dpdf', sprintf('%s_ido1_relative_abundances_10_15_24.pdf', t));
end


%% Fit the model with lasso predictors to get R2

mdl1 = fitlm(tblLasso{:, 11:end},tblLasso.ido1, 'linear', 'PredictorVars',predictors)
% p < 0.0001, R2 = 0.807, R2 adjusted 0.7865 for relative abundances
% p < 0.0001, R2 = 0.699, R2 adjusted 0.6741 for 16S normalized abundances


% predict IDO1 levels based on the model
Xnew = tblLasso{:, 11:end};
Ypred = predict(mdl1, Xnew);

plot(tblLasso.ido1, Ypred, 'ko', 'MarkerFaceColor', 'black')

