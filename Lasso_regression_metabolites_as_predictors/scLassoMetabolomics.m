%% SCript to run lasso on metabolomics data to select for predictors of 
% IDO1 expression levels
close all
clear all

% read NORMALIZED data table
data = readtable('DATA_NORM.xlsx','Sheet','ions');

% subset samples of interest
metadata = readtable('comment_samples.txt', 'Delimiter', '\t');

% define groups
metadata.Treatment = strings(size(metadata,1),1);
for i = 1:size(metadata,1)
    sample = metadata.dsSampleCode{i};
    if contains(sample, '_AVN_PM') == 1
        metadata{strcmp(metadata.dsSampleCode, sample) == 1, "Treatment"} = "AVNProb";
    elseif contains(sample, '_AVN_Kleb') == 1
        metadata{strcmp(metadata.dsSampleCode, sample) == 1, "Treatment"} = "AVNKleb";
    elseif contains(sample, '_AVN') == 1
        metadata{strcmp(metadata.dsSampleCode, sample) == 1, "Treatment"} = "AVN";
    elseif contains(sample, 'C') == 1
        metadata{strcmp(metadata.dsSampleCode, sample) == 1, "Treatment"} = "C";
    end
end

% focus on samples you want to compare
samples = unique(metadata(strcmp(metadata.Treatment, 'AVNProb') == 1 | strcmp(metadata.Treatment, 'AVN') == 1 | strcmp(metadata.Treatment, 'C') == 1, :));
samples1 = samples(strcmp(samples.housing, 'single') == 1, :);
samples1 = samples1(strcmp(samples1.comment, 'outlier') ~= 1, :);

% subset metabolites table for samples of interest
%transpose data
dataT = rows2vars(data(:, [1 7:end]), 'VariableNamesSource', 'ionIdx','VariableNamingRule','preserve');
dataT.Properties.VariableNames{1} = 'SampleID';
dataT1 = dataT(ismember(dataT.SampleID, samples1.dsSampleCode) == 1,:);
dataT = dataT1;

% get information about IDO1 levels
metadata = readtable('DATA_NORM.xlsx','Sheet','samples');
metadata1 = metadata(:, 4:16);
metadata2 = unique(metadata1);
metadata2.Properties.VariableNames{1} = 'SampleID';

% join 2 tables
dataT1 = innerjoin(dataT, metadata2(:,{'SampleID','comment','ido1_expression_level'}));
dataT1 = dataT1(:, [1 end-1 end 2:end-2]);

%% Standardize each taxa with minmax function so that minimum becomes 0 
% and maximum 1

dataT2 = normalize(dataT1(:,4:end), 'range');
tblLasso = dataT1;
tblLasso(:,4:end) = dataT2;

%% Run lasso
nam = tblLasso.Properties.VariableNames(4:end);
opt = statset('UseParallel',true);
alpha = 1;
CV = 5;
Lambda = 100;
LambdaRatio = 1e-2;
filename = strcat('alpha', sprintf('%0.3f', alpha), '_', 'CV', sprintf('%0.0f', CV), '_', 'lambda', sprintf('%0.0f', Lambda), '_', 'lambdaRatio', sprintf('%0.0e', LambdaRatio));
rng(42) %set the seed
[B, FitInfo] = lasso(tblLasso{:, 4:end},tblLasso.ido1_expression_level,'NumLambda',Lambda, 'LambdaRatio',LambdaRatio,'CV', CV ,'Options',opt, 'Alpha',alpha, 'MCReps', 100, 'PredictorNames', nam);

lassoPlot(B,FitInfo,'PlotType','Lambda','XScale','log');
filename1 = 'Trace_plot_coefficients_metabolomics_minmax_scale';
filename2 = strcat(filename1, '_', filename, '_', '10_17_24');
print(figure(1), strcat(filename2, '.pdf'), '-dpdf')

lassoPlot(B,FitInfo,'PlotType','CV');
filename3 = 'Lambda_plot_metabolomics_minmax_scale';
filename4 = strcat(filename3, '_', filename, '_', '10_17_24');
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

% replace ion with top classification
PredNam1 = PredNam;
for i = 1:size(PredNam,2)
    g = str2double(PredNam{i});
    newg = data{data.ionIdx == g, 'ionTopName'};
    PredNam{i} = cell2mat(newg);
end
set(gca, 'XTick', 1:length(PredNam), 'XTickLabel', PredNam);

filename5 = 'Coefficientes_lasso_metabolomics_minmax_scale';
filename6 = strcat(filename5, '_', filename, '_', '10_17_24');
print(figure(3), strcat(filename6, '.pdf'), '-dpdf')

%% Plot selected taxa
for i = 1:length(PredNam1)
    t = str2double(PredNam1{i});
    t1 = cell2mat(data{data.ionIdx == t, 'ionTopName'});
    figure(i)
    scatter(tblLasso.ido1_expression_level, tblLasso{:,t}, 100, 'ob','filled')
    xlabel('IDO1')
    ylabel(sprintf('Scaled abundance %s',t1), 'Interpreter', 'none')
    ylim([0 1])
    xlim([0 1.2])
    h = lsline;
    set(h,'color','g','LineWidth',3)
    ylim([0 1])
    xlim([0 1.2])
    box on
    print('-dpdf', sprintf('%d_ido1_10_17_24.pdf', t));
end
    
