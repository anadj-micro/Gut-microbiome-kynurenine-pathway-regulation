% FIGURE 3 C-F: PCoAs, diversity associations, raw-to-partial family plot.
% Open this file in MATLAB and click Run. See README before editing samples.
here=fileparts(mfilename('fullpath'));
addpath(here);
% Input tables and calculated results remain visible in the workspace.
% Reset these structures so rerunning either script cannot retain old fields.
d=struct;
result=struct;
style=struct;
folder=fullfile(here,'inputs','figure3');

%% 1. Read IDO1, bacterial load, diversity, and ASV read counts.
ido=readtable(fullfile(folder,'ido_tdo_levels.txt'), ...
    'VariableNamingRule','preserve','TextType','string');
q=readtable(fullfile(folder,'16S_qPCR.txt'), ...
    'VariableNamingRule','preserve','TextType','string');
q.Properties.VariableNames{1}='Sample_ID';
diversity=readtable(fullfile(folder,'Diversity.txt'),'TextType','string');
counts=readtable(fullfile(folder,'ASV_unstacked.txt'), ...
    'VariableNamingRule','preserve','TextType','string');
reads=table(counts{:,1},sum(counts{:,2:end},2), ...
    'VariableNames',{'Sample_ID','Reads'});

%% 2. Keep matched mice with >=1,000 reads and usable measurements.
% This reproduces the current cohort; AVN batch selection awaits review.
minReads=1000;
m=innerjoin(ido,q(:,{'Sample_ID','16S_qPCR_per_g'}),'Keys','Sample_ID');
m=innerjoin(m,diversity,'Keys','Sample_ID');
m=innerjoin(m,reads,'Keys','Sample_ID');
keep=m.Reads>=minReads & isfinite(m.ido1) & isfinite(m.Simpson) & ...
    isfinite(m.('16S_qPCR_per_g')) & m.('16S_qPCR_per_g')>0;
m=m(keep,:);
d.samples=table(m.Sample_ID,m.('Sample group'),m.Experiment,m.ido1, ...
    log10(m.('16S_qPCR_per_g')),m.Simpson,m.Reads,'VariableNames', ...
    {'SampleID','Group','Experiment','IDO1','Log10Load','InverseSimpson','Reads'});

%% 3. Convert counts to percentages; read the original family export.
% Counts also supplied read depth, so no duplicate ASV-percentage file.
[x,d.asvNames]=match_abundance(counts,d.samples.SampleID);
d.asvPercent=100*x./sum(x,2);
family=readtable(fullfile(folder,'family_percentage_dada2.txt'), ...
    'VariableNamingRule','preserve','TextType','string');
[d.familyPercent,d.familyNames]=match_abundance(family,d.samples.SampleID);

assert(numel(unique(d.samples.SampleID))==height(d.samples),'Duplicate mouse IDs.');
assert(all(isfinite(d.samples.IDO1)),'Invalid IDO1 normalization.');

m=d.samples;
out=fullfile(here,'eps');

%% 4. Styling: change colors, labels, point size, and font size here.
% style.groups=["C","AVN","AVN_FMT","AVN_recovery","A","V","M","Cp"];
% style.labels=["Control","AVN","AVN + FMT","Recovery", ...
%     "Ampicillin","Vancomycin","Metronidazole","Ciprofloxacin"];
% style.colors=[.35 .35 .35; .78 .16 .14; .25 .60 .40; .48 .27 .67; ...
%     .20 .55 .85; .25 .35 .75; .90 .45 .15; .75 .55 .15];
% style.pointSize=42;
% style.fontSize=10;
% style.colorMap=turbo(256);
% style.controlCorner=[-1 1]; % orient controls left/up; axis signs are arbitrary
% rawColor=[.45 .45 .45];
% partialColor=[.10 .45 .70];
% highlightColors=[.82 .10 .10; .95 .55 .05]; % Lachno, Oscillo

%% 4. Styling: change colors, labels, point size, and font size here.
style.groups=["C","AVN","AVN_FMT","AVN_recovery","A","V","M","Cp"];
style.labels=["Control","AVN","AVN + FMT","Recovery", ...
    "Ampicillin","Vancomycin","Metronidazole","Ciprofloxacin"];
% Group colors (RGB converted from 0–255 to MATLAB 0–1 scale)
style.colors=[...
    0/255   0/255   255/255;   % Control: blue
    255/255 255/255 0/255;     % AVN: yellow
    255/255 255/255 0/255;     % AVN + FMT: yellow
    255/255 255/255 0/255;     % Recovery: yellow
    204/255 204/255 41/255;    % Ampicillin
    153/255 153/255 61/255;    % Vancomycin
    102/255 102/255 61/255;    % Metronidazole
    51/255  51/255  41/255];   % Ciprofloxacin
style.pointSize=100;
style.fontSize=16;
style.colorMap=turbo(256);
style.controlCorner=[-1 1]; % orient controls left/up; axis signs are arbitrary
% Plot background: RGB 230,230,230
style.backgroundColor=[230 230 230]/255;
% Border color
style.borderColor=[0 0 0];
% Border LineWidth (points). 0 = no border.
% C   AVN   FMT   Recovery   A   V   M   Cp
style.borderSize=[0; 0; 1; 2; 0; 0; 0; 0];
rawColor=[.45 .45 .45];
partialColor=[.10 .45 .70];
highlightColors=[.82 .10 .10; .95 .55 .05]; % Lachno, Oscillo

%% 5. Two views of the same unfiltered ASV Bray-Curtis coordinates.
assert(height(m)==54 && all(m.Reads>=1000),'Check Figure 3 sample inclusion.');
[result.pcoa,result.axisPercent]=draw_pcoa(d.asvPercent,m, ...
    style,out,["fig3C_pcoa_groups","fig3D_pcoa_IDO1"]);

%% 6. Rank adjustment shared by the family and diversity correlations.
% Spearman = correlation of ranks. Residualize ranked abundance and IDO1
% against ranked log10 load + experiment indicators, then correlate residuals.
experiment=dummyvar(categorical(m.Experiment));
design=[ones(height(m),1),tiedrank(m.Log10Load),experiment(:,2:end)];
yRank=tiedrank(m.IDO1);
yResidual=yRank-design*(design\yRank);

%% 7. Family screen: prevalence >=20% and mean abundance >=0.1%.
keep=mean(d.familyPercent>0,1)>=.20 & mean(d.familyPercent,1)>=.10;
x=d.familyPercent(:,keep);
names=d.familyNames(keep)';
assert(numel(names)==42,'Expected 42 families in the unchanged cohort.');
raw=corr(x,m.IDO1,'Type','Spearman');
xRank=tiedrank(x);
xResidual=xRank-design*(design\xRank);
partial=corr(xResidual,yResidual); % Pearson on residual ranks = partial Spearman
[partial,order]=sort(partial,'descend');
raw=raw(order);
names=names(order);
result.family=table(names,raw,partial,'VariableNames', ...
    {'Family','RawSpearman','PartialSpearman'});

%% 8. Horizontal connected-dot plot, ordered by partial correlation.
fig=figure('Color',[230 230 230]/255,'Units','inches','Position',[1 .3 11 20.36]);   % <-- CHANGED (height)
ax=axes(fig,'Position',[.04 .10 .45 .85]);   % <-- CHANGED (mirrored horizontal position)
ax.Color=[230 230 230]/255;   % <-- CHANGED (axes background to match)
hold(ax,'on');
y=(1:numel(names))';
rawColor=[82 3 252]/255;       % <-- CHANGED (green for raw Spearman)
partialColor=[3 219 252]/255;   % <-- CHANGED (blue for partial)
change=plot(ax,[raw partial]',[y y]','-','Color',[.3 .3 .3],'LineWidth',.8);
hRaw=scatter(ax,raw,y,50,rawColor,'filled');
hPartial=scatter(ax,partial,y,60,partialColor,'filled');
for k=1:2
    target=["Lachnospiraceae","Oscillospiraceae"];
    use=names==target(k);
    scatter(ax,partial(use),y(use),80,highlightColors(k,:),'filled', ...
        'MarkerEdgeColor','k','LineWidth',.7);
end
xline(ax,0,'k-','HandleVisibility','off');
set(ax,'YTick',y,'YTickLabel',strrep(names,'_',' '),'YDir','reverse', ...
    'TickLabelInterpreter','none','XTick',-1:.5:1);
ax.YAxisLocation='right';   % <-- CHANGED (moves y tick labels to right side)
xlim(ax,[-1 1]);
ylim(ax,[0 numel(names)+1]);
box(ax,'on');   % <-- CHANGED (added)
xlabel(ax,'Spearman association with colonic IDO1');
title(ax,{'42-family associations: raw versus', ...
    'load/experiment-adjusted'},'FontWeight','bold');
legend(ax,[change(1) hRaw hPartial],{'Raw-to-partial change','Raw','Partial'}, ...
    'Location','southoutside','Orientation','horizontal','Box','off','FontSize',16);
export_eps(fig,out,'fig3F_family_raw_partial',style.fontSize);

%% 9. Inverse Simpson: raw and adjusted scatter panels (Figure 3E).
divRank=tiedrank(m.InverseSimpson);
divResidual=divRank-design*(design\divRank);
[rho,p]=corr(m.InverseSimpson,m.IDO1,'Type','Spearman');
partialR=corr(divResidual,yResidual);
df=height(m)-rank(design)-1;
partialP=2*tcdf(-abs(partialR)*sqrt(df/(1-partialR^2)),df);
result.diversity=[rho p partialR partialP];
xx={m.InverseSimpson,divResidual};
yy={m.IDO1,yResidual};
for panel=1:2
    fig=figure('Color','w','Units','inches','Position',[1 1 7.8 5.8]);
    ax=axes(fig);
    hold(ax,'on');
    if panel==1
        h=gobjects(numel(style.groups),1);
        for g=1:numel(style.groups)
            use=m.Group==style.groups(g);
            h(g)=scatter(ax,xx{panel}(use),yy{panel}(use),style.pointSize, ...
                style.colors(g,:),'filled','MarkerEdgeColor','k','LineWidth',.4);
        end
        legend(ax,h,style.labels,'Location','eastoutside','Box','off', ...
            'Interpreter','none','FontSize',style.fontSize-1);
        xlabel(ax,'Inverse Simpson diversity');
        ylabel(ax,'Colonic IDO1');
        label=sprintf('Spearman \\rho = %.3f; P = %.3f',rho,p);
        stem='fig3E_diversity_raw';
    else
        scatter(ax,xx{panel},yy{panel},style.pointSize,partialColor,'filled');
        xlabel(ax,'Inverse Simpson residual rank');
        ylabel(ax,'IDO1 residual rank');
        label=sprintf('Partial Spearman \\rho = %.3f; P = %.3f',partialR,partialP);
        stem='fig3E_diversity_partial';
    end
    % Straight lines aid visualization only; the annotations use rank tests.
    lineX=linspace(min(xx{panel}),max(xx{panel}),100);
    plot(ax,lineX,polyval(polyfit(xx{panel},yy{panel},1),lineX),'k-', ...
        'LineWidth',1,'HandleVisibility','off');
    % Keep the statistics above the axes so no mouse is hidden by a label.
    title(ax,label,'FontWeight','normal','FontSize',style.fontSize);
    axis(ax,'square');
    grid(ax,'on');
    export_eps(fig,out,stem,style.fontSize);
end
disp(result.family(1:2,:));
fprintf('Figure 3: 54 mice; 42 families; five EPS files saved in %s\n',out);
