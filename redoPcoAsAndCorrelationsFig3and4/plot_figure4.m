% FIGURE 4 E-I and supporting S4A: ion correlations and recolonization PCoA.
% Open this file in MATLAB and click Run. All inputs are in this folder.
here=fileparts(mfilename('fullpath'));
addpath(here);
% Input tables and calculated results remain visible in the workspace.
% Reset these structures so rerunning either script cannot retain old fields.
d=struct;
result=struct;
style=struct;
folder=fullfile(here,'inputs','figure4');

%% 1. Start with Ana's original 43-mouse workbook and add five AVN+K mice.
original=readtable(fullfile(folder, ...
    'ido_tdo_levels_10_16_24_original_43_normalized.xlsx'), ...
    'VariableNamingRule','preserve','TextType','string');
m=table(original.Sample_ID,original.('Sample group'),original.Experiment, ...
    original.ido1,'VariableNames',{'SampleID','Group','Experiment','RawIDO1'});
metadata=readtable(fullfile(folder,'DATA_NORM.xlsx'),'Sheet','samples', ...
    'VariableNamingRule','preserve','TextType','string');
k=metadata(metadata.experiment==15 & metadata.dsUser1=="AVNKleb",:);
k=k(:,{'dsSampleCode','experiment','ido1_expression_level'});
% Injection repeats are not extra mice. Reject conflicting IDO1 records.
k=unique(k,'rows','stable');
assert(height(k)==numel(unique(k.dsSampleCode)),'Conflicting AVN+K records.');
m=[m; table(k.dsSampleCode,repmat("AVN_K",height(k),1), ...
    k.experiment,k.ido1_expression_level,'VariableNames',m.Properties.VariableNames)];

%% 2. Normalize IDO1 to each experiment's control mean and label responders.
% Recalculate from raw IDO1; do not use the workbook's normalized column.
m.IDO1=zeros(height(m),1);
for e=unique(m.Experiment)'
    use=m.Experiment==e;
    control=use & m.Group=="C";
    assert(any(control),'An experiment has no control for normalization.');
    m.IDO1(use)=m.RawIDO1(use)/mean(m.RawIDO1(control));
end
% Ana's descriptive rule uses the FULL cohort, not the metabolomics subset.
cutoff=max(m.IDO1(m.Group=="AVN"));
m.Responder=m.Group=="AVN_P" & m.IDO1>cutoff;
d.samples=removevars(m,'RawIDO1');

%% 3. Match the metabolomics columns to these mice (36 matches).
ions=readtable(fullfile(folder,'DATA_NORM.xlsx'),'Sheet','ions', ...
    'VariableNamingRule','preserve','TextType','string');
[hasIon,columns]=ismember(m.SampleID,string(ions.Properties.VariableNames(7:end)));
d.metabolomicsRows=find(hasIon);
d.peakArea=ions{:,6+columns(hasIon)}';
assert(all(isfinite(d.peakArea),'all'),'Missing ion measurements.');
d.ions=table(ions.ionIdx,ions.ionTopName,'VariableNames',{'IonID','Annotation'});
d.ions.Annotation(ismissing(d.ions.Annotation))="";

%% 4. Supporting PCoA: C/AVN/AVN+P only, sorted by mouse ID (43 mice).
rows=find(m.Group~="AVN_K");
[~,order]=sort(m.SampleID(rows));
d.pcoaRows=rows(order);
asv=readtable(fullfile(folder,'ASV_unstacked_percentage_dada2_filtered.txt'), ...
    'VariableNamingRule','preserve','TextType','string');
[d.asvPercent,d.asvNames]=match_abundance(asv,m.SampleID(d.pcoaRows));
assert(numel(unique(d.samples.SampleID))==height(d.samples),'Duplicate mouse IDs.');
assert(all(isfinite(d.samples.IDO1)),'Invalid IDO1 normalization.');

out=fullfile(here,'eps');

%% 5. Styling: same colors for PCoAs and metabolite scatter panels.
style.groups=["C","AVN","AVN_K","AVN_P"];
style.labels=["Control","AVN","AVN + Klebsiella","AVN + mixture"];
style.colors=[.36 .36 .36; .80 .19 .16; .15 .52 .70; .48 .28 .68];
style.responderEdge=[.98 .48 .08];
style.pointSize=46;
style.fontSize=10;
style.colorMap=parula(256);
style.controlCorner=[1 -1]; % orient controls right/down; axis signs are arbitrary
candidateIDs=[888 930 948 69];
candidateLabels={'Gamma-tocotrienol','C29:3', ...
    sprintf('4a-Carboxy-4b-methyl-5a-\ncholesta-8,24-dien-3b-ol'),'Proline'};
candidateColors=[.65 .05 .08; .95 .56 .05; .50 .22 .62; .06 .32 .58];

%% 6. Supplementary PCoA: preserve the saved 43-mouse C/AVN/AVN+P subset.
% The 48-mouse main cohort includes five additional AVN+K mice. They were
% not in the saved S4A ordination. Do not relabel this as a 48-mouse panel.
pcStyle=style;
pcStyle.groups=style.groups([1 2 4]);
pcStyle.labels=style.labels([1 2 4]);
pcStyle.colors=style.colors([1 2 4],:);
% Recompute ASV distances and coordinates; no PERMANOVA label is displayed.
[result.pcoa,result.axisPercent]=draw_pcoa(d.asvPercent, ...
    d.samples(d.pcoaRows,:),pcStyle,out, ...
    ["figS4A_pcoa_groups","figS4A_pcoa_IDO1"]);

%% 7. Recompute the 1,106-ion Pearson screen in the 36 matched mice.
m=d.samples(d.metabolomicsRows,:);
x=m.IDO1;
assert(height(m)==36 && size(d.peakArea,2)==1106);
assert(sum(m.Responder)==3,'Expected three historical responders in this subset.');
[r,p]=corr(d.peakArea,x,'Type','Pearson');
q=bh_fdr(p);
[rSorted,order]=sort(r);
qSorted=q(order);
result.ions=table(d.ions.IonID,r,p,q,'VariableNames',{'IonID','PearsonR','P','Q'});
[yes,candidateRows]=ismember(candidateIDs,d.ions.IonID);
assert(all(yes));
position=zeros(numel(r),1);
position(order)=1:numel(r);

%% 8. All-ion waterfall with BH-FDR colors and the four selected ions.
fig=figure('Color','w','Units','inches','Position',[1 1 11.7 5.8]);
ax=axes(fig);
hold(ax,'on');
masks={qSorted>=.05, qSorted<.05 & rSorted<0, qSorted<.05 & rSorted>0};
colors=[.76 .76 .76; .10 .43 .70; .82 .20 .16];
h=gobjects(3,1);
for k=1:3
    h(k)=scatter(ax,find(masks{k}),rSorted(masks{k}),14,colors(k,:),'filled');
end
yline(ax,0,'k-','HandleVisibility','off');
labelY=[.67 .53 .80 -.75];
markers={'d','s','p','v'};
for k=1:4
    row=candidateRows(k);
    px=position(row);
    scatter(ax,px,r(row),85,candidateColors(k,:),markers{k},'filled', ...
        'MarkerEdgeColor','k','LineWidth',.7,'HandleVisibility','off');
    if k<4
        tx=1144;
    else
        tx=145;
    end
    plot(ax,[px tx-10],[r(row) labelY(k)],'-','Color',candidateColors(k,:), ...
        'HandleVisibility','off');
    text(ax,tx,labelY(k),candidateLabels{k},'Interpreter','none', ...
        'FontSize',9,'Color',candidateColors(k,:),'VerticalAlignment','middle');
end
xlim(ax,[0 1550]);
ylim(ax,[-.82 .89]);
set(ax,'XTick',[]);
grid(ax,'on');
xlabel(ax,'1,106 ions ordered by Pearson correlation');
ylabel(ax,'Pearson r with normalized colonic IDO1');
title(ax,sprintf('36 mice; %d of 1,106 ions pass BH-FDR q < 0.05',sum(q<.05)), ...
    'FontWeight','normal');
legend(ax,h,{'q >= 0.05','q < 0.05, negative','q < 0.05, positive'}, ...
    'Location','northwest','Box','off','FontSize',9);
export_eps(fig,out,'fig4E_ion_waterfall',style.fontSize);

%% 9. Four per-mouse correlation plots; triangles retain historical R labels.
% Peak areas are normalized source values, divided by 1,000 for readability.
% IDO1 was normalized to the control mean within each experiment (15 or 16).
% Pearson r uses all 36 mice; q is from all 1,106 ions, not these four alone.
stems=["fig4F_gammaT3","fig4G_C29_3","fig4H_sterol","fig4I_proline"];
for k=1:4
    row=candidateRows(k);
    y=d.peakArea(:,row)/1000;
    fig=figure('Color','w','Units','inches','Position',[1 1 6.1 6.2]);
    ax=axes(fig);
    hold(ax,'on');
    h=gobjects(5,1);
    for g=1:4
        use=m.Group==style.groups(g) & ~m.Responder;
        h(g)=scatter(ax,x(use),y(use),style.pointSize,style.colors(g,:), ...
            'filled','MarkerEdgeColor','k','LineWidth',.45);
    end
    h(5)=scatter(ax,x(m.Responder),y(m.Responder),1.7*style.pointSize, ...
        style.colors(4,:),'^','filled','MarkerEdgeColor',style.responderEdge, ...
        'LineWidth',1.6);
    axis(ax,'square');
    grid(ax,'on');
    xlim(ax,[max(0,min(x)-.05*range(x)) max(x)+.08*range(x)]);
    ylim(ax,[max(0,min(y)-.08*range(y)) max(y)+.17*range(y)]);
    xlabel(ax,'Normalized colonic IDO1');
    ylabel(ax,'Normalized peak area (\times10^3)');
    title(ax,candidateLabels{k},'FontWeight','normal','Interpreter','none');
    text(ax,.03,.97,sprintf('Pearson r = %.3f; q = %.2g',r(row),q(row)), ...
        'Units','normalized','VerticalAlignment','top','BackgroundColor','w', ...
        'Margin',2,'FontSize',style.fontSize);
    legend(ax,h,[style.labels "AVN+P responder"],'Location','southoutside', ...
        'NumColumns',2,'Interpreter','none','Box','off','FontSize',9);
    export_eps(fig,out,stems(k),style.fontSize);
end
disp(result.ions(candidateRows,:));
fprintf('Figure 4/S4: seven EPS files saved in %s\n',out);

function q=bh_fdr(p)
% Benjamini-Hochberg: sort P values, adjust by number of ions, restore order.
[s,order]=sort(p);
n=numel(p);
adjusted=flipud(cummin(flipud(s.*n./(1:n)')));
q=zeros(size(p));
q(order)=min(1,adjusted);
end
