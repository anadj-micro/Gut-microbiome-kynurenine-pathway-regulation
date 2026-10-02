function [xy,percent] = draw_pcoa(abundance,mice,style,folder,names)
% Compute one ASV Bray-Curtis PCoA and export its two color views.
% Input rows are already matched to mice. No prevalence filtering or CLR.
x=full(abundance);
d=squareform(pdist(x,@(a,b)sum(abs(b-a),2)./sum(b+a,2)));
[xy,eigenvalues]=cmdscale(d);
percent=100*eigenvalues(1:2)/sum(eigenvalues(eigenvalues>0));
xy=xy(:,1:2);
% Orient the control centroid consistently without loading saved coordinates.
% Axis reflection has no biological meaning; it only sets the display direction.
center=mean(xy(mice.Group=="C",:),1);
xy=xy.*(style.controlCorner./sign(center));
responder=false(height(mice),1);
if ismember('Responder',mice.Properties.VariableNames)
    responder=mice.Responder;
end

% Each EPS is a separate, square panel with its own legend or color scale.
for view=1:2
    fig=figure('Color',style.backgroundColor, ...
        'Units','inches','Position',[1 1 7.8 5.8]);
    ax=axes(fig);
    ax.Color=style.backgroundColor;
    hold(ax,'on');
    if view==1
        h=gobjects(numel(style.groups),1);
        for g=1:numel(style.groups)
            use=mice.Group==style.groups(g) & ~responder;
            if style.borderSize(g)>0
                % Bordered groups: AVN+FMT and Recovery.
                % borderSize(g) is now a LineWidth (points), not an area delta.
                h(g)=scatter(ax, ...
                    xy(use,1),xy(use,2), ...
                    style.pointSize, ...
                    style.colors(g,:), ...
                    'filled', ...
                    'Marker','o', ...
                    'MarkerEdgeColor',style.borderColor, ...
                    'LineWidth',style.borderSize(g));
            else
                % Normal groups: no border.
                h(g)=scatter(ax, ...
                    xy(use,1),xy(use,2), ...
                    style.pointSize, ...
                    style.colors(g,:), ...
                    'filled', ...
                    'Marker','o', ...
                    'MarkerEdgeColor','none');
            end
        end
        labels=style.labels;
        if any(responder)
            g=find(style.groups=="AVN_P");
            h(end+1)=scatter(ax,xy(responder,1),xy(responder,2), ...
                1.7*style.pointSize,style.colors(g,:),'^','filled', ...
                'MarkerEdgeColor',style.responderEdge,'LineWidth',1.6);
            labels(end+1)="AVN+P responder";
        end
        legend(ax,h,labels,'Location','eastoutside','Box','off', ...
            'Interpreter','none','FontSize',style.fontSize-1);
        title(ax,'ASV Bray-Curtis PCoA by mouse group','FontWeight','bold');
    else
        scatter(ax,xy(~responder,1),xy(~responder,2),style.pointSize, ...
            mice.IDO1(~responder),'filled','MarkerEdgeColor','k','LineWidth',0.4);
        if any(responder)
            scatter(ax,xy(responder,1),xy(responder,2),1.7*style.pointSize, ...
                mice.IDO1(responder),'^','filled', ...
                'MarkerEdgeColor',style.responderEdge,'LineWidth',1.6);
        end
        colormap(ax,style.colorMap);
        c=colorbar(ax);
        c.Label.String='Colonic IDO1';
        title(ax,'Same PCoA colored by colonic IDO1','FontWeight','bold');
    end
    xlabel(ax,sprintf('PCoA1 (%.1f%%)',percent(1)),'FontWeight','bold');
    ylabel(ax,sprintf('PCoA2 (%.1f%%)',percent(2)),'FontWeight','bold');
    axis(ax,'square');
    grid(ax,'on');
    box(ax,'on');
    % Use identical limits for both views, independent of legend placement.
    span=max(xy)-min(xy);
    xlim(ax,[min(xy(:,1)) max(xy(:,1))]+[-1 1]*0.12*span(1));
    ylim(ax,[min(xy(:,2)) max(xy(:,2))]+[-1 1]*0.12*span(2));
    export_eps(fig,folder,names(view),style.fontSize);
end
end
