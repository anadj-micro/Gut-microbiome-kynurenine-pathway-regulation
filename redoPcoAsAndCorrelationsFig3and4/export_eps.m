function export_eps(fig,folder,name,fontSize)
% Save editable vector artwork. Re-running replaces this EPS, not the data.
if ~isfolder(folder)
    mkdir(folder);
end
ax_all=findall(fig,'Type','axes');
set(ax_all,'FontName','Arial','FontSize',fontSize, ...
    'Box','off','LineWidth',0.8,'TickDir','out','GridAlpha',0.12);
for k=1:numel(ax_all)   % <-- CHANGED (new loop: manual 4-sided border, no mirrored ticks)
    a=ax_all(k);
    rectangle(a,'Position',[a.XLim(1) a.YLim(1) diff(a.XLim) diff(a.YLim)], ...
        'EdgeColor','k','LineWidth',0.8);
end
drawnow;
exportgraphics(fig,fullfile(folder,[char(name) '.eps']),'ContentType','vector');   % <-- CHANGED (replaced print/PaperPosition entirely)
end
