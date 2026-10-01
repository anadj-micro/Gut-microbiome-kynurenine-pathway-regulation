function export_eps(fig,folder,name,fontSize)
% Save editable vector artwork. Re-running replaces this EPS, not the data.
if ~isfolder(folder)
    mkdir(folder);
end
set(findall(fig,'Type','axes'),'FontName','Arial','FontSize',fontSize, ...
    'Box','off','LineWidth',0.8,'TickDir','out','GridAlpha',0.12);
set(fig,'Color','w','Renderer','painters','PaperPositionMode','auto');
drawnow;
print(fig,fullfile(folder,[char(name) '.eps']),'-depsc','-painters');
end
