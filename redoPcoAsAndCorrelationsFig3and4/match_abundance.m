function [x,names]=match_abundance(t,ids)
% Match rows by ID. Blank abundance cells mean zero in these source exports.
assert(numel(unique(t{:,1}))==height(t),'Duplicate abundance sample IDs.');
[found,rows]=ismember(ids,t{:,1});
assert(all(found),'A mouse lacks abundance data.');
x=t{rows,2:end};
x(isnan(x))=0;
assert(all(isfinite(x) & x>=0,'all'),'Invalid abundance.');
keep=any(x>0,1);
x=x(:,keep); % all-zero taxa contribute nothing to distances
names=string(t.Properties.VariableNames(2:end));
names=names(keep);
end
