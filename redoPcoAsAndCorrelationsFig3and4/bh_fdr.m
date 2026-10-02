function q=bh_fdr(p)
% Benjamini-Hochberg: sort P values, adjust by number of ions, restore order.
[s,order]=sort(p);
n=numel(p);
adjusted=flipud(cummin(flipud(s.*n./(1:n)')));
q=zeros(size(p));
q(order)=min(1,adjusted);
end
