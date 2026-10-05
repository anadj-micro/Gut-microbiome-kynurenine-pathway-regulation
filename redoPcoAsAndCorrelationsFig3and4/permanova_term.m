function out = permanova_term(D, Xreduced, Xterm, nPermutations, seed, strata)
%PERMANOVA_TERM Test a model term in a distance-based linear model.
% Term residuals after projection on the reduced covariate design are
% permuted jointly by row, optionally within strata. This preserves the
% fitted covariate component while testing the additional term.

arguments
    D double
    Xreduced double
    Xterm double
    nPermutations (1,1) double {mustBePositive, mustBeInteger} = 9999
    seed (1,1) double = 42
    strata = []
end

n = size(D, 1);
assert(size(D, 2) == n && max(abs(D - D'), [], 'all') < 1e-8, ...
    'D must be a symmetric square distance matrix.');
assert(size(Xreduced, 1) == n && size(Xterm, 1) == n, ...
    'Design matrices must have one row per sample.');

J = eye(n) - ones(n) / n;
G = -0.5 * J * (D .^ 2) * J;
[Fobserved, ssTerm, ssResidual, dfTerm, dfResidual, r2] = localF(G, Xreduced, Xterm);
X0 = localFullRank(Xreduced);
H0 = X0 * pinv(X0);
termFitted = H0 * Xterm;
termResidual = Xterm - termFitted;

rng(seed, 'twister');
Fperm = nan(nPermutations, 1);
for b = 1:nPermutations
    idx = (1:n)';
    if isempty(strata)
        idx = idx(randperm(n));
    else
        strataString = string(strata);
        for level = unique(strataString, 'stable')'
            members = find(strataString == level);
            idx(members) = members(randperm(numel(members)));
        end
    end
    permutedTerm = termFitted + termResidual(idx, :);
    Fperm(b) = localF(G, Xreduced, permutedTerm);
end

out = table(Fobserved, (1 + sum(Fperm >= Fobserved)) / (nPermutations + 1), ...
    r2, ssTerm, ssResidual, dfTerm, dfResidual, nPermutations, ...
    'VariableNames', {'PseudoF','PValue','PartialR2','SSTerm','SSResidual', ...
    'DFTerm','DFResidual','Permutations'});
end

function [F, ssTerm, ssResidual, dfTerm, dfResidual, r2] = localF(G, X0, Xt)
n = size(G, 1);
X0 = localFullRank(X0);
X1 = localFullRank([X0 Xt]);
H0 = X0 * pinv(X0);
H1 = X1 * pinv(X1);
ssTerm = trace((H1 - H0) * G);
ssResidual = trace((eye(n) - H1) * G);
dfTerm = rank(X1) - rank(X0);
dfResidual = n - rank(X1);
F = (ssTerm / dfTerm) / (ssResidual / dfResidual);
r2 = ssTerm / (ssTerm + ssResidual);
end

function X = localFullRank(X)
if isempty(X)
    X = ones(0, 1);
    return
end
[~, R, piv] = qr(X, 'econ', 'vector');
tolerance = max(size(R)) * eps(norm(R, inf));
keep = piv(1:sum(abs(diag(R)) > tolerance));
X = X(:, sort(keep));
end
