function d = bray_curtis_distance(x, y)
%BRAY_CURTIS_DISTANCE Bray-Curtis dissimilarity for pdist.

denominator = sum(abs(x) + abs(y), 2);
numerator = sum(abs(x - y), 2);
d = numerator ./ denominator;
d(denominator == 0) = 0;
end
