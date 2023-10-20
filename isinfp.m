%% project: GNSSR_coverage_simulator
%
% Function that decides if the specular reflection point SP is located
% the footprint of the antenna beam descibed by X, Y and Z.
% described in G0 with an angle resolution dt.
%
% Returns true if SP is within the antenna footprint and flase otherwise.

function [isin] = isinfp(sp, X, Y, Z)

[~, aa] = min(abs(sp(3) - Z(1,:)));
bb = abs(X(:,aa) - sp(1));
[~, bb] = sort(bb);
yy = sort(Y(bb(1:2),aa));
isin = yy(1)<sp(2) && yy(2)>sp(2) && ~all(X(:,aa)<sp(1)) && ~all(X(:,aa)>sp(1));

end

