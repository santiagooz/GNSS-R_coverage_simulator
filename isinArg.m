%% project: GNSSR_coverage_simulator
%
% Function that checks if the point P is over Argentina. To represent the
% territory we use a rectangular patch defined by the limit values of
% latitude and longitude of Argentina.

function [stat] = isinArg(P)

lla = ecef2lla(P, 'WGS84');

PS = polyshape([-21.783297 -21.783297 -55.116443 -55.116443], [-73.554345 -53.637448 -53.637448 -73.554345]);
stat = isinterior(PS, lla(:, 1), lla(:, 2));

end

