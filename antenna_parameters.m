%% project: GNSSR_coverage_simulator
%
% Function that extracts relevant parameters from the antenna pattern
% described in G0 with an angle resolution dt.
%
% el: elevation dimension [degs] - az: azimuth dimension [degs]
% theta_x: aperture angle along track [degrees]
% theta_y: aperture angle across track [degrees]


function [el, az, theta_x, theta_y] = antenna_parameters(G0, dt)

[a, b] = size(G0);
el = (0 : a-1)*dt;
az = (0 : b-1)*dt;

x0 = (az == 0); x1 = (az == 180);
y0 = (az == 90); y1 = (az == 270);

G_max = G0(1,1);
el_aux = el(1) : 0.1 : el(end);
G_aux = interp1(el,G0(:, x0),el_aux,'cubic');
[~, tx0] = min(abs(G_aux - G_max + 3));
G_aux = interp1(el,G0(:, x1),el_aux,'cubic');
[~, tx1] = min(abs(G_aux - G_max + 3));
G_aux = interp1(el,G0(:, y0),el_aux,'cubic');
[~, ty0] = min(abs(G_aux - G_max + 3));
G_aux = interp1(el,G0(:, y1),el_aux,'cubic');
[~, ty1] = min(abs(G_aux - G_max + 3));
theta_x = (el_aux(tx0)+el_aux(tx1))/2;
theta_y = (el_aux(ty0)+el_aux(ty1))/2;
end

