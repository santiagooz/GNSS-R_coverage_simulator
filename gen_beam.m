%% project: GNSSR_coverage_simulator
%
% Function that generates an antenna beam model using an elliptical cone.
% theta_x is the aperture along track and theta_y across track, both in
% degrees. L1 and L2 are the limits of the distance to the antenna
% dimension. L1=0 generates the beam starting from the LEO sat to L2.

function [X, Y, Z] = gen_beam(theta_x, theta_y, L1, L2)


r = linspace(L1,L2,round((L2-L1)/100));
th = linspace(0,2*pi) ;
[R,T] = meshgrid(r,th) ;

a_x = sqrt(1/cos(theta_x*pi/180)^2-1);
a_y = sqrt(1/cos(theta_y*pi/180)^2-1);

X = R.*cos(T)*a_x;
Y = R.*sin(T)*a_y;
Z = -R;
end

