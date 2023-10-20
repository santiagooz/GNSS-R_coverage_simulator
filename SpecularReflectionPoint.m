function [refl_ECEF] = SpecularReflectionPoint(rGPS,rLEO)
% Lucas Sanz, 2021.
% UIDET-SENyT, Facultad de Ingenieria, UNLP.
%
% Funcion para calcular el punto de reflexion especular dada la posicion
% de los satélites GPS y LEO en un instante.
%
% ENTRADAS:
%   - rGPS (1x3): en coordenadas ECEF y en m.
%   - rLEO (1x3): en coordenadas ECEF y en m.
% SALIDA:
%   - llaRefl (1x3). Punto de reflexion en formato LLA

% Paso 1) Transformo las posiciones rGPS y rLEO hacia una esfera unitaria
E = referenceEllipsoid(7030);
a = E.SemimajorAxis;
b = E.SemiminorAxis;
Tmatrix = [1/a 0 0; 0 1/a 0; 0 0 1/b];

v = Tmatrix * rGPS(:);
w = Tmatrix * rLEO(:);

% Paso 2) Ortogonalizacion de Gram-Schmidt. Base que genera el plano XY que
% contiene a Tx, Rx, y el origen.
u1 = v;
u2 = w - (sum(w.*u1)/sum(u1.*u1))*u1;
% Normalizo
e1 = u1/norm(u1);
e2 = u2/norm(u2);

% Paso 3) Busco las coordenadas x_plano, y_plano (sabiendo que z_plano = 0)
% (el plano es el que contiene a Tx, Rx y Punto de Refl)
xt = dot(v,e1); % Coordenada x en la nueva base, del transmisor
yt = dot(v,e2); % Coordenada y en la nueva base, del transmisor

xr = dot(w,e1); % Coordenada x en la nueva base, del receptor
yr = dot(w,e2); % Coordenada y en la nueva base, del receptor

% Paso 4) Busqueda del punto de reflexion con la "ecuacion de espejo esférico"

RE = 1; % Se deja en 1 porque varie el método original hacia uno
% que escala las coordenadas a una esfera unitaria.

c0 = xt*yr + yt*xr - RE * (yt+yr);
c1 = -4*(xt*xr - yt*yr) + 2*RE*(xt+xr);
c2 = -6*(xt*yr + yt*xr);
c3 = 4*(xt*xr - yt*yr) + 2*RE*(xt+xr);
c4 = xt*yr + yt*xr + RE * (yt+yr);

polcoef = [c4 c3 c2 c1 c0];

% Por la base de coordenadas elegida, el "eje x" apunta hacia el satélite
% GPS y el angulo phi donde esta el punto de reflexion (medido desde x
% hacia y antihorario) nunca estara fuera del intervalo -90 a +90 (a
% checkear)

% Por la base coordenada elegida el satélite GPS esta a 0° del eje x.
% El satélite LEO está a phiRx del eje x hacia el eje y antihorario. Por
% lo tanto el punto de reflexión debe estar a un angulo phi entre 0° y
% PhiRx

PhiRx = atan2(yr,xr);
Phi = 2*atan(roots(polcoef));
cnt = 0;
Phi_out = nan;
refl_ECEF = nan;

if (PhiRx < 0)
    for i = 1:4
        if (Phi(i) < 0 && Phi(i) > PhiRx)
            Phi_out = Phi(i);
            cnt = cnt + 1;
        end
    end
end

if (PhiRx > 0)
    for i = 1:4
        if (Phi(i) > 0 && Phi(i) < PhiRx)
            Phi_out = Phi(i);
            cnt = cnt + 1;
        end
    end
end

% if (cnt ~= 1)
%     disp("Advertencia: se encontro más de una posible solución\n");
% end

if ~isnan(Phi_out)
    % Ahora paso de Phi_out a coordenadas x,y del punto de reflexión sobre el
    % plano definido
    
    x_refl = cos(Phi_out);
    y_refl = sin(Phi_out);
    
    % Traslado de nuevo a R3 y en el sistema escalado a una esfera unitaria.
    
    refl_Esf = x_refl*e1 + y_refl*e2;
    
    % Paso 5) Vuelvo al sistema ECEF escalado a la Tierra Eliptica
    Tinv = [a 0 0; 0 a 0; 0 0 b];
    refl_ECEF = Tinv * refl_Esf;
end

end

