%% leer_rcs.m  -  Lee e interpreta un resultado de simular.py (.mat)
%
% Contenido del .mat:
%   az_deg (1 x Naz), el_deg (1 x Nel)      ángulos [grados]
%   rcs_m2, rcs_dbsm (Nel x Naz)            RCS co-polar (filas = elevación, columnas = azimut)
%   cross_m2, cross_dbsm (Nel x Naz)        polarización cruzada (vacía si el método es PO)
%   freq_Hz, lambda_m, pol, metodo, forma, dim_m
%   rcs_max_dbsm, az_max_deg, el_max_deg
%   mesh_vertices_m (N x 3), mesh_faces (M x 3, base 1)
%   descripcion                             convenciones en texto
%
% Convención: azimut en el plano xy desde +x; elevación desde el plano xy (+90° = +z).

clear; close all;
archivo = fullfile('resultados', 'cubo_1m_400MHz_TMTM_mom.mat');   % <-- cambia aquí
r = load(archivo);

%% 1. Resumen
fprintf('%s  %s m | f = %g MHz | lambda = %.3f m | %s | %s\n', r.forma, ...
        mat2str(r.dim_m), r.freq_Hz/1e6, r.lambda_m, r.metodo, r.pol);
fprintf('RCS máxima: %.2f dBsm en az = %g°, el = %g°\n', r.rcs_max_dbsm, r.az_max_deg, r.el_max_deg);

% RCS media sobre la esfera (ponderada con cos(el), en m^2)
w = cosd(r.el_deg(:));
sigma_media = sum(mean(r.rcs_m2, 2) .* w) / sum(w);
fprintf('RCS media (esfera): %.2f dBsm\n', 10*log10(sigma_media));

%% 2. Consultar la RCS en un ángulo cualquiera (interpolación)
az_q = 30; el_q = 20;
s_q = interp2(r.az_deg, r.el_deg, r.rcs_dbsm, az_q, el_q);
fprintf('RCS en az = %g°, el = %g°: %.2f dBsm\n', az_q, el_q, s_q);

%% 3. Figuras
vmax = max(r.rcs_dbsm(:));  vmin = vmax - 40;
figure('Name', 'RCS', 'Position', [80 80 1300 800]);

% 3a. Mapa azimut x elevación
subplot(2,2,1);
imagesc(r.az_deg, r.el_deg, r.rcs_dbsm, [vmin vmax]); axis xy;
colormap(gca, turbo); cb = colorbar; cb.Label.String = 'RCS [dBsm]';
xlabel('azimut [°]'); ylabel('elevación [°]');
title(sprintf('%s  %s @ %g MHz', r.forma, r.pol, r.freq_Hz/1e6));

% 3b. Corte en azimut (el = 0°)
subplot(2,2,2);
[~, i0] = min(abs(r.el_deg));
plot(r.az_deg, r.rcs_dbsm(i0,:), 'LineWidth', 1.5); hold on;
if ~isempty(r.cross_dbsm), plot(r.az_deg, r.cross_dbsm(i0,:), 'Color', [.6 .6 .6]); end
grid on; ylim([vmin vmax+3]); xlim([0 360]);
xlabel('azimut [°]'); ylabel('RCS [dBsm]'); title('Corte en azimut, el = 0°');
legend([{r.pol}, repmat({r.pol_cross}, 1, ~isempty(r.cross_dbsm))], 'Location', 'southoutside', 'Orientation', 'horizontal');

% 3c. Malla del objeto
subplot(2,2,3);
patch('Vertices', r.mesh_vertices_m, 'Faces', r.mesh_faces, ...
      'FaceColor', [.8 .83 .88], 'EdgeColor', [.2 .2 .3], 'LineWidth', 0.3);
axis equal; view(35, 25); camlight; lighting flat; grid on;
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title(sprintf('Malla: %d triángulos', size(r.mesh_faces, 1)));

% 3d. Diagrama 3D (radio = dBsm por encima de vmin)
subplot(2,2,4);
[AZ, EL] = meshgrid(r.az_deg, r.el_deg);
R = max(r.rcs_dbsm - vmin, 0);
surf(R.*cosd(EL).*cosd(AZ), R.*cosd(EL).*sind(AZ), R.*sind(EL), r.rcs_dbsm, 'EdgeColor', 'none');
axis equal; colormap(gca, turbo); caxis([vmin vmax]); view(35, 25); camlight;
title('Diagrama 3D de RCS [dBsm]');
