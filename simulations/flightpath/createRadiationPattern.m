function padrao = createRadiationPattern(abertura_el, abertura_az)
% CREATERADIATIONPATTERN Crea un patrón de radiación con aberturas específicas
%
% Sintaxis:
%   padrao = createRadiationPattern(abertura_el, abertura_az)
%
% Entradas:
%   abertura_el - Abertura en elevación (grados) [por defecto: 70°]
%   abertura_az - Abertura en azimut (grados) [por defecto: 20°]
%
% Salidas:
%   padrao - Matriz 181x181 con patrón de radiación normalizado
%            Filas: elevación de -90° a +90°
%            Columnas: azimut de -90° a +90°
%
% Ejemplo:
%   padrao = createRadiationPattern(70, 20);
%   imagesc(padrao); colorbar; title('Patrón de Radiación');

% Valores por defecto
if nargin < 1
    abertura_el = 70;  % Abertura en elevación (grados)
end
if nargin < 2
    abertura_az = 20;  % Abertura en azimut (grados)
end

% Crear vectores de ángulos de -90° a +90°
angulos = -90:1:90;  % 181 puntos
n_puntos = length(angulos);

% Crear matrices de coordenadas
[az_grid, el_grid] = meshgrid(angulos, angulos);

% Calcular el ancho del haz a -3dB (half-power beamwidth)
% Convertir abertura total a ancho de haz a -3dB
ancho_el_3db = abertura_el / 2;  % ±35° desde el centro
ancho_az_3db = abertura_az / 2;  % ±10° desde el centro

% Crear patrón gaussiano para cada dirección
% Parámetro de forma para gaussiana (ajustado para -3dB en los límites)
sigma_el = ancho_el_3db / sqrt(2 * log(2));  % Para -3dB en los límites
sigma_az = ancho_az_3db / sqrt(2 * log(2));

% Patrón de radiación gaussiano separable
patron_el = exp(-(el_grid.^2) / (2 * sigma_el^2));
patron_az = exp(-(az_grid.^2) / (2 * sigma_az^2));

% Combinar patrones (producto de funciones separables)
padrao = patron_el .* patron_az;

% Normalizar a 1 (0 dB en el máximo)
padrao = padrao / max(padrao(:));

end