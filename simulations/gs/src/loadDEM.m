function [Xg, Yg, Zg] = loadDEM(name)

[Zg, R] = readgeoraster(name);  % Z: matriz de datos, R: referencia geográfica

% Obtener el número de filas y columnas
[nRows, nCols] = size(Zg);

% Crear vectores de latitudes y longitudes
lat = linspace(R.LatitudeLimits(1), R.LatitudeLimits(2), nRows);
lon = linspace(R.LongitudeLimits(1), R.LongitudeLimits(2), nCols);

% Crear el meshgrid
[Lon, Lat] = meshgrid(lon, lat);

% Determinar la zona UTM automáticamente
meanLat = mean(lat(:));
meanLon = mean(lon(:));
zone = utmzone(meanLat, meanLon);

% Obtener el sistema de coordenadas UTM
utmstruct = defaultm('utm');
utmstruct.zone = zone;
utmstruct.geoid = wgs84Ellipsoid;
utmstruct = defaultm(utmstruct);

% Convertir lat/lon a coordenadas UTM (metros)
[Xg, Yg] = mfwdtran(utmstruct, Lat, Lon);

% Centrando en el origen de coordenadas
X_mean = (max(Xg,[],'all') + min(Xg,[],'all'))/2;
Y_mean = (max(Yg,[],'all') + min(Yg,[],'all'))/2;

Xg = Xg - X_mean;
Yg = Yg - Y_mean;

end