function atten = calculateSoilAttenuation(pathLength, freq, n, lossTangent)
% CALCULATESOILATTENUATION  Factor de atenuacion de campo (amplitud) por
% perdidas dielectricas en un medio no magnetico con indice de refraccion
% n y tangente de perdida lossTangent, a lo largo de una distancia
% geometrica pathLength recorrida dentro de dicho medio.
%
% Para un dielectrico con permitividad compleja
% eps_c = eps' * (1 - j*lossTangent), mu = mu0, la constante de
% propagacion es gamma = alpha + j*beta, con la constante de atenuacion
% (exacta, valida para cualquier lossTangent >= 0):
%
%   alpha = (omega/c) * sqrt( (n^2/2) * ( sqrt(1+lossTangent^2) - 1 ) )   [Np/m]
%
% que se anula cuando lossTangent = 0 (medio sin perdidas), reproduciendo
% el comportamiento previo del modelo.
%
% SYNOPSIS:
%   atten = calculateSoilAttenuation(pathLength, freq, n, lossTangent)
%
% INPUT:
%   pathLength  - Distancia geometrica recorrida dentro del suelo [m]
%                 (p.ej. r2Tx o r2Rx de calculateSlantRange), escalar o
%                 array.
%   freq        - Frecuencia portadora [Hz] (escalar).
%   n           - Indice de refraccion del suelo (escalar), n = sqrt(eps_r).
%   lossTangent - Tangente de perdida del suelo tan(delta) = eps''/eps'
%                 (escalar, >= 0).
%
% OUTPUT:
%   atten - Factor de atenuacion de campo (lineal, no en dB), mismo tamano
%           que pathLength. Se multiplica directamente sobre la amplitud
%           compleja de la senal.

c = 299792458;
omega = 2*pi*freq;

alpha = (omega/c) * sqrt( (n^2/2) * ( sqrt(1 + lossTangent^2) - 1 ) ); % Np/m

atten = exp(-alpha .* pathLength);

end
