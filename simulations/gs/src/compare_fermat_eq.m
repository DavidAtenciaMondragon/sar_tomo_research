clc
clear
close all

addpath(genpath('src'))
addpath(genpath(strcat('..',filesep,'tools')))
addpath(genpath(strcat('..',filesep,'common')))

%% 1. Configuración del Entorno y Carga del DEM
% --- PARÁMETROS FÍSICOS ---
n1 = 1.0;      % Índice de refracción del medio 1 (aire)
n2 = 1.5;      % Índice del medio 2 (vidrio/agua/otro)
c  = 3e8;      % Velocidad de la luz en el vacío (m/s)

% --- CARGA DEL MODELO DE ELEVACIÓN DIGITAL (DEM) ---

fprintf('Creando una superficie de ejemplo.\n');
[X_DEM, Y_DEM] = meshgrid(-150:30:150, -150:30:150);
Z_DEM = 50 * sin(X_DEM/100) + 30 * cos(Y_DEM/100);


% --- CREACIÓN DE LA MALLA DE COORDENADAS (Grid) ---
resolution = 30; % Discretización de tu DEM en metros
[rows, cols] = size(Z_DEM);
x_origen = 0; 
y_origen = 0;
X_vec = x_origen + (0:cols-1) * resolution;
Y_vec = y_origen + (0:rows-1) * resolution;

strDEM.X_vec = X_vec;
strDEM.Y_vec = Y_vec;
strDEM.Z_DEM = Z_DEM;

%% 2. Definición de Puntos de Interés
% Tx = [ 100,  50, 120;];  % Punto A (Transmisor)
   Tx = [ 100,  50, 120;
       120,  30, 130;];  % Punto A (Transmisor)
% Rx = [ 200, 200, 110];  % Punto B (Receptor)
   Rx = [ 200, 200, 110;
       280, 30, 100;];  % Punto B (Receptor)
P  = [ 150, 140, -40];  % Punto C (Debajo de la superficie)

bPlotVerbose = false;

[strReflexao, strRefraccoes]   = calculaSlantRangeFermat(strDEM,Tx,Rx,P,n1,n2,bPlotVerbose);
[strReflexao2, strRefraccoes2] = calculaSlantRangeFermat_eff(strDEM,Tx,Rx,P,n1,n2,5);

% Mostrar errores 

disp(strcat("Error en el punto de reflexion: ",num2str(mean(strReflexao.P_reflex - strReflexao2.P_reflex,2))));
disp(strcat("Error en la refraccion (ida): ",num2str(mean(strRefraccoes.P_refrac_ida - strRefraccoes2.P_refrac_ida,2))));
disp(strcat("Error en la refraccion (vuelta): ",num2str(mean(strRefraccoes.P_refrac_volta - strRefraccoes2.P_refrac_volta,2))));
