% INICIALIZA  Rutas y carpetas de salida comunes a todos los experimentos.
RAIZ = fileparts(mfilename('fullpath'));
addpath(fullfile(RAIZ, 'lib'));
DIR_RES = fullfile(RAIZ, 'resultados');
DIR_FIG = fullfile(DIR_RES, 'figs');
if ~exist(DIR_FIG, 'dir'), mkdir(DIR_FIG); end
set(groot, 'defaultFigureVisible', 'off', 'defaultAxesFontSize', 9, ...
    'defaultLineLineWidth', 1.2);
