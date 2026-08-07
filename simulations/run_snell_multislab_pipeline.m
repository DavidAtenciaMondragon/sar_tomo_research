clc
clear
close all

projectDir = fileparts(mfilename('fullpath'));

% Shared parameter suffix for GS and PROC: '_multislab'
paramSuffix = '_multislab';

fprintf('Ejecutando generador de senales (GS)...\n');
addpath(fullfile(projectDir,'gs'));
addpath(fullfile(projectDir,'proc'));

GS_snell_multislab_test_script(paramSuffix);

fprintf('Ejecutando procesador (PROC)...\n');
PROC_snell_multislab_test_script(paramSuffix);

fprintf('Pipeline GS->PROC (multicapas) finalizado correctamente.\n');