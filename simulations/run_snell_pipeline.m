clc
clear
close all

projectDir = fileparts(mfilename('fullpath'));

% Shared parameter suffix for GS and PROC: '_espiral', '_espiral_GRSS', '_espiral_EMIRADOS'
paramSuffix = '_espiral';

fprintf('Ejecutando generador de senales (GS)...\n');
addpath(fullfile(projectDir,'gs'));
addpath(fullfile(projectDir,'proc'));

GS_snell_test_script(paramSuffix);

fprintf('Ejecutando procesador (PROC)...\n');
PROC_snell_test_script(paramSuffix);

fprintf('Pipeline GS->PROC finalizado correctamente.\n');
