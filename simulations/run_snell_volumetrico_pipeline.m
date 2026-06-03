clc
clear
close all

projectDir = fileparts(mfilename('fullpath'));

% Shared parameter suffix for GS and PROC volumetrico: '_espiral', '_espiral_GRSS', '_espiral_EMIRADOS'
paramSuffix = '_espiral_EMIRADOS';

fprintf('Ejecutando generador de senales volumetrico (GS)...\n');
addpath(fullfile(projectDir,'gs'));
addpath(fullfile(projectDir,'proc'));

GS_snell_volumetrico_test_script(paramSuffix);

fprintf('Ejecutando procesador volumetrico (PROC)...\n');
PROC_snell_volumetrico_test_script(paramSuffix);

fprintf('Pipeline volumetrico GS->PROC finalizado correctamente.\n');
