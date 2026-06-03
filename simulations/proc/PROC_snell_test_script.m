function PROC_snell_test_script(paramSuffix)
clc
close all

scriptDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(scriptDir);

addpath(genpath(fullfile(scriptDir,'src')))
addpath(genpath(fullfile(projectDir,'tools')))
addpath(genpath(fullfile(projectDir,'common')))

% Select parameter set suffix: '_espiral', '_espiral_GRSS', '_espiral_EMIRADOS'
if nargin < 1 || strlength(string(paramSuffix)) == 0
    paramSuffix = '_espiral';
else
    paramSuffix = char(string(paramSuffix));
end

figOutputDir = fullfile(projectDir,'io','snell','figure');
if ~exist(figOutputDir,'dir')
    mkdir(figOutputDir);
end

% Parametros del transmisor

radarJSON  = json2struct(fullfile(projectDir,'parametros',sprintf('radarTx%s.json',paramSuffix)));
strRadarTx = radarJSON.radar; clear radarJSON;

% Parametros del receptor

radarJSON  = json2struct(fullfile(projectDir,'parametros',sprintf('radarRx%s.json',paramSuffix)));
strRadarRx = radarJSON.radar; clear radarJSON;

% Parametros do sistema

systemJSON = json2struct(fullfile(projectDir,'parametros',sprintf('system%s.json',paramSuffix)));
strSystem  = systemJSON.system; clear systemJSON;

% Target 

targetJSON = json2struct(fullfile(projectDir,'parametros',sprintf('target%s.json',paramSuffix)));
strTarget  = targetJSON.target; clear targetJSON;

x0  = strTarget.pos(1);
y0  = strTarget.pos(2);

% xy0 = strTarget.pos(1); dxy = 0.04;  lxy = 0.12;
z0  = strTarget.pos(3); %  dz = 0.08;   lz = 0.24;

% Calculate params

strRadarTx.lamb = strSystem.VelocidadeLuz/strRadarTx.FreqPortadora; % Comprimento de onda


%% Load 

rawDataFile = fullfile(projectDir,'io','snell','raw',sprintf('bistatic_sar_snell_raw%s.mat', paramSuffix));
if ~exist(rawDataFile,'file')
    error('No se encontro el archivo de entrada: %s. Ejecuta primero GS_snell_test_script.m', rawDataFile);
end
load(rawDataFile)

%% Preproc 

[Nrng,Nazm] = size(rootData);
ratio = 1;

rcompData = interpft(rootData,Nrng*ratio,1);
rcompData = rcompData.';
fs = ratio*strRadarTx.fs;

clear rootData

%% Processamiento 

tic

wb = waitbar(0);
threshold = 1e-11;
c = strSystem.VelocidadeLuz;
outputData = zeros(size(X));

n2 = n;

for m = 1:numel(X)

    % Show progress
%     if coder.target('MATLAB')
%         fprintf('--- PROCESANDO %d de %d ---\n', m, numel(X));
%     end
    
    % Slant range calculation
    [R1t,R2t] = slantRange(PxT,PyT,PzT,X(m),Y(m),Z(m),n2,threshold);
    [R1r,R2r] = slantRange(PxR,PyR,PzR,X(m),Y(m),Z(m),n2,threshold);
    
    % Fractional range bin sample
    t = (1/c)*(R1t + n2*R2t + n2*R2r + R1r);
    rngBin = 1 + t*fs;
    
    % Phase compensation term
    phi = (2*pi/strRadarTx.lamb)*(R1t + n2*R2t + n2*R2r + R1r);
    
    % Data accumulation
    allPower = 0;

    for idxRng = 1:length(rngBin)
        idxNext = ceil(rngBin(idxRng));
        idxLast = floor(rngBin(idxRng));

        q = rngBin(idxRng) - idxLast;

        % Calcular datos interpolados y contribución de potencia
        dataLast = rcompData(idxRng,idxLast);
        dataNext = rcompData(idxRng,idxNext);
        interpData = (1-q)*dataLast + q*dataNext;
        powerContrib = interpData * exp(1i*phi(idxRng));
        
        allPower = allPower + powerContrib;
    end

    outputData(m) = allPower;
    
    waitbar(m/numel(X),wb,'Ejecutando algoritmo Back-Projection...');
end
close(wb)
toc

%% Save processed output
processedDir = fullfile(projectDir,'io','snell','processed');
if ~exist(processedDir,'dir')
    mkdir(processedDir);
end

processedFile = fullfile(processedDir,sprintf('bistatic_sar_snell_bp_output%s.mat', paramSuffix));
save(processedFile,'outputData','X','Y','Z','x0','y0','z0');
fprintf('Salida de backprojection guardada en: %s\n', processedFile);

%% Plots 

figSlice = figure;

s = slice(X,Y,Z,abs(outputData/max(outputData(:))),x0,y0,-5);
set(s,'EdgeColor','none')
% title({'Output Data','abs()'})
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
axis equal
colormap parula
saveas(figSlice, fullfile(figOutputDir,'slice_abs_output_snell.png'));

% -------------------------------------------------------------------

figIso = figure;
dataDB = 20*log10(abs(outputData/max(outputData(:))));
% patch(isosurface(X,Y,Z,dataDB,-13),...
%     'EdgeColor','none','CData',-13,'FaceColor','flat','FaceAlpha',0.15);
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
% contourslice(X,Y,Z,dataDB,[],[],-1.1)
patch(isosurface(X,Y,Z,dataDB,-3),...
    'EdgeColor','none','CData',-3,'FaceColor','flat');
axis equal
% xlim([-0.15 0.15])
% ylim([-0.15 0.15])
% zlim([-2 -0.5])

xlim([min(X,[],'all') max(X,[],'all')])
ylim([min(Y,[],'all') max(Y,[],'all')])
zlim([min(Z,[],'all') max(Z,[],'all')])

grid on
grid minor
view([-45 45])
camlight right
% camlight left
camlight headlight
lighting gouraud
colormap jet
caxis([-40 0])
saveas(figIso, fullfile(figOutputDir,'isosurface_output_snell.png'));

disp('Resolucion automatica a -3 dB:');
resolutionCsvFile = fullfile(processedDir, sprintf('bistatic_sar_snell_resolution_3dB%s.csv', paramSuffix));
report3dBResolution(X, Y, Z, outputData, resolutionCsvFile);

end