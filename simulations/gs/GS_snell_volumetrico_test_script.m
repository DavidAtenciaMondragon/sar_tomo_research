function GS_snell_volumetrico_test_script(paramSuffix)
clc
close all

scriptDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(scriptDir);

addpath(genpath(fullfile(scriptDir,'src')))
addpath(genpath(fullfile(projectDir,'tools')))
addpath(genpath(fullfile(projectDir,'common')))

% Select parameter set suffix: '_espiral', '_espiral_GRSS', '_espiral_EMIRADOS'
if nargin < 1 || strlength(string(paramSuffix)) == 0
    paramSuffix = '_espiral_EMIRADOS';
else
    paramSuffix = char(string(paramSuffix));
end

figOutputDir = fullfile(projectDir, 'io', 'snell_volumetrico', 'figures');
if ~exist(figOutputDir, 'dir')
    mkdir(figOutputDir);
end

rawOutputDir = fullfile(projectDir, 'io', 'snell_volumetrico', 'raw');
if ~exist(rawOutputDir, 'dir')
    mkdir(rawOutputDir);
end

% Parametros del transmisor

radarJSON  = json2struct(fullfile(projectDir,'parametros',sprintf('radarTx%s.json',paramSuffix)));
strRadarTx = radarJSON.radar; clear radarJSON;

% Parametros del receptor

radarJSON  = json2struct(fullfile(projectDir,'parametros',sprintf('radarRx%s.json',paramSuffix)));
strRadarRx = radarJSON.radar; clear radarJSON;

% Parametros del target

targetJSON = json2struct(fullfile(projectDir,'parametros',sprintf('target%s.json',paramSuffix)));
strTarget  = targetJSON.target; clear targetJSON;

% Parametros do sistema

systemJSON = json2struct(fullfile(projectDir,'parametros',sprintf('system%s.json',paramSuffix)));
strSystem  = systemJSON.system; clear systemJSON;

% Calculate params

strRadarTx.lamb = strSystem.VelocidadeLuz/strRadarTx.FreqPortadora; % Comprimento de onda

%% Create trajectory

% Transmissor
[PxT, PyT, PzT] = funcao_espiral(strRadarTx.NumVoltasEsp, strRadarTx.RaioMenorEsp, strRadarTx.RaioMaiorEsp,...
    strRadarTx.AltMaiorEsp,strRadarTx.AltMenorEsp, strRadarTx.Vt, strRadarTx.PRF,0);

% Receptor
[PxR, PyR, PzR] = funcao_espiral(strRadarRx.NumVoltasEsp, strRadarRx.RaioMenorEsp, strRadarRx.RaioMaiorEsp,...
    strRadarRx.AltMaiorEsp,strRadarRx.AltMenorEsp, strRadarRx.Vt, strRadarRx.PRF,0);

% Plot Trajectory

[Xg,Yg] = meshgrid([-strRadarTx.RaioMaiorEsp strRadarTx.RaioMaiorEsp]);
Zg = zeros(size(Xg));

figTrajectory = figure;
hold on
plot3(PxT, PyT, PzT,'k')
plot3(PxR, PyR, PzR,'r')
% plot3(strTarget.pos(1),strTarget.pos(2),strTarget.pos(3),'o','MarkerFaceColor','k');
surf(Xg,Yg,Zg,'EdgeColor','none','FaceColor',[1 0 1],'FaceAlpha',0.1)

if strSystem.EscenarioExterno
    load('utils\linePoints_targetMask.mat');
    [X,Y,Z] = meshgrid(xAxis,yAxis,zAxis);
    
else
    % Create grid
    [X,Y,Z,xAxis,yAxis,zAxis] = createProcessingGridFromSystem(strSystem, strTarget, 'gridSnellVolumetrico');
    
    % Create target points
    p1 = strTarget.lineSegment.p1;
    p2 = strTarget.lineSegment.p2;
    dxy = strSystem.gridSnellVolumetrico.dxy;
    
    fprintf('Creating 3D mask for line points...\n');
    [linePoints, targetMask] = createLinePointsAndMask(p1, p2, dxy, xAxis, yAxis, zAxis);
    fprintf('Mask created with %d non-zero points\n', sum(targetMask(:)));
end

plot3(linePoints(:,1),linePoints(:,2),linePoints(:,3),'.');

xlabel('X')
ylabel('Y')
zlabel('Z')
view([-35 28]);
axis vis3d
camproj perspective
grid minor
hold off
legend("Tx","Rx","Target")
saveas(figTrajectory, fullfile(figOutputDir, 'trayectorias.png'));

%% Create mask of line points in the grid
if ~strSystem.EscenarioExterno
    visualizeLineMask3D(targetMask, linePoints, xAxis, yAxis, zAxis, fullfile(figOutputDir, 'mascara_objetivo_3d.png'));
else
    figure,plot3(linePoints(:,1),linePoints(:,2),linePoints(:,3),'kx');
end
%% Create raw data

if ~isfield(strSystem, 'IndiceRefracaoSolo')
    error('Falta el parametro "IndiceRefracaoSolo" en system%s.json', paramSuffix);
end
n = strSystem.IndiceRefracaoSolo;
c = strSystem.VelocidadeLuz;
threshold = 1e-10;

if ~isfield(strSystem, 'TangentePerdaSolo')
    error('Falta el parametro "TangentePerdaSolo" en system%s.json', paramSuffix);
end
lossTangent = strSystem.TangentePerdaSolo;

% Create dataset for all points in the linePoints as targets

% Initialize accumulation matrix
auxData = zeros(strSystem.IndiceMaximo,length(PxT));

% RCS for each point (fraction of original target RCS)
pointRCS = strTarget.rcs / size(linePoints,1);

fprintf('Processing %d points in line...\n', size(linePoints,1));
tic

SNRtoProc = [strRadarRx.SNRrange(1):5:strRadarRx.SNRrange(2)];

% Process each point in linePoints as a small reflector (done once)
for pointIdx = 1:size(linePoints,1)
    if mod(pointIdx, 50) == 0 || pointIdx == size(linePoints,1)
        fprintf('Processing point %d/%d\n', pointIdx, size(linePoints,1));
    end

    % Current point coordinates
    x0 = linePoints(pointIdx,1);
    y0 = linePoints(pointIdx,2);
    z0 = linePoints(pointIdx,3);

    % Calculate slant ranges for current point
    [r1Tx,r2Tx,~] = calculateSlantRange(PxT,PyT,PzT,x0,y0,z0,c,n,threshold);
    [r1Rx,r2Rx,~] = calculateSlantRange(PxR,PyR,PzR,x0,y0,z0,c,n,threshold);

    % Calculate timing and phase for current point
    t = (1/c)*(r1Tx + n*r2Tx + r1Rx + n*r2Rx);
    rngBin = 1 + round(t*strRadarRx.fs);
    phi = -(2*pi/strRadarTx.lamb)*(r1Tx + n*r2Tx + r1Rx + n*r2Rx);

    % Atenuacion dielectrica por el trayecto recorrido dentro del suelo
    % (ida en Tx + ida en Rx), no considerada en el timing/fase anteriores
    atenTx = calculateSoilAttenuation(r2Tx, strRadarTx.FreqPortadora, n, lossTangent);
    atenRx = calculateSoilAttenuation(r2Rx, strRadarTx.FreqPortadora, n, lossTangent);
    atenSolo = atenTx .* atenRx;

    % Accumulate response for current point
    validIndices = rngBin > 0 & rngBin <= strSystem.IndiceMaximo;
    IND = sub2ind(size(auxData), rngBin(validIndices), find(validIndices));

    auxData(IND) = auxData(IND) + pointRCS.*atenSolo(validIndices).*exp(1i*phi(validIndices));
end

chirp   = createCHRIP(strRadarTx);
pulsoTx = zeros(1,strSystem.IndiceMaximo);
pulsoTx(1:length(chirp)) = chirp;

K = fix(length(chirp)/2);

reference    = circshift([pulsoTx],K).';
rawData      = ifft(fft(auxData).*fft(reference));

% Compression (computed once)
rootData = ifft(fft(rawData).*conj(fft(reference)));

toc

for j = 1:length(SNRtoProc)
    % Add noise to compressed data
    noisePower = 10^(SNRtoProc(j)/10) ; % Adjust noise power as needed
    noise = sqrt(noisePower/2) * (randn(size(rootData)) + 1i*randn(size(rootData)));
    rootDataNoisy = rootData + noise;
    
    % Plot results after processing
    figResult = figure;
    hold on
    mesh(abs(rootDataNoisy));
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    saveas(figResult, fullfile(figOutputDir, strcat('resultado_snr_', num2str(SNRtoProc(j)), '.png')));
    
    
    % Save raw data for volumetric Snell processing
    outputFile = fullfile(rawOutputDir, sprintf('bistatic_sar_snell_volumetrico_raw%s_snr_%ddB.mat', paramSuffix, SNRtoProc(j)));
    
    % Save line points information as well
    save(outputFile, 'rootData', 'rootDataNoisy', 'linePoints', 'targetMask', 'n', 'PxT', 'PyT', 'PzT', 'PxR', 'PyR', 'PzR', 'X', 'Y', 'Z','noisePower');
    fprintf('Raw volumetrico guardado en: %s\n', outputFile);
    
end

end