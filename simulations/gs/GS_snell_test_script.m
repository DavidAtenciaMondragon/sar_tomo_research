function GS_snell_test_script(paramSuffix)

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
                              strRadarTx.AltMaiorEsp,strRadarTx.AltMenorEsp, strRadarTx.Vt, strRadarTx.PRF, strRadarTx.NorthOffset);

% Receptor
[PxR, PyR, PzR] = funcao_espiral(strRadarRx.NumVoltasEsp, strRadarRx.RaioMenorEsp, strRadarRx.RaioMaiorEsp,...
                              strRadarRx.AltMaiorEsp,strRadarRx.AltMenorEsp, strRadarRx.Vt, strRadarRx.PRF, strRadarRx.NorthOffset);

% Plot Trajectory

[Xg,Yg] = meshgrid([-strRadarTx.RaioMaiorEsp strRadarTx.RaioMaiorEsp]);
Zg = zeros(size(Xg));

figTrajectory = figure;
hold on 
plot3(PxT, PyT, PzT,'k')
plot3(PxR, PyR, PzR,'r')
plot3(strTarget.pos(1),strTarget.pos(2),strTarget.pos(3),'o','MarkerFaceColor','k');
surf(Xg,Yg,Zg,'EdgeColor','none','FaceColor',[1 0 1],'FaceAlpha',0.1)
xlabel('X')
ylabel('Y')
zlabel('Z')
grid minor 
hold off
legend("Tx","Rx","Target")
saveas(figTrajectory, fullfile(figOutputDir,'trayectorias_snell.png'));

%% Create grid 
[X,Y,Z,xAxis,yAxis,zAxis] = createProcessingGridFromSystem(strSystem, strTarget, 'gridSnell');

%% Create raw data

if ~isfield(strSystem, 'IndiceRefracaoSolo')
    error('Falta el parametro "IndiceRefracaoSolo" en system%s.json', paramSuffix);
end
n = strSystem.IndiceRefracaoSolo;
c = strSystem.VelocidadeLuz;
threshold = 1e-10;

x0 = strTarget.pos(1);
y0 = strTarget.pos(2);
z0 = strTarget.pos(3);

tic
[r1Tx,r2Tx,angTxTg] = calculateSlantRange(PxT,PyT,PzT,x0,y0,z0,c,n,threshold);
[r1Rx,r2Rx,angTgRx] = calculateSlantRange(PxR,PyR,PzR,x0,y0,z0,c,n,threshold);
toc

figRanges = figure;
hold on 
plot(r1Tx)
plot(r2Tx)
plot(r1Rx)
plot(r2Rx)
legend("r1Tx","r2Tx","r1Rx","r2Rx")
grid minor
hold off 
saveas(figRanges, fullfile(figOutputDir,'slant_ranges_snell.png'));


% Create raw matrix 

t       = (1/c)*(r1Tx + n*r2Tx + r1Rx + n*r2Rx);
rngBin  = 1 + round(t*strRadarRx.fs);

phi     = -(2*pi/strRadarTx.lamb)*(r1Tx + n*r2Tx + r1Rx + n*r2Rx);

auxData = zeros(strSystem.IndiceMaximo,length(PxT));
IND     = sub2ind(size(auxData),rngBin,1:length(PxT));

auxData(IND) = strTarget.rcs.*exp(1i*phi);

chirp   = createCHRIP(strRadarTx);
pulsoTx = zeros(1,strSystem.IndiceMaximo);

pulsoTx(1:length(chirp)) = chirp;


K = fix(length(chirp)/2);

reference    = circshift([pulsoTx],K).';
rawData      = ifft(fft(auxData).*fft(reference));

% Compression

rootData = ifft(fft(rawData).*conj(fft(reference)));

% Save raw data into dedicated IO folder
rawDataDir = fullfile(projectDir,'io','snell','raw');
if ~exist(rawDataDir,'dir')
    mkdir(rawDataDir);
end

rawDataFile = fullfile(rawDataDir, sprintf('bistatic_sar_snell_raw%s.mat', paramSuffix));
save(rawDataFile, 'rootData', 'x0', 'y0', 'z0', 'n', 'PxT', 'PyT', 'PzT', 'PxR', 'PyR', 'PzR', 'X', 'Y', 'Z');
fprintf('Raw data guardada en: %s\n', rawDataFile);

end
