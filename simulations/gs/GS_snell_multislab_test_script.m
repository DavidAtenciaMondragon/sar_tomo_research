function GS_snell_multislab_test_script(paramSuffix)

scriptDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(scriptDir);

addpath(genpath(fullfile(scriptDir,'src')))
addpath(genpath(fullfile(projectDir,'tools')))
addpath(genpath(fullfile(projectDir,'common')))

% Select parameter set suffix: '_multislab'
if nargin < 1 || strlength(string(paramSuffix)) == 0
    paramSuffix = '_multislab';
else
    paramSuffix = char(string(paramSuffix));
end

figOutputDir = fullfile(projectDir,'io','snell_multislab','figure');
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

% Capas de suelo (multicapa)

if ~isfield(strSystem, 'CapasSuelo')
    error('Falta el parametro "CapasSuelo" en system%s.json', paramSuffix);
end
n_layers = [strSystem.CapasSuelo.n];
d_layers = [strSystem.CapasSuelo.d];

x0 = strTarget.pos(1);
y0 = strTarget.pos(2);
z0 = strTarget.pos(3);

%% Create trajectory

% Transmissor
[PxT, PyT, PzT] = funcao_espiral(strRadarTx.NumVoltasEsp, strRadarTx.RaioMenorEsp, strRadarTx.RaioMaiorEsp,...
                              strRadarTx.AltMaiorEsp,strRadarTx.AltMenorEsp, strRadarTx.Vt, strRadarTx.PRF, strRadarTx.NorthOffset);

% Receptor
[PxR, PyR, PzR] = funcao_espiral(strRadarRx.NumVoltasEsp, strRadarRx.RaioMenorEsp, strRadarRx.RaioMaiorEsp,...
                              strRadarRx.AltMaiorEsp,strRadarRx.AltMenorEsp, strRadarRx.Vt, strRadarRx.PRF, strRadarRx.NorthOffset);

% Plot Trajectory

[Xg,Yg] = meshgrid([-strRadarTx.RaioMaiorEsp strRadarTx.RaioMaiorEsp]);

cumDepth     = cumsum(d_layers);
nInterfaces  = numel(cumDepth);
slabColors   = lines(nInterfaces+1);

figTrajectory = figure;
hold on
plot3(PxT, PyT, PzT,'k')
plot3(PxR, PyR, PzR,'r')
plot3(strTarget.pos(1),strTarget.pos(2),strTarget.pos(3),'o','MarkerFaceColor','k');

legendEntries = {"Tx","Rx","Target"};

% Superficie (interfaz aire / capa 1 de suelo)
surf(Xg,Yg,zeros(size(Xg)),'EdgeColor','none','FaceColor',slabColors(1,:),'FaceAlpha',0.15)
legendEntries{end+1} = sprintf('Superficie (n=1 | n=%.1f)', n_layers(1));

% Interfaces entre capas de suelo (y limite declarado de la ultima capa)
for i = 1:nInterfaces
    Zi = -cumDepth(i)*ones(size(Xg));
    surf(Xg,Yg,Zi,'EdgeColor','none','FaceColor',slabColors(i+1,:),'FaceAlpha',0.15)
    if i < nInterfaces
        legendEntries{end+1} = sprintf('Interfaz capa %d|%d (n=%.1f | n=%.1f)', i, i+1, n_layers(i), n_layers(i+1));
    else
        legendEntries{end+1} = sprintf('Limite declarado capa %d (n=%.1f, semi-infinita)', i, n_layers(i));
    end
end

xlabel('X')
ylabel('Y')
zlabel('Z')
grid minor
hold off
legend(legendEntries)
saveas(figTrajectory, fullfile(figOutputDir,'trayectorias_snell_multislab.png'));

% Perfil de capas (esquema del medio multicapa)

figProfile = figure;
hold on

profileWidth = max(cumDepth)*2;
xProfile = [-profileWidth profileWidth];

airHeight = max(d_layers);

% Aire
patch([xProfile(1) xProfile(2) xProfile(2) xProfile(1)], [0 0 airHeight airHeight], slabColors(1,:), 'FaceAlpha',0.25,'EdgeColor','none');
text(xProfile(1), airHeight/2, sprintf('  Aire (n=1)'), 'VerticalAlignment','middle');

% Capas de suelo
zTop = 0;
for i = 1:numel(n_layers)
    zBot = -cumDepth(i);
    patch([xProfile(1) xProfile(2) xProfile(2) xProfile(1)], [zTop zTop zBot zBot], slabColors(i+1,:), 'FaceAlpha',0.25,'EdgeColor','none');
    text(xProfile(1), (zTop+zBot)/2, sprintf('  Capa %d: n=%.1f, d=%.2f m', i, n_layers(i), d_layers(i)), 'VerticalAlignment','middle');
    zTop = zBot;
end

% Limite declarado de la ultima capa (semi-infinita en el modelo)
plot(xProfile, [zTop zTop], '--k');
text(xProfile(1), zTop - 0.1*airHeight, sprintf('  Capa %d semi-infinita en el modelo', numel(n_layers)), 'VerticalAlignment','top');

% Target
plot(0, z0, 'o', 'MarkerFaceColor','k', 'MarkerEdgeColor','k');
text(0, z0, '  Target', 'VerticalAlignment','middle');

xlabel('x (m, ilustrativo)')
ylabel('z (m)')
title('Perfil de capas del medio (escenario multicapa)')
grid minor
hold off
saveas(figProfile, fullfile(figOutputDir,'perfil_capas_multislab.png'));

%% Create grid
[X,Y,Z,xAxis,yAxis,zAxis] = createProcessingGridFromSystem(strSystem, strTarget, 'gridSnell');

%% Create raw data

c = strSystem.VelocidadeLuz;
threshold = 1e-10;

tic
[RsegTx,nSeg,thetaIncTx,L] = calculateSlantRangeMultislab(PxT,PyT,PzT,x0,y0,z0,n_layers,d_layers,threshold);
[RsegRx,~,   thetaIncRx]   = calculateSlantRangeMultislab(PxR,PyR,PzR,x0,y0,z0,n_layers,d_layers,threshold);
toc

fprintf('Target ubicado en la capa %d (n = %.3f).\n', L, n_layers(L));

% --- Recorrido optico total (sum n_i * R_i) por enlace ---
pathTx = sum(nSeg.*RsegTx, 1);
pathRx = sum(nSeg.*RsegRx, 1);

figRanges = figure;
hold on
plot(RsegTx(1,:))
plot(RsegTx(end,:))
plot(RsegRx(1,:))
plot(RsegRx(end,:))
legend("R_{aire,Tx}","R_{suelo,Tx}","R_{aire,Rx}","R_{suelo,Rx}")
grid minor
hold off
saveas(figRanges, fullfile(figOutputDir,'slant_ranges_snell_multislab.png'));

%% Coeficientes de transmision TM (multicapa) por enlace
% Capas atravesadas hasta llegar a la capa que contiene al target (slabs
% finitos); la capa del target actua como medio de salida semi-infinito.
n_slabs = n_layers(1:L-1);
d_slabs = d_layers(1:L-1);
n_out   = n_layers(L);

[~,Tau_Tx,~,T_Tx] = multislabTMcoef(strRadarTx.FreqPortadora, thetaIncTx, 1, n_slabs, d_slabs, n_out);
[~,Tau_Rx,~,T_Rx] = multislabTMcoef(strRadarTx.FreqPortadora, thetaIncRx, 1, n_slabs, d_slabs, n_out);

figTransm = figure;
hold on
plot(T_Tx)
plot(T_Rx)
legend("T_{TM,Tx}","T_{TM,Rx}")
ylabel('Transmitancia de potencia')
grid minor
hold off
saveas(figTransm, fullfile(figOutputDir,'transmitancia_TM_multislab.png'));

% Create raw matrix

t       = (1/c)*(pathTx + pathRx);
rngBin  = 1 + round(t*strRadarRx.fs);

phi     = -(2*pi/strRadarTx.lamb)*(pathTx + pathRx);

auxData = zeros(strSystem.IndiceMaximo,length(PxT));
IND     = sub2ind(size(auxData),rngBin,1:length(PxT));

auxData(IND) = strTarget.rcs.*Tau_Tx.*Tau_Rx.*exp(1i*phi);

chirp   = createCHRIP(strRadarTx);
pulsoTx = zeros(1,strSystem.IndiceMaximo);

pulsoTx(1:length(chirp)) = chirp;


K = fix(length(chirp)/2);

reference    = circshift([pulsoTx],K).';
rawData      = ifft(fft(auxData).*fft(reference));

% Compression

rootData = ifft(fft(rawData).*conj(fft(reference)));

% Save raw data into dedicated IO folder
rawDataDir = fullfile(projectDir,'io','snell_multislab','raw');
if ~exist(rawDataDir,'dir')
    mkdir(rawDataDir);
end

rawDataFile = fullfile(rawDataDir, sprintf('bistatic_sar_snell_multislab_raw%s.mat', paramSuffix));
save(rawDataFile, 'rootData', 'x0', 'y0', 'z0', 'n_layers', 'd_layers', 'L', 'PxT', 'PyT', 'PzT', 'PxR', 'PyR', 'PzR', 'X', 'Y', 'Z');
fprintf('Raw data guardada en: %s\n', rawDataFile);

end
