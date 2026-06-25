function PROC_snell_volumetrico_test_script(paramSuffix)
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

figOutputDir = fullfile(projectDir, 'io', 'snell_volumetrico', 'figures', 'proc');
if ~exist(figOutputDir, 'dir')
    mkdir(figOutputDir);
end

rawInputDir = fullfile(projectDir, 'io', 'snell_volumetrico', 'raw');
annDataDir  = fullfile(projectDir, 'io', 'snell_volumetrico', 'ANNdata');

if ~exist(annDataDir, 'dir')
    mkdir(annDataDir);
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

% Calculate params

strRadarTx.lamb = strSystem.VelocidadeLuz/strRadarTx.FreqPortadora; % Comprimento de onda


%

SNRtoProc = [strRadarRx.SNRrange(1):5:strRadarRx.SNRrange(2)];


%% Load

for j = 1:length(SNRtoProc)
    inputFile = fullfile(rawInputDir, sprintf('bistatic_sar_snell_volumetrico_raw%s_snr_%ddB.mat', paramSuffix, SNRtoProc(j)));
    if ~exist(inputFile, 'file')
        error('No se encontro el archivo de entrada: %s. Ejecuta primero GS_snell_volumetrico_test_script.m', inputFile);
    end
    load(inputFile)
    
    x0  = mean(linePoints(:,1));
    y0  = mean(linePoints(:,1));
    
    % xy0 = strTarget.pos(1); dxy = 0.04;  lxy = 0.12;
    z0  = mean(linePoints(:,3)); %  dz = 0.08;   lz = 0.24;
    
    %% Preproc
    
    [Nrng,Nazm] = size(rootDataNoisy);
    ratio = 1;
    
    rcompData = interpft(rootDataNoisy,Nrng*ratio,1);
    rcompData = rcompData.';
    fs = ratio*strRadarTx.fs;
    
    clear rootData rootDataNoisy
    
    %% Processamiento
    
    tic
    
    wb = waitbar(0);
    threshold = 1e-11;
    c = strSystem.VelocidadeLuz;
    outputData = zeros(size(X));
    
    n2 = n;
    totalPoints = numel(X);
    logStep = max(1, round(totalPoints / 50)); % Log every 2% of progress
    
    for m = 1:totalPoints
        
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
        
        if mod(m, logStep) == 0 || m == totalPoints
            fprintf('  Back-Projection: %d/%d puntos (%.0f%%) - SNR %d dB\n', m, totalPoints, 100*m/totalPoints, SNRtoProc(j));
        end
        waitbar(m/totalPoints, wb, sprintf('Back-Projection SNR %d dB: %.0f%%', SNRtoProc(j), 100*m/totalPoints));
    end
    close(wb)
    toc
    
    %% Plots

    timestamp = datestr(now,'yyyymmdd_HHMMSS');
    snrSuffix = num2str(SNRtoProc(j));
    
    figIso = figure;
    dataDB = 20*log10(abs(outputData/max(outputData(:))));
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    % contourslice(X,Y,Z,dataDB,[],[],-1.1)
    patch(isosurface(X,Y,Z,dataDB,-3),...
        'EdgeColor','none','CData',-3,'FaceColor','flat');
    axis equal
    
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
    saveas(figIso, fullfile(figOutputDir, strcat(timestamp, '_isosurface_', snrSuffix, '.png')))
    
    % -------------------------------------------------------------------
    
    figSlice = figure;
    
    s = slice(X,Y,Z,abs(outputData/max(outputData(:))),x0,y0,z0);
    set(s,'EdgeColor','none')
    % title({'Output Data','abs()'})
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    axis equal
    colormap parula
    view(0,90)
    saveas(figSlice, fullfile(figOutputDir, strcat(timestamp, '_slice_', snrSuffix, '.png')))
    
    
    % -------------------------------------------------------------------
    
    %% Save data to ANN
    
    % Save outputData on .mat file with version -v7
    save(fullfile(annDataDir, strcat(timestamp,'_outputData_',snrSuffix,'.mat')),'outputData','-v7')
    
    % Save X,Y,Z matrizes on .mat file with version -v7
    save(fullfile(annDataDir, strcat(timestamp,'_XYZ_',snrSuffix,'.mat')),'X','Y','Z','-v7')
    
    % Save targetMask on .mat file with version -v7
    save(fullfile(annDataDir, strcat(timestamp,'_targetMask_',snrSuffix,'.mat')),'targetMask','-v7')
    
    % Save platform positions on .mat file with version -v7
    save(fullfile(annDataDir, strcat(timestamp,'_platformPositions_',snrSuffix,'.mat')),'PxT','PyT','PzT','PxR','PyR','PzR','-v7')
    
    % Save metadata on .json file
    metadata.timestamp = timestamp;
    metadata.n = n;
    metadata.noisePower = noisePower;
    metadata.fs = fs;
    metadata.FreqPortadora = strRadarTx.FreqPortadora;
    
    % Save metadata to JSON file using built-in MATLAB functions
    jsonStr = jsonencode(metadata);
    fileName = fullfile(annDataDir, strcat(timestamp,'_metadata_',snrSuffix,'.json'));
    fid = fopen(fileName, 'w');
    fwrite(fid, jsonStr, 'char');
    fclose(fid);
    
    % Elapsed time is 5460.676299 seconds.
    
end

