function report3dBResolution(X, Y, Z, outputData, csvFilePath)
thresholdDB = -3;
normalizedMagnitude = abs(outputData) ./ max(abs(outputData(:)));
dataDB = 20*log10(max(normalizedMagnitude, eps));

[~, peakLinearIdx] = max(normalizedMagnitude(:));
[peakYIdx, peakXIdx, peakZIdx] = ind2sub(size(normalizedMagnitude), peakLinearIdx);

xAxis = squeeze(X(1,:,1));
yAxis = squeeze(Y(:,1,1));
zAxis = squeeze(Z(1,1,:));

xProfileDB = squeeze(dataDB(peakYIdx,:,peakZIdx));
yProfileDB = squeeze(dataDB(:,peakXIdx,peakZIdx));
zProfileDB = squeeze(dataDB(peakYIdx,peakXIdx,:));

[xResolution, xLimits] = compute3dBWidth(xAxis, xProfileDB, peakXIdx, thresholdDB);
[yResolution, yLimits] = compute3dBWidth(yAxis, yProfileDB, peakYIdx, thresholdDB);
[zResolution, zLimits] = compute3dBWidth(zAxis, zProfileDB, peakZIdx, thresholdDB);

disp(sprintf('  Pico detectado en: X = %.4f m, Y = %.4f m, Z = %.4f m', ...
    X(peakLinearIdx), Y(peakLinearIdx), Z(peakLinearIdx)));
disp(sprintf('  Resolucion X (-3 dB): %.4f m [%.4f, %.4f]', xResolution, xLimits(1), xLimits(2)));
disp(sprintf('  Resolucion Y (-3 dB): %.4f m [%.4f, %.4f]', yResolution, yLimits(1), yLimits(2)));
disp(sprintf('  Resolucion Z (-3 dB): %.4f m [%.4f, %.4f]', zResolution, zLimits(1), zLimits(2)));

if nargin >= 5 && strlength(string(csvFilePath)) > 0
    axisLabel = {'X'; 'Y'; 'Z'};
    resolutionM = [xResolution; yResolution; zResolution];
    lowerLimitM = [xLimits(1); yLimits(1); zLimits(1)];
    upperLimitM = [xLimits(2); yLimits(2); zLimits(2)];
    peakX = repmat(X(peakLinearIdx), 3, 1);
    peakY = repmat(Y(peakLinearIdx), 3, 1);
    peakZ = repmat(Z(peakLinearIdx), 3, 1);
    thresholdDBColumn = repmat(thresholdDB, 3, 1);

    resolutionTable = table(axisLabel, resolutionM, lowerLimitM, upperLimitM, ...
        peakX, peakY, peakZ, thresholdDBColumn, ...
        'VariableNames', {'Axis', 'Resolution_m', 'LowerLimit_m', 'UpperLimit_m', ...
        'PeakX_m', 'PeakY_m', 'PeakZ_m', 'Threshold_dB'});

    writetable(resolutionTable, csvFilePath);
    fprintf('  Resoluciones guardadas en CSV: %s\n', csvFilePath);
end
end