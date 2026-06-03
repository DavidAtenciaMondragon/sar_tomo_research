function [width, limits] = compute3dBWidth(axisValues, profileDB, peakIdx, thresholdDB)
axisValues = axisValues(:);
profileDB = profileDB(:);

if profileDB(peakIdx) < thresholdDB
    width = NaN;
    limits = [NaN NaN];
    return
end

leftBelowIdx = find(profileDB(1:peakIdx) < thresholdDB, 1, 'last');
if isempty(leftBelowIdx)
    leftCrossing = axisValues(1);
else
    leftCrossing = interpolateThresholdCrossing(axisValues(leftBelowIdx), axisValues(leftBelowIdx + 1), ...
        profileDB(leftBelowIdx), profileDB(leftBelowIdx + 1), thresholdDB);
end

rightBelowRelativeIdx = find(profileDB(peakIdx:end) < thresholdDB, 1, 'first');
if isempty(rightBelowRelativeIdx)
    rightCrossing = axisValues(end);
else
    rightBelowIdx = peakIdx + rightBelowRelativeIdx - 1;
    rightCrossing = interpolateThresholdCrossing(axisValues(rightBelowIdx - 1), axisValues(rightBelowIdx), ...
        profileDB(rightBelowIdx - 1), profileDB(rightBelowIdx), thresholdDB);
end

limits = [leftCrossing rightCrossing];
width = abs(rightCrossing - leftCrossing);
end
