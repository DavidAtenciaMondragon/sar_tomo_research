function [linePoints, targetMask] = createLinePointsAndMask(p1, p2, dxy, xAxis, yAxis, zAxis)
% Generate 3D line sample points and build the corresponding grid mask.
if numel(p1) ~= 3 || numel(p2) ~= 3
    error('p1 y p2 deben tener exactamente 3 componentes [x y z].');
end
if dxy <= 0
    error('dxy debe ser un valor positivo.');
end

% Normalize vectors so downstream math is shape-safe across MATLAB versions.
p1 = reshape(p1, 1, 3);
p2 = reshape(p2, 1, 3);

lineDirection = p2 - p1;
lineLength = norm(lineDirection);

numPoints = round(lineLength / dxy) + 1;
t = linspace(0, 1, numPoints).';
linePoints = p1 + t * lineDirection;

targetMask = zeros(length(yAxis), length(xAxis), length(zAxis));
for i = 1:size(linePoints,1)
    [~, xIdx] = min(abs(xAxis - linePoints(i,1)));
    [~, yIdx] = min(abs(yAxis - linePoints(i,2)));
    [~, zIdx] = min(abs(zAxis - linePoints(i,3)));
    targetMask(yIdx, xIdx, zIdx) = 1;
end
end
