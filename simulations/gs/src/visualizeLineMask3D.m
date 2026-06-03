function visualizeLineMask3D(targetMask, linePoints, xAxis, yAxis, zAxis, outputFile)
% Visualize 3D mask and target line, then save figure and print stats.
figMask = figure;
hold on

% Create coordinate arrays for all grid points (matching meshgrid order)
[meshY, meshX, meshZ] = ndgrid(yAxis, xAxis, zAxis);

% Find indices for background (zeros) and target (ones) points
backgroundIdx = targetMask(:) == 0;
targetIdx = targetMask(:) == 1;

% Convert to coordinate arrays
allX = meshX(:);
allY = meshY(:);
allZ = meshZ(:);

% Plot background points (zeros) in blue with high transparency
if sum(backgroundIdx) > 0
    % Sample background points to avoid overloading (every 8th point)
    sampleIdx = 1:8:sum(backgroundIdx);
    backgroundSample = find(backgroundIdx);
    backgroundSample = backgroundSample(sampleIdx);

    scatter3(allX(backgroundSample), allY(backgroundSample), allZ(backgroundSample), ...
        15, [0.2 0.2 1], 'filled', 'MarkerFaceAlpha', 0.1, 'MarkerEdgeAlpha', 0.2);
end

% Plot target points (ones) in red with low transparency
if sum(targetIdx) > 0
    targetPoints = find(targetIdx);
    scatter3(allX(targetPoints), allY(targetPoints), allZ(targetPoints), ...
        80, [1 0 0], 'filled', 'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 1.0);
end

% Plot original line for reference
plot3(linePoints(:,1), linePoints(:,2), linePoints(:,3), 'ko-', 'LineWidth', 2, 'MarkerSize', 4);

xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
title('3D Target Mask Visualization')
legend('Background grid (sampled)', 'Target points', 'Original line', 'Location', 'best')
grid on
axis equal
view(3)
hold off

saveas(figMask, outputFile);

% Show mask statistics
fprintf('Mask dimensions: %dx%dx%d\n', size(targetMask,1), size(targetMask,2), size(targetMask,3));
fprintf('Total grid points: %d\n', numel(targetMask));
fprintf('Target points: %d (%.2f%%)\n', sum(targetMask(:)), 100*sum(targetMask(:))/numel(targetMask));
end
