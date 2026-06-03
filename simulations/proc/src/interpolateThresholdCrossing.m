function crossing = interpolateThresholdCrossing(x1, x2, y1, y2, threshold)
if y1 == y2
    crossing = (x1 + x2) / 2;
    return
end

crossing = x1 + (threshold - y1) * (x2 - x1) / (y2 - y1);
end