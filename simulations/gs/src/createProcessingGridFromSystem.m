function [X, Y, Z, xAxis, yAxis, zAxis] = createProcessingGridFromSystem(strSystem, strTarget, gridFieldName)
% Build processing grid centered on target using grid parameters in system JSON.
if ~isfield(strSystem, gridFieldName)
    error('No se encontro la configuracion de grid ''%s'' en system JSON.', gridFieldName);
end

gridCfg = strSystem.(gridFieldName);
requiredFields = {'dxy', 'lxy', 'dz', 'lz'};
for k = 1:numel(requiredFields)
    if ~isfield(gridCfg, requiredFields{k})
        error('Falta el parametro ''%s'' en la configuracion ''%s'' del system JSON.', requiredFields{k}, gridFieldName);
    end
end

x0 = strTarget.pos(1);
y0 = strTarget.pos(2);
z0 = strTarget.pos(3);

xAxis = x0-gridCfg.lxy:gridCfg.dxy:x0+gridCfg.lxy;
yAxis = y0-gridCfg.lxy:gridCfg.dxy:y0+gridCfg.lxy;
zAxis = z0-gridCfg.lz:gridCfg.dz:z0+gridCfg.lz;

[X, Y, Z] = meshgrid(xAxis, yAxis, zAxis);
end
