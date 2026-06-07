function tg = createGridTarget(xlim,ylim,zlimMin,zlimMax,elemPerDimX, elemPerDimY, elemPerDimZ)

x = linspace(-xlim,xlim,elemPerDimX);
y = linspace(-ylim,ylim,elemPerDimY);
z = linspace(zlimMin,zlimMax,elemPerDimZ);

[X,Y,Z] = meshgrid(x,y,z);

tg = [X(:),Y(:),Z(:)];

end