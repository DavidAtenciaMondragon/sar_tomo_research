function [rngAbove,rngBelow] = slantRange(xRadar,yRadar,zRadar,...
                                          xPixel,yPixel,zPixel,...
                                          refrIndx,threshold)
% slantRange    Calculates the slant ranges above and below ground between
% the radar and a pixel. 
%
% Inputs:
% xRadar, yRadar, zRadar --> Radar position in cartesian coordinates
% xPixel,yPixel,zPixel   --> Pixel position in cartesian coordinates
% refrIndx               --> The ground's refraction index
% threshold              --> Convergence threshold in terms of return time
%
% Outputs:
% rngAbove --> Slant range above ground
% rngBelow --> Slant range below ground

% Initial guess ------------------------------------------------------------
u = sqrt((xRadar-xPixel).^2 + (yRadar-yPixel).^2);
u1  = u./(1-zPixel./(refrIndx*zRadar));
u2 = u - u1;
% --------------------------------------------------------------------------

% Numeric solution ---------------------------------------------------------
c = physconst('LightSpeed');
error = 1;
while error > threshold
    rngAbove = sqrt(u1.^2+zRadar.^2);
    rngBelow = sqrt(u2.^2+zPixel.^2);
    dRng = rngBelow - rngAbove.*refrIndx.*u2./u1;
    error = norm(2*refrIndx*dRng/c);
    u1  = u./(1+rngBelow./(refrIndx*rngAbove));
    u2 = u - u1;
end
% --------------------------------------------------------------------------
end