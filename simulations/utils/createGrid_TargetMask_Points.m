clc
clear 
close all 

% Load data 
load('ground_truth_centrada_flip_h.mat');

% Quedarme solo con los x desde 225 hasta el final y para los y desde 390 hasta el final
z_axis = 0:-0.5:-5;

% x_axis = x_axis(225:end);
% y_axis = y_axis(390:end);
% X = X(390:end,225:end,:);
% Y = Y(390:end,225:end,:);

% ground_truth = ground_truth(390:end,225:end);

targetMask = zeros([size(X), 11]);
targetMask(:,:,8) = ground_truth;

[X,Y,Z] = meshgrid(x_axis,y_axis,z_axis);

% linePoints debe contener los puntos X,Y,Z donde targetMask es 1
linePoints = [X(targetMask==1), Y(targetMask==1), Z(targetMask==1)];

xAxis = x_axis;
yAxis = y_axis;
zAxis = z_axis;

% Save linePoints and targetMask on .mat file with version -v7
save('linePoints_targetMask.mat','linePoints','targetMask','xAxis','yAxis','zAxis','-v7')

