function [scriptDir, projectDir] = setupGSWorkspace(currentScriptFullPath)
% Configure common paths for GS scripts using the caller script location.
scriptDir = fileparts(currentScriptFullPath);
projectDir = fileparts(scriptDir);

addpath(genpath(fullfile(scriptDir,'src')))
addpath(genpath(fullfile(projectDir,'tools')))
addpath(genpath(fullfile(projectDir,'common')))
end
