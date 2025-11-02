%% run_demo.m
% Quick launcher for the FBG SHM GUI

rng(1);  % for reproducibility

% Add source code to path
addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'..','src')));

% Check environment (optional)
if exist('checkEnvironment','file')
    checkEnvironment();
end

% Launch main GUI
if exist('fbg_shm_gui.m','file')
    disp('Launching FBG SHM GUI...');
    fbg_shm_gui;   % <---- this is your real main entry
else
    warning('fbg_shm_gui.m not found; running headless simulation instead.');
    runSimulation_v2; % fallback
end
