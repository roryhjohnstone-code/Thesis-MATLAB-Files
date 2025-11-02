function varargout = runSimulation(injectFaults,useML,useKalman,faultType,axSensors,axConsensus,ax3D,doAnimate,doExport,fig)
% Thin wrapper: call v2 with a fixed signature only.

    if exist('runSimulation_v2','file') ~= 2
        error('Run failed: runSimulation_v2 not found on path.');
    end

    % Call-through (matches your GUI/script usage)
    [varargout{1:nargout}] = runSimulation_v2( ...
        injectFaults, useML, useKalman, faultType, ...
        axSensors, axConsensus, ax3D, doAnimate, doExport, fig);
end

