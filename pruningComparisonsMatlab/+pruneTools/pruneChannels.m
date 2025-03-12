function [nirs, qualityMatrices, windowInfo] = pruneChannels(nirs, params)

% Master function for pruning low-quality channels from the NIRS data. 
% Detects motion artifacts, applies initial channel pruning using Homer, 
% then either applies QT-NIRS or coefficient of variation (CV) pruning 
% methods based on the specified parameters.

    % Initialise windowInfo structure
    windowInfo = struct();
    
    % Detect motion artifacts using amp and std thresholds (by channel)
    [nirs.SD.tInc, nirs.SD.tIncCh] = hmrMotionArtifactByChannel(nirs.d, nirs.t, nirs.SD, [], params.tMotion, params.tMask, params.STDEVthresh, params.AMPthresh);

    % Prune channels using Homer prior to QT-NIRS
    nirs.SD = enPruneChannels(nirs.d, nirs.SD, nirs.SD.tInc, params.dRange, params.SNRthresh, params.SDrange, 0);

    % Get number of channels (both chromophores)
    numChan = length(nirs.SD.MeasListAct)/2;
    arrayChans = linspace(1, 2*numChan, 2*numChan);
    
    if params.pruneQT == 1
        % QT-NIRS pruning
        [nirs, qualityMatrices] = pruneTools.pruneWithQT(nirs, params, numChan);
        
        % Extract window information for SNR calculation
        windowInfo.windowSizeSamples = qualityMatrices.sampPerWindow;
        windowInfo.windowCount = qualityMatrices.n_windows;
        windowInfo.windowMotion = qualityMatrices.valid_windows;
        windowInfo.windowTimes = qualityMatrices.window_times;
    else
        % CV pruning
        [nirs, windowMotion, windowTimes, windowCount, windowSizeSamples] = pruneTools.pruneWithCV(nirs, params, numChan, arrayChans);
        
        % Set window information for SNR calculation
        windowInfo.windowSizeSamples = windowSizeSamples;
        windowInfo.windowCount = windowCount;
        windowInfo.windowMotion = windowMotion;
        windowInfo.windowTimes = windowTimes;
        
        % Set empty qualityMatrices for CV method
        qualityMatrices = [];
    end
end