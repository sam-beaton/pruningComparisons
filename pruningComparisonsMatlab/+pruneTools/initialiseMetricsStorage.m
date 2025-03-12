function [homerPruneMat, sciAvgs, pspAvgs, motionWindows, channelsRetained, snrMat, saveFlags] = initialiseMetricsStorage(statsOutLoc, cohort, task, timepoint, sub, params)
    
% Initialises data structures for storing quality metrics, pruning results,
% and SNR values across subjects. Creates empty matrices sized 
% appropriately for the dataset and tracks which metrics need to be saved.

    % Load first subject to get channel count
    firstData = load(sub(1).name, '-mat');
    if ~isfield(firstData.SD, 'MeasListAct')
        firstData.SD.MeasListAct = firstData.SD.MeasList(:,3);
    end
    numChansBothChroms = length(firstData.SD.MeasListAct);
    
    % Initialise save flags structure
    saveFlags = struct();
    
    % Initialise homer pruning matrix if needed
    if ~exist(fullfile(statsOutLoc, cohort, 'overall', task, 'prune', [task timepoint 'HomerPrunedChannels.mat']), 'file') && params.pruneQT==1
        homerPruneMat = zeros([numChansBothChroms length(sub)]);
        saveFlags.saveHomerPrune = 1;
    else
        homerPruneMat = [];
    end
    
    % Initialise SCI averages matrix if needed
    if ~exist(fullfile(statsOutLoc, cohort, 'overall', task, 'prune', [task timepoint 'avgSCImat.mat']), 'file') && params.pruneQT==1
        sciAvgs = nan([numChansBothChroms/2 length(sub)]);
        saveFlags.saveAvgSCI = 1;
    else
        sciAvgs = [];
    end
    
    % Initialise PSP averages matrix if needed
    if ~exist(fullfile(statsOutLoc, cohort, 'overall', task, 'prune', [task timepoint 'AvgPSPmat.mat']), 'file') && params.pruneQT==1
        pspAvgs = nan([numChansBothChroms/2 length(sub)]);
        saveFlags.saveAvgPSP = 1;
    else
        pspAvgs = [];
    end
    
    % Initialise motion windows matrix if needed
    if ~exist(fullfile(statsOutLoc, cohort, 'overall', task, 'prune', [task timepoint 'motionWindowsMat.mat']), 'file') && params.pruneQT==1
        motionWindows = nan([numChansBothChroms length(sub)]);
        saveFlags.saveMotion = 1;
    else
        motionWindows = [];
    end
    
    % Initialise channels retained and SNR matrices
    channelsRetained = nan([length(sub) 1]);
    
    % Get ROI channels for SNR calculation
    roiChans = pruneTools.getROIChannels(task, timepoint);
    snrMat = nan([length(sub) 2*length(roiChans)]); % 2x for both chromophores
end