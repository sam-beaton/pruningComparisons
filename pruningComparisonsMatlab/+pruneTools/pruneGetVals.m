function [meanChans, meanSNR] = pruneGetVals(params, varargin)
% prePPrune -   Prune fNIRS data for the BRIGHT study, storing key metrics
%               for analysis.
% 
% =========================================================================
%
% params:           provide data locations and distinguish data population 
% (Required)        Assumes data is stored in nested folders:
%                       cohort > timepoint > task
%
% dataRawLoc:       folder containing raw .nirs files to be processed
%
% dataOutLoc        folder in which the user wants the processed files to
%                   be stored
%
% statsOutLoc:      folder where parameter specific outcome measures will
%                   be stored for collation
%
% cohort:           'gm' or 'uk'
%
% timepoint:        defines age of participants; must include 'mo'
%                   INDiGO:
%                   '1mo', '6mo', '12mo'
%                   
% task:             name of task 
%                   INDiGO: 'hand', 'fc1' or 'fc2'
%
% sci_threshold:    SCI threshold parameter for QT-NIRS
%
% psp_threshold:    PSP threshold parameter for QT-NIRS
%
% -------------------------------------------------------------------------
% Arguments:
% (Optional)
%
% pruneSCI:         1: prune using QT-NIRS
%                   0: prune using coefficient of variation
%
% =========================================================================
% Outputs: summary statistics per participant for the number of channels
% included and the data quality (SNR) in ROI channels
%
% SLB 17/1/2024 - original implementation
% SLB 11/3/25 - modularised version


%% Input processing. Check input arguments and assign where necessary


    % Validate inputs and set default parameters
    params = pruneTools.validateInputsAndSetDefaults(params);

    % Extract parameters from the structure
    dataRawLoc = params.dataRawLoc;
    dataOutLoc = params.dataOutLoc;
    statsOutLoc = params.statsOutLoc;
    cohort = params.cohort;
    timepoint = params.timepoint;
    task = params.task;
    sciThreshold = params.sci_threshold;
    pspThreshold = params.psp_threshold;
    pruneQT = params.QT;

    
    % configure input/output paths
    [inputPath, outputsFolderPath] = pruneTools.configurePaths(dataRawLoc, dataOutLoc, cohort, timepoint, task, params);
    
    % Set method parameters (change in function itself - just keeps things tidier if modular)
    params = pruneTools.configureMethodParameters(params, timepoint);
    
    % Get subject files
    cd(inputPath);
    sub = dir;
    sub = sub(~ismember({sub.name}, {'.', '..', '.DS_Store'}));
    
    % Initialise metrics storage
    [homerPruneMat, sciAvgs, pspAvgs, motionWindows, channelsRetained, snrMat, saveFlags] = initialiseMetricsStorage(statsOutLoc, cohort, task, timepoint, sub, params);
    
    % Get ROI channels for SNR calculation
    roiChans = pruneTools.getROIChannels(task, timepoint);
    
    % Process subjects
    for nsub = 1:length(sub)

        % Load subject data
        [nirs, partName] = pruneTools.loadSubjectData(sub(nsub).name);
        
        % Calculate DPF for age
        params.dpf = pruneTools.calculateDPF(timepoint, nirs);
        
        % Detect motion artifacts and prune channels
        [nirs, qualityMatrices, windowInfo] = pruneTools.pruneChannels(nirs, params);
        
        % Update metrics
        [homerPruneMat, sciAvgs, pspAvgs, motionWindows] = pruneTools.updateQualityMetrics(nirs, qualityMatrices, nsub, homerPruneMat, sciAvgs, pspAvgs, motionWindows, saveFlags, params);
        
        % Count retained channels
        channelsRetained(nsub, 1) = sum(nirs.SD.MeasListAct == 1);
        
        % Calculate SNR for ROI channels
        snrMat = pruneTools.calculateSNR(nirs, roiChans, windowInfo, snrMat, nsub, params);
        
        % save processed data
        pruneTools.saveProcessedData(nirs, partName, outputsFolderPath, params);
    end
    
    % Save summary stats
    pruneTools.saveStatistics(statsOutLoc, cohort, task, timepoint, params, homerPruneMat, sciAvgs, pspAvgs, motionWindows, channelsRetained, snrMat, saveFlags);
    
    % calculate and return means
    meanChans = mean(channelsRetained, "omitnan");
    meanSNR = mean(snrMat, "all", "omitnan");
end