function saveStatistics(statsOutLoc, cohort, task, timepoint, params, homerPruneMat, sciAvgs, pspAvgs, motionWindows, channelsRetained, snrMat, saveFlags)
    
% Saves the calculated statistics and quality metrics from all subjects to 
% files, including Homer pruning matrices, SCI averages, PSP averages, and 
% motion window information.

    % Create required directories
    if ~exist(fullfile(statsOutLoc, cohort, 'overall', task, 'prune'), 'dir')
        mkdir(fullfile(statsOutLoc, cohort, 'overall', task, 'prune'));
    end
    
    if ~exist(fullfile(statsOutLoc, cohort, timepoint, task, 'prune'), 'dir')
        mkdir(fullfile(statsOutLoc, cohort, timepoint, task, 'prune'));
    end
    
    % Save Homer pruning matrix if needed
    if isfield(saveFlags, 'saveHomerPrune')
        save(fullfile(statsOutLoc, cohort, 'overall', task, 'prune', [task timepoint 'HomerPrunedChannels']), 'homerPruneMat');
    end
    
    % Save SCI averages if needed
    if isfield(saveFlags, 'saveAvgSCI')
        save(fullfile(statsOutLoc, cohort, 'overall', task, 'prune', [task timepoint 'AvgSCImat']), 'sciAvgs');
    end
    
    % Save PSP averages if needed
    if isfield(saveFlags, 'saveAvgPSP')
        save(fullfile(statsOutLoc, cohort, 'overall', task, 'prune', [task timepoint 'AvgPSPmat']), 'pspAvgs');
    end
    
    % Save motion windows if needed
    if isfield(saveFlags, 'saveMotion')
        save(fullfile(statsOutLoc, cohort, 'overall', task, 'prune', [task timepoint 'MotionWindowsMat']), 'motionWindows');
    end
    
    % Set file names based on pruning method
    if params.pruneQT == 1
        % Remove periods from paramsAppend to save files
        paramsAppend = strrep(params.pruneName, '.', '');
        paramsAppend = [paramsAppend '_SCI' num2str(params.sci_threshold) '_PSP' num2str(params.psp_threshold)];
        paramsAppend = strrep(paramsAppend, '.',