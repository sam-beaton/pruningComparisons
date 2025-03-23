function saveStatistics(statsOutLoc, cohort, task, timepoint, params, homerPruneMat, sciAvgs, pspAvgs, motionWindows, channelsRetained, snrMat, saveFlags)
    
% Saves the calculated statistics and quality metrics from all subjects to 
% files, including Homer pruning matrices, SCI averages, PSP averages, and 
% motion window information.

    overallDir = fullfile(statsOutLoc, cohort, 'overall', task, 'prune');
    timepointDir = fullfile(statsOutLoc, cohort, timepoint, task, 'prune');

    % Create required directories
    if ~exist(overallDir, 'dir')
        mkdir(overallDir);
    end
    
    if ~exist(timepointDir, 'dir')
        mkdir(timepointDir);
    end
    
    % Save Homer pruning matrix if needed
    if isfield(saveFlags, 'saveHomerPrune')
        save(fullfile(overallDir, [task timepoint 'HomerPrunedChannels']), 'homerPruneMat');
    end
    
    % Save SCI averages if needed
    if isfield(saveFlags, 'saveAvgSCI')
        save(fullfile(overallDir, [task timepoint 'AvgSCImat']), 'sciAvgs');
    end
    
    % Save PSP averages if needed
    if isfield(saveFlags, 'saveAvgPSP')
        save(fullfile(overallDir, [task timepoint 'AvgPSPmat']), 'pspAvgs');
    end
    
    % Save motion windows if needed
    if isfield(saveFlags, 'saveMotion')
        save(fullfile(overallDir, [task timepoint 'MotionWindowsMat']), 'motionWindows');
    end
    
    % Set file names based on pruning method
    if params.pruneQT == 1
        % Remove periods from paramsAppend to save files
        paramsAppend = strrep(params.pruneName, '.', '');
        paramsAppend = [paramsAppend '_SCI' num2str(params.sciThreshold) '_PSP' num2str(params.pspThreshold)];
    else
        paramsAppend = 'CV';
    end

    % Save retained channels
    save(fullfile(timepointDir, [task timepoint '_ChannelsRetained_' paramsAppend '.mat']), 'channelsRetained');
    
    % Save SNR matrix
    save(fullfile(timepointDir, [task timepoint '_SNRMat_' paramsAppend '.mat']), 'snrMat');
end