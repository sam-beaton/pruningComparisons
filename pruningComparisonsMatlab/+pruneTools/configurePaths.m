function [inputPath, outputsFolderPath] = configurePaths(dataRawLoc, dataOutLoc, cohort, timepoint, task, params)

% Configures the input and output file paths based on the cohort, 
% timepoint, task, and pruning method. Creates output directories if they 
% don't exist, ensuring proper file organization throughout the processing 
% pipeline.

    % Configure input path
    inputPath = fullfile(dataRawLoc, cohort, timepoint, task);
    
    % Configure output path based on pruning method
    if params.pruneQT == 1
        paramsAppend = strcat(params.pruneName, '_SCI', num2str(params.sciThreshold), '_PSP', num2str(params.pspThreshold));
        outputsFolderPath = fullfile(dataOutLoc, cohort, timepoint, task, 'prune', paramsAppend);
    else
        outputsFolderPath = fullfile(dataOutLoc, cohort, timepoint, task, 'prune/pruneCV08');
    end
    
    % Create output directory if it doesn't exist
    if ~exist(outputsFolderPath, 'dir')
        mkdir(outputsFolderPath);
    end
end