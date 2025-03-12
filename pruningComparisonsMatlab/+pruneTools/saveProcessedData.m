function saveProcessedData(nirs, partName, outputsFolderPath, params)

% Saves the processed NIRS data after pruning to an output file with a 
% name that indicates the pruning method and parameters used.

    if params.pruneQT == 1
        paramsAppend = strcat(params.pruneName, '_SCI', num2str(params.sci_threshold), '_PSP', num2str(params.psp_threshold));
        pruneFileName = fullfile([partName(1:end-5) '_' paramsAppend '.nirs']);
    else
        pruneFileName = fullfile([partName(1:end-5) '_PruneCV.nirs']);
    end
    pruneFileName = fullfile(outputsFolderPath, pruneFileName);
    
    % Uncomment to save actual processed data
    % save(pruneFileName, '-struct', 'nirs');
end
