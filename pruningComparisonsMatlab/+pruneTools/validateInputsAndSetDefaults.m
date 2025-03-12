function params = validateInputsAndSetDefaults(params)

% Validates required inputs for data processing, checking for correct data 
% types. Creates a params structure with default values for optional 
% parameters. 

    % Check required fields exist
    requiredFields = {'dataRawLoc', 'dataOutLoc', 'statsOutLoc', 'cohort', 'timepoint', 'task', 'sciThreshold', 'pspThreshold'};
    
    for i = 1:length(requiredFields)
        if ~isfield(params, requiredFields{i})
            error('Required field "%s" is missing from params structure', requiredFields{i});
        end
    end
    
    % Check input types
    if ~ischar(params.dataRawLoc)
        error('Field "dataRawLoc" must be a character input');
    end
    
    if ~ischar(params.dataOutLoc)
        error('Field "dataOutLoc" must be a character input');
    end
    
    if ~ischar(params.cohort)
        error('Field "cohort" must be a character input');
    end
    
    if ~ischar(params.timepoint)
        error('Field "timepoint" must be a character input');
    end
    
    if ~ischar(params.task)
        error('Field "task" must be a character input');
    end
    
    % Set default values for optional parameters
    if ~isfield(params, 'pruneQT')
        params.pruneQT = true;
    end
    
    % Set pruning method name
    if params.pruneQT
        params.pruneName = 'PruneQT';
    else 
        params.pruneName = 'PruneCV';
    end
    
    % Return the validated and completed params structure
end