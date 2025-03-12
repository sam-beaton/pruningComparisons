function [nirs, partName] = loadSubjectData(subjectName)

% Loads NIRS data for a specific subject from file and ensures proper data 
% structure formatting. Initialises measurement activity list if not 
% present in the original data.

    % Load .nirs data
    partName = subjectName;
    nirs = load(partName, '-mat');
    
    % Check and initialise MeasListAct if it doesn't exist
    if ~isfield(nirs.SD, 'MeasListAct')
        nirs.SD.MeasListAct = nirs.SD.MeasList(:,3);
    end
end