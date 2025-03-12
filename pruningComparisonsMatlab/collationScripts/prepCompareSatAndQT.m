
clear; close all;

homerOnlyDir = '[path-to...]/pruningComparisonsMatlab/stats/brightHomerOnly';
qtOnlyDir = '[path-to...]/pruningComparisonsMatlab/stats/brightAllParams';


%specify ages to include in analysis
ages = {'05', '08', '12', '18', '24'};

%specify other info about data
tasks = {'hand', 'social'};
cohorts = {'gm', 'uk'};

%% ---------- collate data -----------------

numGroupCombinations = length(ages)*length(tasks)*length(cohorts);
qtCurrentDir = strcat(qtOnlyDir, '/', cohorts{1}, '/', ages{1}, 'mo/', tasks{1}, '/prune');
numParamCombinations = length(dir(fullfile(qtCurrentDir, '*SNR*.mat')));

% initilaise variables for saturation table
satCohort = cell(numGroupCombinations, 1);
satTask = cell(numGroupCombinations, 1);
satAge = cell(numGroupCombinations, 1);
satROI_SNR = zeros([numGroupCombinations 1]);

qtCohort = cell(numGroupCombinations*numParamCombinations, 1);
qtTask = cell(numGroupCombinations*numParamCombinations, 1);
qtAge = cell(numGroupCombinations*numParamCombinations, 1);
qtSCIParam = zeros([numGroupCombinations*numParamCombinations 1]);
qtPSPParam = zeros([numGroupCombinations*numParamCombinations 1]);
qtROI_SNR = zeros([numGroupCombinations*numParamCombinations 1]);

satCount = 1;
qtCount = 1;

for iCohort = 1:length(cohorts)

    cohort = cohorts{iCohort};
    
    for iTask = 1:length(tasks)
    
        task = tasks{iTask};
    
        for iAge = 1:length(ages)

            age = ages{iAge};
            
            %saturated table
            satCohort{satCount} = cohort;
            satTask{satCount} = task;
            satAge{satCount} = age;

            load(strcat(homerOnlyDir, '/', cohort, '/', age, 'mo/', task, '/prune/', task, 'FullPruneSNR_Homer.mat'));
            snrMat = real(snrMat);
            meanSatSNR = mean(snrMat, 2, 'omitnan');
            meanSatSNR = mean(meanSatSNR, 'omitnan');
            satROI_SNR(satCount) = meanSatSNR;

            % qt table
            qtCurrentDir = strcat(qtOnlyDir, '/', cohort, '/', age, 'mo/', task, '/prune');
            mats = dir(fullfile(qtCurrentDir, '*SNR*.mat'));

            for iMat = 1:numParamCombinations

                qtCohort{qtCount} = cohort;
                qtTask{qtCount} = task;
                qtAge{qtCount} = age;

                currMat = mats(iMat);

                % get SCI parameter
                startIndex = strfind(currMat.name, 'SCI');
                startIndex = startIndex + length('SCI');
                endIndex = strfind(currMat.name, '_PSP');
                endIndex = endIndex - 1;
                sciParam = currMat.name(startIndex:endIndex);
                dotPos = 2;
                part1 = sciParam(1:dotPos-1);
                part2 = sciParam(dotPos:end);
                sciParam = [part1 '.' part2];
                qtSCIParam(qtCount) = str2double(sciParam);
                
                %get psp parameter
                startIndex = strfind(currMat.name, 'PSP');
                startIndex = startIndex + length('PSP');
                endIndex = strfind(currMat.name, '.mat');
                endIndex = endIndex - 1;
                pspParam = currMat.name(startIndex:endIndex);
                dotPos = 2;
                part1 = pspParam(1:dotPos-1);
                part2 = pspParam(dotPos:end);
                pspParam = [part1 '.' part2];
                qtPSPParam(qtCount) = str2double(pspParam);
    
                load(fullfile(qtCurrentDir, currMat.name));
                snrMat = real(snrMat);
                meanSatSNR = mean(snrMat, 2, 'omitnan');
                qtROI_SNR(qtCount) = mean(meanSatSNR, 'omitnan');

                qtCount = qtCount+1;
            end

            satCount = satCount +1;
           
        end
    end
end

satTable = table(satCohort, satTask, satAge, satROI_SNR, ...
    'VariableNames', {'Cohort', 'Task', 'Age', 'ROI_SNR'});

qtTable = table(qtCohort, qtTask, qtAge, qtSCIParam, qtPSPParam, qtROI_SNR, ...
    'VariableNames', {'Cohort', 'Task', 'Age', 'SCIParam', 'PSPParam', 'ROI_SNR'});


%% ---- generate heatmaps ---------------

% Get unique combinations of Cohort, Age, and Task in qtTable
uniqueCombinations = unique(qtTable(:, {'Cohort', 'Age', 'Task'}));

% Loop over each unique combination
for i = 1:height(uniqueCombinations)
    
    % Extract current combination details
    cohort = uniqueCombinations.Cohort{i};
    age = uniqueCombinations.Age{i};
    task = uniqueCombinations.Task{i};
    
    % Filter data for the current combination in qtTable
    subsetData = qtTable(strcmp(qtTable.Cohort, cohort) & ...
                         strcmp(qtTable.Age, age) & ...
                         strcmp(qtTable.Task, task), :);
                     
    % Get the corresponding threshold value from satTable
    satThreshold = satTable.ROI_SNR(strcmp(satTable.Cohort, cohort) & ...
                                    strcmp(satTable.Age, age) & ...
                                    strcmp(satTable.Task, task));
    
    % Prepare data for heatmap
    [x, ~, xIdx] = unique(subsetData.PSPParam);
    [y, ~, yIdx] = unique(subsetData.SCIParam);
    heatmapData = accumarray([yIdx, xIdx], subsetData.ROI_SNR, [], @mean, NaN);
    
    % Apply greying out: set values less than satThreshold to NaN
    heatmapData(heatmapData < satThreshold) = NaN;
    
    % Generate the heatmap
    figure;
    h = heatmap(x, y, heatmapData, 'Colormap', jet, 'ColorbarVisible', 'on');
    h.MissingDataLabel = 'Greyed out';
    h.MissingDataColor = [0.8, 0.8, 0.8]; % light grey for values below threshold
    xlabel('PSPParam');
    ylabel('SCIParam');
    % Set the title with satThreshold included
    title(sprintf('Cohort: %s, Age: %s mo, Task: %s | satTable ROI\\_SNR: %.2f', cohort, age, task, satThreshold));    
end


% ---- generate heatmaps with 4 per figure ---------------

% Get unique combinations of Cohort, Age, and Task in qtTable
uniqueCombinations = unique(qtTable(:, {'Cohort', 'Age', 'Task'}));

% Number of heatmaps per figure
heatmapsPerFigure = 4;

% Total number of combinations
numCombinations = height(uniqueCombinations);

% Number of figures required
numFigures = ceil(numCombinations / heatmapsPerFigure);

% Loop over each unique combination
for figIdx = 1:numFigures
    
    % Create a new figure
    figure;
    
    % Determine the heatmaps to display in this figure
    startIdx = (figIdx - 1) * heatmapsPerFigure + 1;
    endIdx = min(figIdx * heatmapsPerFigure, numCombinations);
    numHeatmapsInFigure = endIdx - startIdx + 1;
    
    % Loop over heatmaps for this figure
    for subplotIdx = 1:numHeatmapsInFigure
        
        % Current combination index
        comboIdx = startIdx + subplotIdx - 1;
        
        % Extract current combination details
        cohort = uniqueCombinations.Cohort{comboIdx};
        age = uniqueCombinations.Age{comboIdx};
        task = uniqueCombinations.Task{comboIdx};
        
        % Filter data for the current combination in qtTable
        subsetData = qtTable(strcmp(qtTable.Cohort, cohort) & ...
                             strcmp(qtTable.Age, age) & ...
                             strcmp(qtTable.Task, task), :);
                         
        % Get the corresponding threshold value from satTable
        satThreshold = satTable.ROI_SNR(strcmp(satTable.Cohort, cohort) & ...
                                        strcmp(satTable.Age, age) & ...
                                        strcmp(satTable.Task, task));
        
        % Prepare data for heatmap
        [x, ~, xIdx] = unique(subsetData.PSPParam);
        [y, ~, yIdx] = unique(subsetData.SCIParam);
        heatmapData = accumarray([yIdx, xIdx], subsetData.ROI_SNR, [], @mean, NaN);
        
        % Apply greying out: set values less than satThreshold to NaN
        heatmapData(heatmapData < satThreshold) = NaN;
        
        % Generate the heatmap in a subplot
        subplot(2, 2, subplotIdx); % 4 heatmaps per figure (2x2 grid)
        h = heatmap(x, y, heatmapData, 'Colormap', jet, 'ColorbarVisible', 'on');
        h.MissingDataLabel = 'Greyed out';
        h.MissingDataColor = [0.8, 0.8, 0.8]; % light grey for values below threshold
        xlabel('PSPParam');
        ylabel('SCIParam');
        % Set the title with satThreshold included
        title(sprintf('Cohort: %s, Age: %s mo, Task: %s | ROI\\_SNR: %.2f', cohort, age, task, satThreshold));
    end
end


