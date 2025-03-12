clear; close;

%specify ages to include in analysis
ages = {'05', '08', '12', '18', '24'};

%specify other info about data
tasks = {'hand', 'social'};
cohorts = {'gm', 'uk'};

% create empty table
columnNames = {'ID', 'Cohort', 'Task', 'Age', 'Channel', 'Homer_Pruned_Total', 'Average_SCI', 'Average_PSP', 'Percentage_Motion_Windows', 'AP_Orientation'};
variableTypes = {'string', 'string', 'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double'};
data = table('Size', [0, length(columnNames)], ...
                   'VariableTypes', variableTypes, ...
                   'VariableNames', columnNames);

for iCohort = 1:length(cohorts)

    cohort = cohorts{iCohort};
    statsDir = strcat('[path-to...]/pruningComparisonsMatlab/stats/bright/', cohort, '/overall/');
    
    for iTask = 1:length(tasks)
    
        task = tasks{iTask};
        taskDir = strcat(task, '/prune/');
    
        for iAge = 1:length(ages)
        
            age = ages{iAge};
            namesDir = strcat('[path-to...]/bright/', cohort, '/', age, 'mo/', task);
            cd(strcat(statsDir, taskDir));
            
            % get participant names
            if iCohort == 1 %hand
                searchNames = fullfile(namesDir, '*.nirs');
                namesSub = dir(searchNames);
                partIDs = {namesSub.name};
                partIDs = cellfun(@(x) x(1:min(4, end)), partIDs, 'UniformOutput', false);
                partIDs = partIDs(:);  % Ensures partIDs is a column cell array
            elseif iCohort ==2 %social
                searchNames = fullfile(namesDir, '*.nirs');
                namesSub = dir(searchNames);
                partIDs = {namesSub.name};
                partIDs = cellfun(@(x) x(5:min(8, end)), partIDs, 'UniformOutput', false);
                partIDs = partIDs(:);  % Ensures partIDs is a column cell array
            end

            % get total number of channels in array for ROI/Excluded channels
            % calculation
            load(fullfile([namesDir, '/', namesSub(1).name]), '-mat', 'SD');
            numChans = length(SD.MeasListAct);

            %create variable containing channel numbers
            chanNums = cell(numChans, 1);
            for iChan = 1:numChans
                chanNums{iChan} = iChan;
            end

            %create variable for number of parts to prevent having to 
            % recalculate every time
            numParts = length(partIDs);
        
            % load tables containing channel characteristics
            % previously pruned channels using Homer
            load(fullfile([statsDir, taskDir, strcat(task, age, 'moHomerPrunedChannels.mat')]));
            % avg SCI
            load(fullfile([statsDir, taskDir, strcat(task, age, 'moAvgSCImat.mat')]));            
            % avgPSP
            load(fullfile([statsDir, taskDir, strcat(task, age, 'moAvgPSPmat.mat')]));
            % percentage of motion Windows
            load(fullfile([statsDir, taskDir, strcat(task, age, 'moMotionWindowsMat.mat')]));

            %convert age to numeric
            age = str2double(age);
            
            % create table for this age group
            ID = repmat(partIDs, numChans, 1);
            Cohort = repmat(cohort, [numChans*numParts 1]);
            Task = repmat(task, [numChans*numParts 1]);
            Age = repmat(age,[numChans*numParts 1]);
            Channel = zeros([numChans*numParts 1]);
            Homer_Pruned_Total = zeros([numChans*numParts 1]);
            Average_SCI = zeros([numChans*numParts 1]);
            Average_PSP = zeros([numChans*numParts 1]);
            Percentage_Motion_Windows = zeros([numChans*numParts 1]);
            AP_Orientation = zeros([numChans*numParts 1]);
            ageTaskTable = table(ID, Cohort, Task, Age, Channel, Homer_Pruned_Total, Average_SCI, Average_PSP, Percentage_Motion_Windows, AP_Orientation);
            
            % counter for rows
            iRow = 1;
            
            for i=1:numChans
            
                for j = 1:numParts
                    
                    % assign channel number
                    ageTaskTable.Channel(iRow) = i;
                    
                    % generate a second i value using modulo to account for
                    % sciAvgs and pspAvgs only having 1 value per channel
                    iMod = mod(i - 1, size(sciAvgs, 1)) + 1;
            
                    ageTaskTable.Homer_Pruned_Total(iRow) = numChans - sum(homerPruneMat(:,j));
                    ageTaskTable.Average_SCI(iRow) = mean(sciAvgs(iMod, j), 'omitnan'); 
                    ageTaskTable.Average_PSP(iRow) = mean(pspAvgs(iMod, j), 'omitnan'); 
                    ageTaskTable.Percentage_Motion_Windows(iRow) = motionWindows(i, j);

                    switch i
                        case {2, 7}
                            ageTaskTable.AP_Orientation(iRow) = 1;
                        case {1, 5, 6, 10}
                            ageTaskTable.AP_Orientation(iRow) = 2;
                        case {4, 9}
                            ageTaskTable.AP_Orientation(iRow) = 3;
                        case {3, 8, 13, 19}
                            ageTaskTable.AP_Orientation(iRow) = 4;
                        case {12, 18}
                            ageTaskTable.AP_Orientation(iRow) = 5;
                        case {11, 16, 17, 22}
                            ageTaskTable.AP_Orientation(iRow) = 6;
                        case {15, 21}
                            ageTaskTable.AP_Orientation(iRow) = 7;
                        case {14, 20, 25, 28}
                            ageTaskTable.AP_Orientation(iRow) = 8;
                        case {24, 27}
                            ageTaskTable.AP_Orientation(iRow) = 9;
                        case {23, 26, 31, 34}
                            ageTaskTable.AP_Orientation(iRow) = 10;
                        case {30, 33}
                            ageTaskTable.AP_Orientation(iRow) = 11;
                        case {29, 32}
                    end

                    iRow = iRow + 1;
                end
            end
        
            % Add ageTable to dataTable
            data = [data; ageTaskTable];
        
        end
    
    end

end

% save table
writetable(data, fullfile(['[path-to...]/pruningComparisonsMatlab/stats/bright/overall/', 'pruneChannelInfoTableQT.csv']));
