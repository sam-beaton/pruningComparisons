clear; close;

%specify ages to include in analysis
ages = {'05', '08', '12', '18', '24'};

%specify other info about data
tasks = {'hand', 'social'};
cohorts = {'gm', 'uk'};

statsDir = '[path-to...]/pruningComparisonsMatlab/stats/bright/'; %CHANGE

% create empty table
columnNames = {'ID', 'Cohort', 'Task', 'Age', 'SCI', 'PSP', 'ROI_SNR', 'Channels_Retained', 'Homer_Pruned_Total', 'Average_SCI', 'Average_PSP', 'Percentage_Motion_Windows', 'SNR_Mult_CR', 'Infant_Excluded'};
variableTypes = {'string', 'string', 'string', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'};
data = table('Size', [0, length(columnNames)], ...
                   'VariableTypes', variableTypes, ...
                   'VariableNames', columnNames);

for iCohort = 1:length(cohorts)

    cohort = cohorts{iCohort};
    cohortStatsDir = strcat('/Users/sambe/Documents/MATLAB/matlabProjects/PreP/stats/bright/', cohort, '/');
    
    for iTask = 1:length(tasks)
    
        task = tasks{iTask};
    
        for iAge = 1:length(ages)
        
            age = ages{iAge};
            ageStatsDir = strcat(cohortStatsDir, age, 'mo/', task, '/prune');
            namesDir = strcat('[path-to...]/bright/', cohort, '/', age, 'mo/', task); %CHANGE
            
            cd(ageStatsDir)
            chanSub = dir(strcat(task, 'FullPruneChan_P*'));
            snrSub = dir(strcat(task, 'FullPruneSNR_P*'));
            
            % get participant names
            if iCohort == 1
                searchNames = fullfile(namesDir, '*.nirs');
                namesSub = dir(searchNames);
                partIDs = {namesSub.name};
                partIDs = cellfun(@(x) x(1:min(4, end)), partIDs, 'UniformOutput', false);
                partIDs = partIDs(:);  % Ensures partIDs is a column cell array
            elseif iCohort ==2
                searchNames = fullfile(namesDir, '*.nirs');
                namesSub = dir(searchNames);
                partIDs = {namesSub.name};
                partIDs = cellfun(@(x) x(5:min(8, end)), partIDs, 'UniformOutput', false);
                partIDs = partIDs(:);  % Ensures partIDs is a column cell array
            end

            %create variable for length of participants to prevent having to recalculate
            %every time
            numParts = length(partIDs);
        
            % load tables containing channel characteristics
            % previously pruned channels using Homer
            load(fullfile([cohortStatsDir, 'overall/', task, '/prune/', strcat(task, age, 'moHomerPrunedChannels.mat')]));
            % avg SCI
            load(fullfile([cohortStatsDir, 'overall/', task, '/prune/', strcat(task, age, 'moAvgSCImat.mat')]));            
            % avgPSP
            load(fullfile([cohortStatsDir, 'overall/', task, '/prune/', strcat(task, age, 'moAvgPSPmat.mat')]));
            % percentage of motion Windows
            load(fullfile([cohortStatsDir, 'overall/', task, '/prune/', strcat(task, age, 'moMotionWindowsMat.mat')]));

            %convert age to numeric
            age = str2double(age);
            
            % get total number of channels in array for ROI/Excluded channels
            % calculation
            load(fullfile([namesDir, '/', namesSub(1).name]), '-mat', 'SD');
            numChans = length(SD.MeasListAct);
            
            % create table for this age group
            ID = repmat(partIDs, length(chanSub), 1);
            Cohort = repmat(cohort, [length(partIDs)*length(chanSub) 1]);
            Task = repmat(task, [length(partIDs)*length(chanSub) 1]);
            SCI = zeros([length(partIDs)*length(chanSub) 1]);
            PSP = zeros([length(partIDs)*length(chanSub) 1]);
            Age = zeros([length(partIDs)*length(chanSub) 1]);
            ROI_SNR = zeros([length(partIDs)*length(chanSub) 1]);
            Channels_Retained = zeros([length(partIDs)*length(chanSub) 1]);
            Homer_Pruned_Total = zeros([length(partIDs)*length(chanSub) 1]);
            Average_SCI = zeros([length(partIDs)*length(chanSub) 1]);
            Average_PSP = zeros([length(partIDs)*length(chanSub) 1]);
            Percentage_Motion_Windows = zeros([length(partIDs)*length(chanSub) 1]);
            SNR_Mult_CR= zeros([length(partIDs)*length(chanSub) 1]);
            Infant_Excluded = zeros([length(partIDs)*length(chanSub) 1]);
            ageTaskTable = table(ID, Cohort, Age, Task, SCI, PSP, ROI_SNR, Channels_Retained, Homer_Pruned_Total, Average_SCI, Average_PSP, Percentage_Motion_Windows, SNR_Mult_CR, Infant_Excluded);
            
            % counter for rows
            iRow = 1;
            
            for i=1:length(chanSub)
            
                % get SCI parameter
                startIndex = strfind(chanSub(i).name, 'SCI');
                startIndex = startIndex + length('SCI');
                endIndex = strfind(chanSub(i).name, '_PSP');
                endIndex = endIndex - 1;
                sciParam = chanSub(i).name(startIndex:endIndex);
                dotPos = 2;
                part1 = sciParam(1:dotPos-1);
                part2 = sciParam(dotPos:end);
                sciParam = [part1 '.' part2];
                sciParam = str2double(sciParam);
                
                %get psp parameter
                startIndex = strfind(chanSub(i).name, 'PSP');
                startIndex = startIndex + length('PSP');
                endIndex = strfind(chanSub(i).name, '.mat');
                endIndex = endIndex - 1;
                pspParam = chanSub(i).name(startIndex:endIndex);
                dotPos = 2;
                part1 = pspParam(1:dotPos-1);
                part2 = pspParam(dotPos:end);
                pspParam = [part1 '.' part2];
                pspParam = str2double(pspParam);
                
                % load chan and snr data
                load(snrSub(i).name);
                load(chanSub(i).name);
            
                for j = 1:numParts
            
                    ageTaskTable.Age(iRow) = age;
                    ageTaskTable.SCI(iRow) = sciParam;
                    ageTaskTable.PSP(iRow) = pspParam;
                    ageTaskTable.ROI_SNR(iRow) = mean(snrMat(j, :), 'omitnan');
                    ageTaskTable.Channels_Retained(iRow) = channelsRetained(j);
                    ageTaskTable.Homer_Pruned_Total(iRow) = numChans - sum(homerPruneMat(:,j));
                    ageTaskTable.Average_SCI(iRow) = mean(sciAvgs(:, j), 'omitnan'); %sciAvgs  mean(sciAvgs(:, j), 'omitnan')
                    ageTaskTable.Average_PSP(iRow) = mean(pspAvgs(:, j), 'omitnan'); %pspAvgs
                    ageTaskTable.Percentage_Motion_Windows(iRow) = mean(motionWindows(:, j), 'omitnan'); %motionWindows
                    ageTaskTable.SNR_Mult_CR(iRow) = (mean(snrMat(j, :), 'omitnan')) * channelsRetained(j);
                    if channelsRetained(j) < 0.6*(numChans)
                        ageTaskTable.Infant_Excluded(iRow) = 1;
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
writetable(data, fullfile([statsDir, 'overall/pruneMLMInputTable.csv']));

% Get info about infants with all channels pruned and those excluded due to
% < 60% of channels 
% ========= No channels remaining ==========================

noChans = find(data.Channels_Retained ==0 & data.Task == "hand");
dataNoChans = data(noChans, :);


dataNoChans = dataNoChans{:, {'ID', 'Age', 'SCI', 'PSP'}};

uniqueNoChans = str2double(unique(dataNoChans, 'rows'));
uniqueNoChans = sortrows(uniqueNoChans, [2 1]);

% ========= Excluded ==========================

exclInfants = find(data.Channels_Retained < (0.6*68) & data.Task == "hand");
dataExcl = data(exclInfants, :);


dataExcl = dataExcl{:, {'ID', 'Age', 'SCI', 'PSP'}};

uniqueExcl = str2double(unique(dataExcl, 'rows'));
uniqueExcl = sortrows(uniqueExcl, [2 1]);