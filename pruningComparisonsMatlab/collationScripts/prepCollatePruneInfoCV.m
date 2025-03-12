clear; close;

%specify ages to include in analysis
ages = {'05', '08', '12', '18', '24'};

%specify other info about data
tasks = {'hand', 'social'};
cohorts = {'gm', 'uk'};

% create empty table
columnNames = {'ID', 'Cohort', 'Task', 'Age', 'Channels', 'ROI_SNR'};
variableTypes = {'string', 'string', 'string', 'double', 'double', 'double'};
data = table('Size', [0, length(columnNames)], ...
                   'VariableTypes', variableTypes, ...
                   'VariableNames', columnNames);

for iCohort = 1:length(cohorts)

    cohort = cohorts{iCohort};
    statsDir = strcat('[path-to...]/pruningComparisonsMatlab/stats/bright/', cohort, '/');
    
    for iAge = 1:length(ages)
    
        age = ages{iAge};
    
        for iTask = 1:length(tasks)

            task = tasks{iTask};
            taskDir = strcat(age, 'mo/', task, '/prune/');
        
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

            %create variable for number of parts to prevent having to 
            % recalculate every time
            numParts = length(partIDs);
        
%             % load tables containing channel characteristics
%             % previously pruned channels using Homer
%             load(fullfile([statsDir, taskDir, strcat(task, age, 'moHomerPrunedChannels.mat')]));
%             % avg SCI
%             load(fullfile([statsDir, taskDir, strcat(task, age, 'moAvgSCImat.mat')]));            
%             % avgPSP
%             load(fullfile([statsDir, taskDir, strcat(task, age, 'moAvgPSPmat.mat')]));
%             % percentage of motion Windows
%             load(fullfile([statsDir, taskDir, strcat(task, age, 'moMotionWindowsMat.mat')]));

            % load SNR matrix
            load(strcat(task,'FullPruneChan_CV_20.mat'));
            %load Channels retained matrix
            load(strcat(task, 'FullPruneSNR_CV_20.mat'));

            %convert age to numeric
            ageDouble = str2double(age);
            
            % create table for this age group
            ID = partIDs;
            Cohort = repmat(cohort, [length(partIDs) 1]);
            Task = repmat(task, [length(partIDs) 1]);
            Age = repmat(ageDouble,[length(partIDs) 1]);
            Channels = zeros([length(partIDs) 1]);
            ROI_SNR = zeros([length(partIDs) 1]);
            ageTaskTable = table(ID, Cohort, Task, Age, Channels, ROI_SNR);
            
            % counter for rows
            iRow = 1;
            
            for j = 1:numParts
                
                % assign channels
                ageTaskTable.Channels(iRow) = channelsRetained(j);
                
                %assign avg. ROI SNR
                ageTaskTable.ROI_SNR(iRow) = nanmean(snrMat(j,:));

                iRow = iRow + 1;
            end
        
            % Add ageTable to dataTable
            data = [data; ageTaskTable];
        
        end
    
    end

end

% save table
writetable(data, fullfile(['[path-to...]/pruningComparisonsMatlab/stats/bright/overall/', 'pruneChannelInfoTableCV.csv']));
