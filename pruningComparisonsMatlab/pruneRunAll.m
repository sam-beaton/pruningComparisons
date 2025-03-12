%% pruneRunAll
% Script which calls pruneGetVals for specifed pruning method (saves
% having to repeat in every iteration of the preprocessing pipeline)
% whilst saving the .nirs files
%
% SLB 15/2/24


clear; close all;

%% Add relevant toolboxes to current path
% Underneath parameter settings for ease of access to change
addpath(genpath('[path-to...]/Homer2'))
addpath(genpath('[path-to...]/qt-nirs'))
addpath(genpath('[path-to...]/pruningComparisons')) [path-to...]

%initialise parameters variable
params = struct();

%Change cohort variables here (won't change for entire run of script)
params.dataRawLoc = '...'; % Main directory with original .nirs data file 
params.dataOutLoc = '...'; %parent directory for saved pruned files
params.statsOutLoc = 'e.g. [path-to...]/stats/'; %parent directory for saved data files

%%% ============ Cohort variables/arguments =================
cohorts = {'gm', 'uk'}; % 'uk' or 'gm'
tasks = {'hand', 'social'}; % e.g. 'hand', 'social'
timepoints = {'05mo', '08mo', '12mo', '18mo', '24mo'}; %'01mo', '05mo', '08mo', '12mo', ..., '60mo'

%%% ======= Pruning variables/arguments ============
params.pruneQT = 0; %1 = QT-NIRS; 0 = CV pruning; 2 = QT-NIRS no Homer pre-pruning; 3 = just Homer pre-pruning
detailLevel = 'fine'; % broad or fine
if pruneQT == 1
    sciThresholdValues = [0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05]; %your param values here
    detailLevel = lower(detailLevel); %just in case!
    if strcmp(detailLevel, 'broad')
        pspThresholdValues = [0.1, 0.075, 0.05, 0.025]; %broader - as a first pass; can be changed
    elseif strcmp(detailLevel, 'fine')
        pspThresholdValues = [0.1, 0.095, 0.09, 0.085, 0.08, 0.075, 0.07, 0.065, 0.06, 0.055, 0.05, 0.045, 0.04, 0.035, 0.03, 0.025, 0.02, 0.015, 0.01, 0.005]; %finer (based on broad first pass, change as necessary)
    else
        error("Please enter a valid detail level - broad or fine");
    end
    detailLevel(1) = upper(detailLevel(1)); %for filename later in script
end

for iCohort = 1:length(cohorts)

    params.cohort = cohorts{iCohort};

    for iTask = 1:length(tasks)

        params.task = tasks{iTask};

        if pruneQT == 1
            % Start a parallel pool (specify number of workers)
            parpool(6);
            
            % Initialize storage for pruneMetrics outside the loop
            pruneMetricsSummary = cell(length(timepoints), 1); % Cell array to store results for each timepoint
            
            parfor i = 1:length(timepoints)

                params.timepoint = timepoints{i};
                pruneMatrix = zeros([2, length(pspThresholdValues) * length(sciThresholdValues)]); %initialise
                colList = 'Columns = SCI&PSP parameters: '; %initialise
                
                %%% ======= Run pruning ==========
                for sciThreshold = sciThresholdValues
                    for pspThreshold = pspThresholdValues

                        params.sciThreshold = sciThreshold
                        params.pspThreshold = pspThreshold
                        
                        %prune:        
                        [meanChans, meanSNR] = pruneGetVals(params);                
                        
                        %update summary data:
                        nextIndices = find(pruneMatrix == 0);
                        pruneMatrix(nextIndices(1)) = meanChans;
                        pruneMatrix(nextIndices(2)) = meanSNR;
                        
                        %add parameter values to list for metadata (for later use)
                        paramsComb = strjoin({num2str(sciThreshold), num2str(pspThreshold)}, '&');
                        colList = strjoin({colList, paramsComb}, '  ');
                    end
                end
                
                % Store the result for this timepoint
                pruneMetrics = struct();
                pruneMetrics.metaData = strjoin({'Rows = (All are mean values, calculated across all participants) (1) Mean of Mean# Channels included; (2) Mean of Mean SNR', colList}, ' || ');
                pruneMetrics.sciThresholds = sciThresholdValues;
                pruneMetrics.pspThresholds = pspThresholdValues;
                pruneMetrics.summaryMatrix = real(pruneMatrix); % just in case! sometimes produces complex numbers
                pruneMetricsSummary{i} = pruneMetrics; % Save this timepoint's results
            end
        
            % save each result after the parfor loop
            for i = 1:length(timepoints)
                pruneMetrics = pruneMetricsSummary{i};
                timepoint = timepoints{i};
                statsDirecName = strcat(statsOutLoc, cohort, '/overall/', task, '/prune/');
                if ~exist(statsDirecName, 'dir')
                    mkdir(statsDirecName)
                end
                pruneMetricsFileName = strcat(statsDirecName, task, timepoint, 'SummaryPruneMetrics.mat');
                save(pruneMetricsFileName, 'pruneMetrics');
            end
        
            % Shut down the parallel pool
            delete(gcp('nocreate'));
        
        else
        
            for i = 1:length(timepoints)
        
                %%% ======= Initialise variables ==========
                timepoint = timepoints{i};
                %just so function runs
                sciThreshold = 0.8;
                pspThreshold = 0.1;
        
                %prune:        
                [meanChans, meanSNR] = pruneGetVals(params);
        
            end
        end
    end
end

fprintf("COMPLETE \n")