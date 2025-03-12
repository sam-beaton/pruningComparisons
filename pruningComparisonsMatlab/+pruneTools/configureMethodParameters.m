function params = configureMethodParameters(params, timepoint)

% Sets up the specific parameters for signal processing and artifact 
% detection methods. Includes thresholds for motion detection, channel 
% pruning, and signal quality assessment based on established research 
% parameters.

    % Motion detection windows threshold
    params.badWindowThresh = 3;
    
    % Channel pruning parameters
    params.dRange = [3e-04 1e+07]; % LLoyd-Fox et al. 2019 & Di Lorenzo et al. 2019
    params.SNRthresh = 0; % Di Lorenzo et al. 2019
    params.SDrange = [0 45]; % Di Lorenzo et al. 2019
    
    % Threshold of Coefficient of Variation
    params.CV = 0.08; % Lloyd-Fox et al. 2019
    
    % QT-NIRS parameters
    params.bpFmin = 1.2; % Minagawa et al. 2023
    params.bpFmax = 3.0; % Minagawa et al. 2023
    params.windowSec = 3; % Better motion detection exclusion using windows from QT-NIRS pruning
    params.windowOverlap = 0; % Pollonini et al 2016
    params.quality_threshold = 0.75; % QT-NIRS quality threshold
    params.gui_flag = 0; % Set to 1 to visualize QT-NIRS quality graphics
    
    % Motion detection parameters
    params.tMotion = 1; % Di Lorenzo et al. 2019
    params.tMask = 0; % Only consider the motion itself for QT-NIRS
    params.STDEVthresh = 15; % Di Lorenzo et al. 2019
    params.AMPthresh = 0.4; % Di Lorenzo et al. 2019
end