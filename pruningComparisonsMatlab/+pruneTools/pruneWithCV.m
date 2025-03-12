function [nirs, windowMotion, windowTimes, windowCount, windowSizeSamples] = pruneWithCV(nirs, params, numChan, arrayChans)
    
% Prunes channels using coefficient of variation (CV) thresholds. Analyses
% data in windows, marks motion artifacts, and calculates CV values to 
% determine which channels should be retained for further analysis.

    % Calculate the sampling frequency
    fs = 1 / mean(diff(nirs.t)); 
    
    % Calculate the number of samples per window
    windowSizeSamples = floor(params.windowSec * fs);
    
    % Initialise a matrix for valid samples
    validSamples = true(size(nirs.d));
    
    % Initialise counter for motion detection
    motionCount = zeros(size(nirs.d));

    % Create variables for later SNR calculation
    windowMotion = zeros([size(arrayChans,2) floor(length(nirs.t)./windowSizeSamples)]);
    windowTimes = zeros([2 floor(length(nirs.t)./windowSizeSamples)]);
    windowCount = 0;
    
    % Loop through the data in windows
    for startIdx = 1:windowSizeSamples:(length(nirs.d)-mod(length(nirs.d), windowSizeSamples))
        endIdx = min(startIdx + windowSizeSamples - 1, size(nirs.d, 1));

        windowCount = windowCount+1;
        
        % Store start/endpoints of window to be used later for SNR calc
        windowTimes(1, windowCount) = startIdx;
        windowTimes(2, windowCount) = endIdx;
        
        % Extract the current window of motion data
        motionWindow = nirs.SD.tIncCh(startIdx:endIdx, :);

        for iChan = 1:length(arrayChans)
            if nirs.SD.MeasListAct(iChan) == 1
                % Check for motion (zeros) in the window
                if any(motionWindow(:, iChan) == 0)
                    validSamples(startIdx:endIdx, iChan) = false; % Mark samples as invalid
                    windowMotion(iChan, windowCount) = true; % Mark window as containing motion
                end
            end
        end
    end
    
    % Initialise CV and 'good channel' matrices
    CV = zeros([1 length(nirs.SD.MeasListAct)]);
    idGoodChan = zeros([1 length(nirs.SD.MeasListAct)]);

    % Calculate CV on channel by channel basis, using only motion-free segments
    for iChan = 1:numChan
        CV(iChan) = nanstd(nirs.d(find(validSamples(:,iChan) == 1), iChan)) / nanmean(nirs.d(find(validSamples(:,iChan) == 1), iChan));
    end

    % Check CV against threshold
    for iChan = 1:numChan
        if abs(CV(iChan) <= params.CV && CV(iChan+numChan)) <= params.CV
            idGoodChan(iChan) = 1;
            idGoodChan(iChan+numChan) = 1;
        end
    end
    
    % Check for channels already pruned
    idGoodChan = intersect(find(idGoodChan), find(nirs.SD.MeasListAct));

    % Array for 'good' channels
    goodChans = zeros(length(nirs.SD.MeasListAct), 1); 

    % Set elements of goodChan to 1 according to CV val. above/below thresh
    goodChans(idGoodChan) = 1;

    % Find channels where both chroms are retained
    goodChans = goodChans(1:numChan) + goodChans(numChan+1:end);
    idGoodChan = find(goodChans == 2);

    % Account for multiple chromophores
    idGoodChan = [idGoodChan; idGoodChan + numChan];

    % Set all channels except those in the 'good' list to 0
    nirs.SD.MeasListAct(setdiff(arrayChans,idGoodChan)) = 0;
end