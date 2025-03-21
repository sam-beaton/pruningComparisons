function snrMat = calculateChannelSNR(nirs, chanIdx, windowInfo, snrMat, nsub, resultIdx, params)

%Calculates SNR for a single channel by excluding motion-affected windows, 
% then computing the ratio of signal mean to standard deviation. Returns 
% NaN if insufficient valid samples remain.

    % Extract channel data
    chanData = nirs.d(:, chanIdx);

    if chanIdx > size(windowInfo.windowMotion,1)
        chanIdxMotion = chanIdx - size(windowInfo.windowMotion,1);
    else
        chanIdxMotion = chanIdx;
    end
    
    % Process each window
    for iWindow = 1:windowInfo.windowCount
        % Remove non-valid (motion) samples from SNR calculation
        if windowInfo.windowMotion(chanIdxMotion, iWindow) == 0
            if params.pruneQT == 1
                windowStart = find(nirs.t == windowInfo.windowTimes(1, iWindow));
                windowEnd = find(nirs.t == windowInfo.windowTimes(2, iWindow));
            else
                windowStart = windowInfo.windowTimes(1, iWindow);
                windowEnd = windowInfo.windowTimes(2, iWindow);
            end
            chanData(windowStart:windowEnd) = 0;
        end
    end
    
    % Remove zero elements before calculating mean
    chanData = nonzeros(chanData);
    
    % Debug to stop single/small number of sample values creating 'infinite' means
    if length(find(chanData ~= 0)) < windowInfo.windowSizeSamples
        snrMat(nsub, resultIdx) = NaN;
    else
        % Calculate SNR for remaining samples
        snrMat(nsub, resultIdx) = 20 * (log10((mean(chanData))/(std(chanData))));
    end
end