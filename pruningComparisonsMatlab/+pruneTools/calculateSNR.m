function snrMat = calculateSNR(nirs, roiChans, windowInfo, snrMat, nsub, params)

% Calculates signal-to-noise ratio (SNR) for all ROI channels (both HbO 
% and HbR) for the current subject, excluding windows with detected motion 
% artifacts.

    for iChan = 1:length(roiChans)
        % Calculate SNR for HbO channels
        if nirs.SD.MeasListAct(roiChans(iChan)) == 1
            snrMat = pruneTools.calculateChannelSNR(nirs, roiChans(iChan), windowInfo, snrMat, nsub, iChan, params);
        end
        
        % Calculate SNR for HbR channels
        hbrChanIdx = roiChans(iChan) + (length(nirs.SD.MeasListAct)/2);
        if nirs.SD.MeasListAct(hbrChanIdx) == 1
            snrMat = pruneTools.calculateChannelSNR(nirs, hbrChanIdx, windowInfo, snrMat, nsub, iChan+length(roiChans), params);
        end
    end
end