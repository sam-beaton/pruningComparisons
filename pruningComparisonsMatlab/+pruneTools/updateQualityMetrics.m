function [homerPruneMat, sciAvgs, pspAvgs, motionWindows] = updateQualityMetrics(nirs, qualityMatrices, nsub, homerPruneMat, sciAvgs, pspAvgs, motionWindows, saveFlags, params)
    
% Updates the quality metrics matrices with data from the current subject, 
% including SCI averages, PSP averages, and motion window information 
% depending on the pruning method used.

    % Update Homer pruning matrix if needed
    if isfield(saveFlags, 'saveHomerPrune')
        homerPruneMat(:, nsub) = nirs.SD.MeasListAct;
    end
    
    % Update quality metrics for QT-NIRS method
    if params.pruneQT == 1
        % Update SCI averages if needed
        if isfield(saveFlags, 'saveAvgSCI')
            sciAvgs(:, nsub) = qualityMatrices.avg_sci_channel;
        end
        
        % Update PSP averages if needed
        if isfield(saveFlags, 'saveAvgPSP')
            pspAvgs(:, nsub) = qualityMatrices.avg_psp_channel;
        end
        
        % Update motion windows if needed
        if isfield(saveFlags, 'saveMotion')
            motionWindows(:, nsub) = sum(qualityMatrices.channelWindowMarked, 2) ./ size(qualityMatrices.channelWindowMarked, 2);
        end
    end
end