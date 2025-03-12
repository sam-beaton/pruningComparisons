function [nirs, qualityMatrices] = pruneTools.pruneWithQT(nirs, params, numChan)

% Applies the QT-NIRS algorithm to prune channels based on scalp coupling 
% index (SCI) and power spectral power (PSP) thresholds. Ensures that both 
% chromophore channels are consistently pruned if one is marked for 
% removal.

    % Ensure both chroms pruned if channel removed by Homer enPrune
    homChans = find(nirs.SD.MeasListAct == 0);
    hboChans = homChans(homChans <= numChan) + numChan;        
    hbrChans = homChans(homChans > numChan) - numChan;
    homChans = [homChans, hboChans, hbrChans];
    homChans = unique(homChans);
    
    % Prune corresponding chrom channels using list 'homChans'
    nirs.SD.MeasListAct(homChans) = 0;

    % Prune channels with QT NIRS
    qualityMatrices = sbPrePqtnirs( ...
        nirs, ...
        'freqCut',[params.bpFmin, params.bpFmax], ...
        'window', params.windowSec, ...
        'overlap', params.windowOverlap, ...
        'qualityThreshold', params.quality_threshold, ...
        'sciThreshold', params.sci_threshold, ...
        'pspThreshold', params.psp_threshold, ...
        'conditionsMask','all', ...
        'dodFlag', 0, ...
        'guiFlag', params.gui_flag);

    % Insert correct values for pruning 
    % qtnirs changes the channel order during computation so need to find
    % correct rows which have been pruned first
    for i = 1:length(nirs.SD.MeasListAct)
        % Only changes if qualityMatrices val is 0 so as to not override
        % previous changes by enPruneChannels
        if qualityMatrices.MeasListAct(i) == 0 
            iRow = find(nirs.SD.MeasList(:,1) == qualityMatrices.MeasList(i, 1) ...
                        & nirs.SD.MeasList(:,2) == qualityMatrices.MeasList(i, 2) ...
                        & nirs.SD.MeasList(:,4) == qualityMatrices.MeasList(i, 4));
            nirs.SD.MeasListAct(iRow) = 0;
        end
    end
end