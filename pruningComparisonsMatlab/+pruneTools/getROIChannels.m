function roiChans = getROIChannels(task, timepoint)

% Selects region of interest (ROI) channels specific to the experimental 
% task and timepoint. Returns an array of channel indices that correspond 
% to the relevant brain regions for the analysis.

    % Select ROI channels based on task and timepoint
    switch task
        case 'hand'
            switch timepoint
                case '05mo'
                    roiChans = [3, 12, 17, 8, 18, 9, 19];
                case '08mo'
                    roiChans = [3, 12, 4, 13, 18];
                case '12mo'
                    roiChans = [3, 12, 4, 21, 8, 18, 9];
                case '18mo'
                    roiChans = [3, 12, 16, 4, 8, 18, 22, 9, 19];
                case '24mo'
                    roiChans = [15, 3, 12, 16, 4, 8, 18, 22];
                case '60mo'
                    roiChans = [3, 12, 4, 13, 17, 8, 18, 22, 19];
            end
        case 'social'
            switch timepoint
                case '05mo'
                    roiChans = [3, 5, 9, 10, 12, 18, 24, 27];
                case '08mo'
                    roiChans = [2, 7, 10, 27];
                case '12mo'
                    roiChans = [3, 4, 11, 12, 13, 17, 18, 19, 20, 21];
                case '18mo'
                    roiChans = [3, 10, 12, 18, 27];
                case '24mo'
                    roiChans = [8, 14, 15, 18, 27];
                case '60mo'
                    roiChans = [2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 17, 18, 19, 20, 22, 24, 26, 27, 31];
            end
    end
end