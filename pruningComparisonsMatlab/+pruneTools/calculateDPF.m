function dpf = calculateDPF(timepoint, nirs)

% Calculates the differential pathlength factor (DPF) based on participant 
% age and wavelength. The DPF accounts for the increased photon path 
% length in tissue compared to the source-detector separation.

    % Extract age from timepoint
    ageYears = str2num(timepoint(1)) / 12;
    
    % Calculate differential pathlength factor (age and gamma dependent)
    [DPF1, DPF2] = sbPrePCalcDPF(ageYears, nirs.SD.Lambda(1), nirs.SD.Lambda(2));
    dpf = [DPF1 DPF2];
end
