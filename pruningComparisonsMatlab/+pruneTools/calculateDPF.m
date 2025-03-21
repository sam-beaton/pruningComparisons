function dpf = calculateDPF(age, nirs, convertMonthsToYears)
% Calculates the differential pathlength factor (DPF) based on participant 
% age and wavelength. The DPF accounts for the increased photon path 
% length in tissue compared to the source-detector separation.
%
% Based on Scholkmann and Wolf (2013)
%
% Inputs:
%   age:            in years (or months if convertToYears is true)
%   nirs:           NIRS data structure containing SD.Lambda
%   convertToYears: (optional) boolean, if true, converts age in months to years

    if nargin < 3  % Check if the third argument is provided
        convertMonthsToYears = false; % Default to false
    end

    if convertMonthsToYears
        age = age(1) / 12;
    end

    Alph = 223.3;
    Bet = 0.05624;
    Gam = 0.8493;
    Delt = -0.0000005723;
    Eps = 0.001245;
    Zet = -0.9025;

    wav1 = nirs.SD.Lambda(1);
    wav2 = nirs.SD.Lambda(2);
    
    DPF1 = Alph + Bet*(age^Gam) + Delt*(wav1^3) + Eps*(wav1^2) + Zet*wav1;
    DPF2 = Alph + Bet*(age^Gam) + Delt*(wav2^3) + Eps*(wav2^2) + Zet*wav2;
    
    dpf = [DPF1 DPF2];
end

