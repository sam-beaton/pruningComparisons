clear; close all

% Read the table from the CSV file
data = readtable('[path-to...]/pruningComparisonsMatlab/stats/bright/overall/pruneChannelInfoTableCV.csv');

%% Average Channels and ROI_SNR by Age
% Calculate the average 'Channels' and 'ROI_SNR' for each unique value of 'Age'
averageData = groupsummary(data, 'Age', 'mean', {'Channels', 'ROI_SNR'});

% Display the results
disp(averageData);