%% RunTrialData
% Loads the trial data and runs all relevant functions for trial data
clc
clear
close all

%% Variables

filepath = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\MATLAB scripts\Data\Trial1\";
imported_data = load(filepath + "data.mat");
dims = size(imported_data.data.meanerror);
numtrials = dims(1, 2);
numtrials = 7;
trialtypes = [0 0.15 0.20 0.25 0.35 0 0.05];

for i=1:numtrials
    try
        last = trialtypes(i);
    catch
        trialtypes(end+1) = -1; %it means that the perturbation was not defined
    end
end

%% Plot Width 

WidthVariability = imported_data.data.straightsvariability(1, :)';
AvgSpeed = imported_data.data.walkingspeed';

% Sort the data by x values
[trialtypes, sortOrder] = sort(trialtypes);
WidthVariability = WidthVariability(sortOrder);

for i=1:numtrials
    trialname = sprintf('Trial%04d', i);
    RMS(i) = sqrt(sum(imported_data.data.step_parameters.(trialname).stepwidth.step_width_straights.^2)/length(imported_data.data.step_parameters.(trialname).stepwidth.step_width_straights));
end
RMS = RMS(sortOrder);

trialtypes = trialtypes(1, 2:end);
WidthVariability = WidthVariability(2:end, 1);
RMS = RMS(2:end);

% RMS Step Width
f = figure();
hold on
plot(trialtypes, RMS, 'Marker', '.','MarkerSize', 40, 'Color', "#0072BD", 'LineWidth', 1)

xlabel('Amplitude (m)', 'FontSize', 15)
ylabel('Mean Step Width (mm)', 'FontSize', 15)
title("RMS Step Width vs Perturbation Amplitude", 'FontSize', 20)
f.Position = [50 50 850 700];
axis("padded")
hold off


% Step Width Variability
f = figure();
hold on
plot(trialtypes, WidthVariability, 'Marker', '.','MarkerSize', 40, 'LineWidth', 1, 'Color', "#0072BD")

xlabel('Amplitude (m)', 'FontSize', 15)
ylabel('Mean Step Variability (mm)', 'FontSize', 15)
title("Step Width Variability vs Perturbation Amplitude", 'FontSize', 20)
f.Position = [50 50 850 700];
axis("padded")
hold off
