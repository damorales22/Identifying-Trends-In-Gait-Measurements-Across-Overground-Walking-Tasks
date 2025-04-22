%% RunTrialData
% Loads the trial data and runs all relevant functions for trial data
clc
clear
close all

%% Variables
filepath = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\MATLAB scripts\Data\Trial1\";
numtrials = 7;
trialtypes = [0 0.15 0.20 0.25 0.35 0 0.05];

%% Create Data Struct
data = struct();
for i = 1:numtrials
    trialname = sprintf('Trial%04d', i);
    trial = load([num2str(filepath) num2str(trialname) '.mat']); % load trial
    
    try
        markers = trial.(trialname).Trajectories.Labeled.Labels;
    catch
        alternative_trialname = fieldnames(trial);
        trialname = alternative_trialname{1};
        markers = trial.(trialname).Trajectories.Labeled.Labels;
    end


    %% Separate markers
    markerdata = struct();
    for j = 1:length(markers)
       markerdata.(char(markers(j))) = trial.(num2str(trialname)).Trajectories.Labeled.Data(j, :, :);
    end
    markerdata.FrameRate = trial.(trialname).FrameRate;
    markerdata.Frames = trial.(trialname).Frames;

    data.markerdata.(trialname) = markerdata;
    
    % %% Locate Floor Markers
    % [points] = LocateMarkers(markerdata);
    % data.points.(trialname) = points;
    
    %% Get Steps
    delaytime = 1.5;
    %threshold = 3;
    try
        threshold = CalculateThreshold(trial, trialname);
    catch
        threshold = CalculateThreshold(trial, alternative_trialname{1}); %if there is a problem with the trial name
    end
    checktime = 0.2;
    steps = LocateSteps(markerdata, delaytime, threshold, checktime);
    data.steps.(trialname) = steps;
    
    % %% Run StepAccuracy
    % matches = StepError(points, steps, markerdata);
    % data.matches.(trialname) = matches;
    % data.meanerror_right(i) = matches.right.overall_mean;
    % data.meanerror_left(i) = matches.left.overall_mean;
    % data.meanerror(i) = matches.overall.mean;
    
    %% Run StepWidthStraights
    pos_limit = 3000;
    %neg_limit = 0;
    neg_limit = -1850;
    stepparameters = StepParametersHuxham(markerdata, steps, pos_limit, neg_limit);
    data.stepparameters.(trialname) = stepparameters;
    data.meanwidth(i) = stepparameters.stepwidth.straights_mean;
    data.variability(i) = stepparameters.stepwidth.straights_variability;
    data.RMS(i) = sqrt(sum(stepparameters.stepwidth.step_width_straights.^2)/length(stepparameters.stepwidth.step_width_straights));
    
    %% Run WalkingSpeed
    try
        walkingspeed = WalkingSpeed(trial, markerdata, trialname);
    catch
        walkingspeed = WalkingSpeed(trial, markerdata, alternative_trialname{1});
    end
    
    data.walkingspeed(i) = walkingspeed.avg;

end

%% Create CSV File
% Trial Numbers 
Trial = (1:numtrials)';

% step width
MeanWidth = data.meanwidth(1, :)';
WidthVariability = data.variability(1, :)';

% speed
AvgSpeed = data.walkingspeed';

% Create a table
T = table(Trial, MeanWidth, WidthVariability, AvgSpeed);

% Specify the file path
filename = filepath + 'data.csv';

% Write table T to CSV file
% writetable(T, filename);

%% Plot Width 

% Sort the data by x values
[trialtypes, sortOrder] = sort(trialtypes);
WidthVariability = WidthVariability(sortOrder);
RMS = data.RMS(sortOrder);

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
