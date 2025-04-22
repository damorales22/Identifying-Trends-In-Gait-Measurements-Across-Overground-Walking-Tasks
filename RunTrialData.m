%% RunTrialData
% Loads the trial data and runs all relevant functions for trial data
clc
clear
close all
%% Variables
trialname = "Trial0012";

% for Locate Steps function
delaytime = 1.5;
checktime = 0.2;

% for Step Parameters Huxham function
%pos_limit = 3750;
%neg_limit = -1850;
pos_limit = 3700;
neg_limit = -1970;

%% Loads Trial Data
%trial = load(['/Users/sarahli/Documents/Ability Lab/7_24_Data/' num2str(trialname) '.mat']); % load trial
trial = load(['C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\Trial1\MC\'+(trialname)+'.mat']); % load trial
static = load(['C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\Trial1\MC\static0001.mat']);
markers = trial.(trialname).Trajectories.Labeled.Labels; %gets the labels used in that trial
floor_markers = static.static0001.Trajectories.Labeled.Labels; %floor and non floor markers; separated on LocateMarkers

%% Separates out markers
markerdata = struct();
floordata = struct();
for i = 1:length(markers)
   markerdata.(char(markers(i))) = trial.(num2str(trialname)).Trajectories.Labeled.Data(i, :, :); %stores in the struct only the data from each label (1x4xFrames)
end
for i = 1:length(floor_markers) 
   floordata.(char(floor_markers(i))) = static.static0001.Trajectories.Labeled.Data(i, :, :);
end

markerdata.FrameRate = trial.(trialname).FrameRate; %stores the frame rate of that experiment
markerdata.Frames = trial.(trialname).Frames; %stores the number of frames of that experiment
floordata.FrameRate = trial.(trialname).FrameRate;
floordata.Frames = trial.(trialname).Frames;

%% Calculate Threshold
try
    threshold = CalculateThreshold(trial, trialname);
catch
    threshold = CalculateThreshold(trial, alternative_trialname{1}); %if there is a problem with the trial name
end

%% Locate Floor Markers
[points] = LocateMarkers(markerdata, floordata); %returns the coordinates of the outer and iner markers (red tape) and their middle point

%% Get Steps
steps = LocateSteps(markerdata, delaytime, threshold, checktime); %returns the timestamps (and frames) of every step taken for both feet.

%% Runs StepAccuracy
[matches]= StepError(points, steps, markerdata);

%Plot the circuit and the steps of the right foot
hold on
for i=1:length(fieldnames(points))
    plot(points.("FL" + num2str(i))(1,1:2), points.("FL" + num2str(i))(2,1:2), 'o')
    %plot(points.("FL" + num2str(i))(1:2,1), points.("FL" + num2str(i))(1:2,2), 'o')
end
legend()
for i=1:length(matches.right.match)
    timestamp = matches.right.match(i,1);
    step = markerdata.RMTH5(:, 1:2, timestamp);
    plot(step(1,1), step(1,2), 'ro', 'MarkerFaceColor', 'r')
end 
axis equal
axis square;
title('Right foot and floor markers distribution')
hold off

%% Runs StepWidth
[step_parameters] = StepParametersHuxham(markerdata, steps, pos_limit, neg_limit);  

%% Runs WalkingSpeed
[speed_avg] = WalkingSpeed(trial,markerdata, trialname);
