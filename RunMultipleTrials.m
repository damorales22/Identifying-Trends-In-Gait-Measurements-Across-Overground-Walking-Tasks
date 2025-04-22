% RunTrialData
% Loads the trial data and runs all relevant functions for trial data
clc
clear
close all
warning('off','all')
warning

%% Variables - BMH10
% filepath = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\BMH10\MC\";
% pos_limit = 3700;  %Positive and negative y axis limits only for the straight sections. Check it manually for each experiment
% neg_limit = -1970;
% height = 178; % In cm
% weight = 78; % kg
% age = 34;
% sex = "M";
% 
% trialspeeds =     [1 1 2 1 0 0 0 2 0 0 1 2 1 1 1 2 2 2 0 1 0 0 1 1 2 2 2 1 0 3 3 3 3 1 1]; % input trial speedprompts (0=slow, 1=medium, 2=fast)     3: idk
% trialaccuracies = [0 0 1 0 1 0 1 2 2 2 1 0 2 0 1 0 0 2 0 0 2 1 2 2 2 1 1 1 0 3 3 3 3 0 0]; % input trial accuracy prompts (0=none, 1=medium, 2=high)
% trialbalances =   [0 0 1 1 0 2 1 2 0 2 0 2 2 2 2 1 0 0 0 0 1 2 1 0 1 2 0 1 1 3 3 3 3 0 0]; % input trial balance conditions (0=no perturbation, 1 = medium, 2=high)
% 
% automatic_detection = 0; % 0 to always ask the user or use its preloaded prompts, 1 to do it automatically
% usedecisions = 1; % 1 to load the stored decisions to clasify steps, 0 to not do it (default)

%% Variables - BMH02
% filepath = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\BMH02\MC\";
% pos_limit = 4200;  %Positive and negative y axis limits only for the straight sections. Check it manually for each experiment
% neg_limit = -2800;
% height = 183; % In cm
% weight = 67; % kg
% age = 25;
% sex = "M";
% 
% trialspeeds =     [1 1 1 1 0 0 0 1 1 2 2 0 1 2 2 0 0 2 1 2 0 2 0 1 2 1 1 2 0 3 3 3 3 1 1]; % input trial speedprompts (0=slow, 1=medium, 2=fast)     3: idk
% trialaccuracies = [0 0 2 0 2 1 0 1 1 0 1 1 2 1 2 0 2 0 0 2 1 0 2 2 1 0 1 2 0 3 3 3 3 0 0]; % input trial accuracy prompts (0=none, 1=medium, 2=high)
% trialbalances =   [0 0 1 1 2 1 0 0 1 1 1 2 0 2 2 1 1 2 2 0 0 0 0 2 0 0 2 1 2 3 3 3 3 0 0]; % input trial balance conditions (0=no perturbation, 1 = medium, 2=high)
% 
% automatic_detection = 0; % 0 to always ask the user or use its preloaded prompts, 1 to do it automatically
% usedecisions = 1; % 1 to load the stored decisions to clasify steps, 0 to not do it (default)

%% Variables - BMH08
% filepath = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\BMH08\MC\";
% pos_limit = 4400;  %Positive and negative y axis limits only for the straight sections. Check it manually for each experiment
% neg_limit = -3050;
% height = 181; % In cm
% weight = 72; % kg
% age = 31;
% sex = "M";
% 
% trialspeeds =     [1 1 2 0 1 2 0 1 2 0 2 0 1 0 0 1 2 0 0 1 1 2 1 0 2 2 1 1 2 3 3 3 3 1 1]; % input trial speedprompts (0=slow, 1=medium, 2=fast) 3: idk
% trialaccuracies = [0 0 0 1 2 0 1 1 1 1 2 0 0 0 2 0 2 0 2 2 2 1 1 2 0 2 0 1 1 3 3 3 3 0 0]; % input trial accuracy prompts (0=none, 1=medium, 2=high)
% trialbalances =   [0 0 1 1 0 0 0 1 2 2 0 0 2 1 2 0 1 2 0 2 1 1 2 1 2 2 1 0 0 3 3 3 3 0 0]; % input trial balance conditions (0=no perturbation, 1 = medium, 2=high)
% 
% automatic_detection = 0; % 0 to always ask the user or use its preloaded prompts, 1 to do it automatically
% usedecisions = 1; % 1 to load the stored decisions to clasify steps, 0 to not do it (default)

%% Variables - BMH20
% filepath = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\BMH20\MC\";
% pos_limit = 4200;  %Positive and negative y axis limits only for the straight sections. Check it manually for each experiment
% neg_limit = -2800;
% height = 180; % In cm
% weight = 88.8; % kg
% age = 25;
% sex = "F";
% 
% trialspeeds =     [1 1 2 0 0 2 1 1 2 0 2 0 0 1 0 0 0 1 2 2 1 2 1 1 2 2 1 1 0 3 3 3 3 1 1]; % input trial speedprompts (0=slow, 1=medium, 2=fast)     3: idk
% trialaccuracies = [0 0 2 0 0 0 0 1 0 2 2 2 1 0 1 2 1 2 1 0 1 1 1 0 2 1 2 2 0 3 3 3 3 0 0]; % input trial accuracy prompts (0=none, 1=medium, 2=high)
% trialbalances =   [0 0 0 2 0 1 1 1 2 0 2 1 1 0 0 2 2 0 1 0 2 0 0 2 1 2 1 2 1 3 3 3 3 0 0]; % input trial balance conditions (0=no perturbation, 1 = medium, 2=high)
% 
% automatic_detection = 0; % 0 to always ask the user or use its preloaded prompts, 1 to do it automatically
% usedecisions = 1; % 1 to load the stored decisions to clasify steps, 0 to not do it (default)

%% Variables - BMH06
% filepath = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\BMH06\MC\";
% pos_limit = 4545;  %Positive and negative y axis limits only for the straight sections. Check it manually for each experiment
% neg_limit = -3125;
% height = 178; % In cm
% weight = 61.2; % kg
% age = 22;
% sex = "M";
% 
% trialspeeds =     [1 1 0 2 0 0 1 0 0 0 1 2 2 0 1 2 2 0 2 1 2 2 1 1 0 1 2 1 1 3 3 3 3 1 1]; % input trial speedprompts (0=slow, 1=medium, 2=fast)     3: idk
% trialaccuracies = [0 0 0 2 1 2 1 2 2 1 1 1 1 1 2 2 1 0 0 0 0 0 0 1 0 2 2 0 2 3 3 3 3 0 0]; % input trial accuracy prompts (0=none, 1=medium, 2=high)
% trialbalances =   [0 0 0 2 2 2 2 0 1 1 0 0 2 0 0 0 1 2 2 2 1 0 1 1 1 2 1 0 1 3 3 3 3 0 0]; % input trial balance conditions (0=no perturbation, 1 = medium, 2=high)
% 
% automatic_detection = 0; % 0 to always ask the user or use its preloaded prompts, 1 to do it automatically
% usedecisions = 1; % 1 to load the stored decisions to clasify steps, 0 to not do it (default)

%% Variables - BMH01
% filepath = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\BMH01\MC\";
% pos_limit = 3785;  %Positive and negative y axis limits only for the straight sections. Check it manually for each experiment
% neg_limit = -2930;
% height = 165; % In cm
% weight = 60; % kg
% age = 28;
% sex = "F";
% 
% trialspeeds =     [1 1 0 0 0 2 2 2 1 0 2 2 2 1 0 2 2 1 1 1 1 1 0 0 0 0 2 1 1 3 3 3 3 1 1]; % input trial speedprompts (0=slow, 1=medium, 2=fast)     3: idk
% trialaccuracies = [0 0 0 2 2 2 2 1 1 1 2 1 0 1 2 1 0 0 2 0 0 1 1 1 0 0 0 2 2 3 3 3 3 0 0]; % input trial accuracy prompts (0=none, 1=medium, 2=high)
% trialbalances =   [0 0 1 2 0 1 0 0 1 2 2 1 2 2 1 2 0 1 2 2 0 0 1 0 0 2 1 1 0 3 3 3 3 0 0]; % input trial balance conditions (0=no perturbation, 1 = medium, 2=high)
% 
% automatic_detection = 0; % 0 to always ask the user or use its preloaded prompts, 1 to do it automatically
% usedecisions = 1; % 1 to load the stored decisions to clasify steps, 0 to not do it (default)

%% Variables - BMH19
% filepath = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\BMH19\MC\";
% pos_limit = 4500;  %Positive and negative y axis limits only for the straight sections. Check it manually for each experiment
% neg_limit = -2850;
% height = 173; % In cm
% weight = 85; % kg
% age = 22;
% sex = "F";
% 
% trialspeeds =     [1 1 2 1 0 1 2 0 0 1 1 0 2 0 1 0 1 1 2 1 2 0 0 2 1 0 2 2 2 3 3 3 3 1 1]; % input trial speedprompts (0=slow, 1=medium, 2=fast)     3: idk
% trialaccuracies = [0 0 1 0 0 2 0 0 1 2 1 2 0 1 0 2 2 0 2 1 0 1 0 2 1 2 1 1 2 3 3 3 3 0 0]; % input trial accuracy prompts (0=none, 1=medium, 2=high)
% trialbalances =   [0 0 2 1 0 2 2 1 2 0 1 2 0 1 2 1 1 0 0 0 1 0 2 1 2 0 0 1 2 3 3 3 3 0 0]; % input trial balance conditions (0=no perturbation, 1 = medium, 2=high)
% 
% automatic_detection = 0; % 0 to always ask the user or use its preloaded prompts, 1 to do it automatically
% usedecisions = 1; % 1 to load the stored decisions to clasify steps, 0 to not do it (default)

%% Variables - BMH17
% filepath = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\BMH17(yo)\MC\";
% pos_limit = 4300;  %Positive and negative y axis limits only for the straight sections. Check it manually for each experiment
% neg_limit = -2400;
% height = 172; % In cm
% weight = 52; % kg
% age = 24;
% sex = "M";
% 
% trialspeeds =     [1 1 2 0 0 0 2 1 0 2 0 1 0 2 1 2 1 1 1 1 0 2 2 2 1 0 2 1 0 3 3 3 3 1 1]; % input trial speedprompts (0=slow, 1=medium, 2=fast)     3: idk
% trialaccuracies = [0 0 1 1 2 1 2 1 1 1 2 0 0 0 0 1 2 2 1 1 0 0 2 2 0 0 0 2 2 3 3 3 3 0 0]; % input trial accuracy prompts (0=none, 1=medium, 2=high)
% trialbalances =   [0 0 0 0 0 2 0 2 1 2 2 2 2 0 0 1 1 0 1 0 0 1 1 2 1 1 2 2 1 3 3 3 3 0 0]; % input trial balance conditions (0=no perturbation, 1 = medium, 2=high)
% 
% automatic_detection = 0; % 0 to always ask the user or use its preloaded prompts, 1 to do it automatically
% usedecisions = 1; % 1 to load the stored decisions to clasify steps, 0 to not do it (default)

%% Variables - BMH09
% filepath = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\BMH09\MC\";
% pos_limit = 4300;  %Positive and negative y axis limits only for the straight sections. Check it manually for each experiment
% neg_limit = -2400;
% height = 159; % In cm
% weight = 80; % kg
% age = 23;
% sex = "F";
% 
% trialspeeds =     [1 1 0 2 0 1 1 2 1 2 0 2 0 0 0 2 2 1 0 1 0 2 2 2 0 1 1 1 1 3 3 3 3 1 1]; % input trial speedprompts (0=slow, 1=medium, 2=fast)     3: idk
% trialaccuracies = [0 0 0 2 1 2 0 1 1 1 2 1 2 1 0 0 2 1 1 1 2 2 0 0 0 0 0 2 2 3 3 3 3 0 0]; % input trial accuracy prompts (0=none, 1=medium, 2=high)
% trialbalances =   [0 0 0 1 2 2 1 0 2 2 1 1 0 0 1 1 0 0 1 1 2 2 0 2 2 2 0 0 1 3 3 3 3 0 0]; % input trial balance conditions (0=no perturbation, 1 = medium, 2=high)
% 
% automatic_detection = 0; % 0 to always ask the user or use its preloaded prompts, 1 to do it automatically
% usedecisions = 1; % 1 to load the stored decisions to clasify steps, 0 to not do it (default)

%% Variables - BMH13
% filepath = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\BMH13\MC\";
% pos_limit = 4500;  %Positive and negative y axis limits only for the straight sections. Check it manually for each experiment
% neg_limit = -2560;
% height = 183; % In cm
% weight = 65; % kg
% age = 26;
% sex = "M";
% 
% trialspeeds =     [1 1 0 2 1 1 1 0 0 2 0 0 1 2 1 0 1 1 0 0 1 2 2 1 2 2 2 2 0 3 3 3 3 1 1]; % input trial speedprompts (0=slow, 1=medium, 2=fast)     3: idk
% trialaccuracies = [0 0 2 0 0 0 2 0 0 0 1 2 2 0 1 1 1 1 2 1 0 1 1 2 2 2 2 1 0 3 3 3 3 0 0]; % input trial accuracy prompts (0=none, 1=medium, 2=high)
% trialbalances =   [0 0 2 1 2 0 0 2 0 0 1 1 2 2 0 2 1 2 0 0 1 2 1 1 1 2 0 0 1 3 3 3 3 0 0]; % input trial balance conditions (0=no perturbation, 1 = medium, 2=high)
% 
% automatic_detection = 0; % 0 to always ask the user or use its preloaded prompts, 1 to do it automatically
% usedecisions = 1; % 1 to load the stored decisions to clasify steps, 0 to not do it (default)

%% Variables - BMH07
filepath = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\BMH07\MC\";
pos_limit = 4170;  %Positive and negative y axis limits only for the straight sections. Check it manually for each experiment
neg_limit = -2775;
height = 175; % In cm
weight = 64; % kg
age = 25;
sex = "M";

trialspeeds =     [1 1 2 0 2 1 1 0 2 1 0 1 1 0 2 1 1 2 2 2 1 0 0 0 0 1 0 2 2 3 3 3 3 1 1]; % input trial speedprompts (0=slow, 1=medium, 2=fast)     3: idk
trialaccuracies = [0 0 0 1 2 1 0 0 0 1 1 2 2 1 1 2 1 1 0 1 0 2 0 0 2 0 2 2 2 3 3 3 3 0 0]; % input trial accuracy prompts (0=none, 1=medium, 2=high)
trialbalances =   [0 0 0 0 1 1 2 2 2 0 1 0 1 2 2 2 2 0 1 1 0 0 0 1 1 1 2 2 0 3 3 3 3 0 0]; % input trial balance conditions (0=no perturbation, 1 = medium, 2=high)

automatic_detection = 0; % 0 to always ask the user or use its preloaded prompts, 1 to do it automatically
usedecisions = 1; % 1 to load the stored decisions to clasify steps, 0 to not do it (default)

%% Load data

% Look for files ending with '.mat'
matFiles = dir(fullfile(filepath, '*Trial*.mat'));

% Count the number of .mat files
numtrials = length(matFiles); 

try
    static = load(filepath + 'static.mat');
catch
    static = load(filepath + 'static0001.mat');
end

counter = 1;
bmh = regexp(filepath, 'BMH\d{2}', 'match');

%% Create Data Struct
data = struct(); %where the relevant data of all trials are stored
missing = [];
extra = [];
for i=1:numtrials
    trialname = sprintf('Trial%04d', i);
    trial = load([num2str(filepath) char(trialname) '.mat']); % load trial
    try
        markers = trial.(trialname).Trajectories.Labeled.Labels;
        floor_markers = static.static0001.Trajectories.Labeled.Labels; %floor and non floor markers; separated on LocateMarkers
    catch
        alternative_trialname = fieldnames(trial);
        trialname = alternative_trialname{1};

        markers = trial.(trialname).Trajectories.Labeled.Labels;
        try
            floor_markers = static.static.Trajectories.Labeled.Labels; %floor and non floor markers; separated on LocateMarkers
        catch
            floor_markers = static.static0001.Trajectories.Labeled.Labels;
        end
    end

    %% Correct mistakes in the data structure of each participant
    if bmh == "BMH20" && floor_markers(32) == "New 0031"
        floor_markers{32} = 'RTH4';
    end
    
    %% Separate markers
    markerdata = struct();

    for j = 1:length(markers)
       markerdata.(char(markers(j))) = trial.(num2str(trialname)).Trajectories.Labeled.Data(j, :, :); %stores in the struct only the data from each label (1x4xFrames)
    end
    
    try % If it is named static0001
        for j = 1:length(floor_markers) 
           floordata.(char(floor_markers(j))) = static.static0001.Trajectories.Labeled.Data(j, :, :); %floor and non floor markers; separated on LocateMarkers
        end
    catch % If it is named static
        for j = 1:length(floor_markers) 
           floordata.(char(floor_markers(j))) = static.static.Trajectories.Labeled.Data(j, :, :);
        end       
    end
    markerdata.FrameRate = trial.(trialname).FrameRate;
    markerdata.Frames = trial.(trialname).Frames;
    floordata.FrameRate = trial.(trialname).FrameRate;
    floordata.Frames = trial.(trialname).Frames;

    trialname = sprintf('Trial%04d', i); %To fix the alternative_trialname, if it was the case
    %data.markerdata.(trialname) = markerdata; % Don't store this in the output file (huge amount of useless data once it's processed)
    
    %% Calculate Threshold
    %try
        %threshold = CalculateThreshold(trial, trialname);
    %catch
        %threshold = CalculateThreshold(trial, alternative_trialname{1}); %if there is a problem with the trial name
    %end
    %all_thesholds(i,1) = threshold;

    %% Butterworth Filter for BMH02, 08, 10
    % For the rest it was applied on Qualysis, but for these do manually (with the same filter parameters)
    if bmh == "BMH02" || bmh == "BMH08" || bmh == "BMH10"
        fs = markerdata.FrameRate;  % Sampling frequency (Hz)
        fc = 6;    % Cutoff frequency (Hz)
        order = 2; % Filter order
        
        [b, a] = butter(order, fc/(fs/2), 'low'); % Butterworth low-pass filter
        
        fieldNames = fieldnames(markerdata);
        
        for m = 1:length(fieldNames)
            fieldName = fieldNames{m};
            
            % Check if the field is a 1x4xn structure
            if isfield(markerdata, fieldName) && isnumeric(markerdata.(fieldName)) ...
                    && ndims(markerdata.(fieldName)) == 3 && size(markerdata.(fieldName), 2) == 4
                
                data_filtered = markerdata.(fieldName);
                
                % Apply filter only to x, y, z (first three rows)
                for n = 1:3
                    raw_signal = squeeze(data_filtered(1, n, :)); 
                    nan_idx = isnan(raw_signal); % Identify NaN positions
                    
                    % Filter only non-NaN values
                    valid_signal = raw_signal(~nan_idx);
                    if ~isempty(valid_signal)
                        filtered_signal = filtfilt(b, a, valid_signal);
                        raw_signal(~nan_idx) = filtered_signal; % Restore filtered data
                    end
                    data_filtered(1, n, :) = raw_signal; % Save back to structure
                end
                markerdata.(fieldName) = data_filtered;
            end
        end
    end

    %% Locate Floor Markers
    [points] = LocateMarkers(markerdata,floordata);
    data.points.(trialname) = points;
    
    %% Get Steps
    delaytime = 1.5;
    %steps = LocateSteps(markerdata, delaytime,threshold, 0.2, trialname); %the last term is checktime. It detects when the ball of the foot is firmly on the ground and not moving (more or less when the foot is fully horizontal)
    steps = LocateSteps(trial, trialname, counter, usedecisions, automatic_detection, filepath);
    counter = steps.counter;
    data.steps.(trialname) = steps;
    close all
    
    %% Run StepError
    matches = StepError(points, steps, markerdata, trialname,filepath);
    data.matches.(trialname) = matches;
    data.meanerror_right(i) = data.matches.(trialname).right.overall_mean; %Longitudinal absolute error
    data.meanerror_left(i) = data.matches.(trialname).left.overall_mean;  
    data.meanerror(i) = data.matches.(trialname).overall.overall_mean_absolute;
    data.meanerror_median(i) = data.matches.(trialname).overall.straight_median;

    %Total, mediolateral and longitudinal errors; Total error = data.matches.right.match(:,4); Longitudinal error = data.matches.right.match(:,3). It is a triangle
    data.matches.(trialname).right.match(:,6) = real(sqrt(data.matches.(trialname).right.match(:,5).^2 - data.matches.(trialname).right.match(:,3).^2)); %Mediolateral error
    data.matches.(trialname).left.match(:,6) = real(sqrt(data.matches.(trialname).left.match(:,5).^2 - data.matches.(trialname).left.match(:,3).^2));

    %For straight and curved sections
    data.matches.(trialname).right.match_straight(:,6) = real(sqrt(data.matches.(trialname).right.match_straight(:,5).^2 - data.matches.(trialname).right.match_straight(:,3).^2)); %Mediolateral error
    data.matches.(trialname).left.match_straight(:,6) = real(sqrt(data.matches.(trialname).left.match_straight(:,5).^2 - data.matches.(trialname).left.match_straight(:,3).^2));

    data.matches.(trialname).right.match_curve(:,6) = real(sqrt(data.matches.(trialname).right.match_curve(:,5).^2 - data.matches.(trialname).right.match_curve(:,3).^2)); %Mediolateral error
    data.matches.(trialname).left.match_curve(:,6) = real(sqrt(data.matches.(trialname).left.match_curve(:,5).^2 - data.matches.(trialname).left.match_curve(:,3).^2));

    %% Run StepParameters
    step_parameters = StepParametersHuxham(markerdata, steps, pos_limit, neg_limit); %Check pos_limit and neg_limit
    data.step_parameters.(trialname) = step_parameters;

    data.meanwidthstraights(i) = step_parameters.stepwidth.straights_mean;
    data.meanwidthstraights_median(i) = step_parameters.stepwidth.straights_median;
    data.meanwidthcurves(i) = step_parameters.stepwidth.curves_mean;
    data.straightsvariability(i) = step_parameters.stepwidth.straights_variability;
    data.curvesvariability(i) = step_parameters.stepwidth.curves_variability;
    
    data.meanlength_straights(i) = step_parameters.steplength.straights_mean;
    data.meanlength_straights_median(i) = step_parameters.steplength.straights_median;
    data.meanlength_curves(i) = step_parameters.steplength.curves_mean;
    data.meanlength_straights_variability(i) = step_parameters.steplength.straights_variability;
    data.meanlength_curves_variability(i) = step_parameters.steplength.curves_variability;
    
    %% Run WalkingSpeed
    if trialname == "Trial0015"
        a=0;
    end
    try
        [walkingspeed, stride_time] = WalkingSpeed(trial, markerdata, trialname,steps);
    catch
        [walkingspeed, stride_time] = WalkingSpeed(trial, markerdata, alternative_trialname{1},steps); %if there is a problem with the trial name
    end
    data.walkingspeed(i) = walkingspeed.avg;
    data.walkingspeed_std(i) = walkingspeed.std;
    data.stride_time(i) = stride_time.avg;
    data.stride_std(i) = stride_time.std;
    data.stride_time_median(i) = stride_time.median;

    %% Run Head_back_angles
    if i>=3 && i<=29
        try
            [head_angles, stj_headset_angles, trunk_angles] = Head_back_angles(trial, markerdata, trialname, static.static0001.Trajectories.Labeled, bmh);
        catch
            [head_angles, stj_headset_angles, trunk_angles] = Head_back_angles(trial, markerdata, trialname, static.static.Trajectories.Labeled, bmh);
        end
        data.head_angles.(trialname) = head_angles;
        data.stj_headset_angles.(trialname) = stj_headset_angles;
        data.trunk_angles.(trialname) = trunk_angles;
    end

    %% Where are they looking at?
    if i>=3 && i<=29
        try
            [dist_gaze] = Dist_gaze(trial, markerdata, trialname, static.static0001.Trajectories.Labeled, data.head_angles.(trialname), data.trunk_angles.(trialname), height, bmh);
        catch
            [dist_gaze] = Dist_gaze(trial, markerdata, trialname, static.static.Trajectories.Labeled, data.head_angles.(trialname), data.trunk_angles.(trialname), height, bmh);
        end
        data.dist_gaze.(trialname) = dist_gaze;
    end

    %% Locate Heel Strikes
    %try
        %strikes = LocateHeelStrikes(trial, trialname, steps);
    %catch
        %strikes = LocateHeelStrikes(trial, alternative_trialname{1}, threshold);
        %fprint("Other file name used at iteration " + i)
    %end
    %data.strikes.(trialname) = strikes;

    %% Create Speed, Accuracy, and Balance Columns. Create labels for each trial
    % Speeds
    try
        if trialspeeds(i) == 0
            Speed(i, 1) = "Slow";
        elseif trialspeeds(i) == 1
            Speed(i, 1) = "Medium";
        elseif trialspeeds(i) == 2
            Speed(i, 1) = "Fast";
        elseif trialspeeds(i) == 3
            Speed(i, 1) = "Unknown";
        end
    catch
        trialspeeds(end+1) = 3;
        Speed(i, 1) = "Unknown";
    end
    % Accuracies
    try
        if trialaccuracies(i) == 0
            Accuracy(i, 1) = "Low";
        elseif trialaccuracies(i) == 1
            Accuracy(i, 1) = "Medium";
        elseif trialaccuracies(i) == 2
            Accuracy(i, 1) = "High";
        elseif trialaccuracies(i) == 3
            Accuracy(i, 1) = "Unknown";
        end
     catch
        trialaccuracies(end+1) = 3;
        Accuracy(i, 1) = "Unknown";
    end
     % Balances
     try
        if trialbalances(i) == 0
            Balance(i, 1) = "None";
        elseif trialbalances(i) == 1
            Balance(i, 1) = "Medium";
        elseif trialbalances(i) == 2
            Balance(i, 1) = "High";
        elseif trialbalances(i) == 3
            Balance(i, 1) = "Unknown";
        end
     catch
        trialbalances(end+1) = 3;
        Balance(i, 1) = "Unknown";
     end

    % Create labels
    data.labels(i) = "s" + string(trialspeeds(i)) + "a" + string(trialaccuracies(i)) + "b" + string(trialbalances(i));

    fprintf('Trial %d / %d analyzed\n', i, numtrials);
end %end of data gathering of each trial

%% Store other parameters
Trial = (1:numtrials)'; % Trial Numbers 

%Step error 
for i=1:numtrials
    trialname = sprintf('Trial%04d', i);
    MeanError_straight(i,1) = data.matches.(trialname).overall.straight_mean_absolute; %Absolute longitudinal error in the straight sections
end
MeanError_r_mm = data.meanerror_right(1, :)'; %Absolute longitudinal error
MeanError_l_mm = data.meanerror_left(1, :)';  %Absolute longitudinal error
MeanError_mm = (MeanError_r_mm + MeanError_l_mm) / 2; %Average of right and left
MeanError_mm_median = data.meanerror_median(1, :)';

%Step width
MeanWidthStraights_mm = data.meanwidthstraights(1, :)';
StraightsWidthVariability_mm = data.straightsvariability(1, :)';
MeanWidthStraights_mm_median = data.meanwidthstraights_median(1, :)';
MeanWidthCurves_mm = data.meanwidthcurves(1, :)';
CurvesWidthVariability_mm = data.curvesvariability(1, :)';

% Step length
MeanLength_straights = data.meanlength_straights(1, :)';
MeanLength_straights_variability = data.meanlength_straights_variability(1, :)';
MeanLength_straights_median = data.meanlength_straights_median(1, :)';
MeanLength_curves = data.meanlength_curves(1, :)';
MeanLength_curves_variability = data.meanlength_curves_variability(1, :)';

%Speed
AvgSpeed_m_per_s = data.walkingspeed(1,:)';

%Prompts
data.prompt.speed = Speed;
data.prompt.speed(:,2) = trialspeeds';
data.prompt.accuracy = Accuracy;
data.prompt.accuracy(:,2) = trialaccuracies';
data.prompt.balance = Balance;
data.prompt.balance(:,2) = trialbalances';

% Participant data
data.participant_data.age = age;
data.participant_data.height = height;
data.participant_data.weight = weight;
data.participant_data.sex = sex;

%% Overall errors across all the trials

% Store all errors in array format so that it is easier to compare them
for i=1:numtrials
    trial_name = sprintf('Trial%04d', i);
    % if trial_name == "Trial0004"
    %     i = i+1;
    %     continue
    % end
    overall_mean_absolute_array (i) = data.matches.(trial_name).overall.overall_mean_absolute;
    overall_std_absolute_array (i) = data.matches.(trial_name).overall.overall_std_absolute;
    overall_mean_signed_array (i) = data.matches.(trial_name).overall.overall_mean_signed;
    overall_std_signed_array (i) = data.matches.(trial_name).overall.overall_std_signed;
    straight_mean_absolute_array (i) = data.matches.(trial_name).overall.straight_mean_absolute;
    straight_std_absolute_array (i) = data.matches.(trial_name).overall.straight_std_absolute;
    curve_mean_absolute_array (i) = data.matches.(trial_name).overall.curve_mean_absolute;
    curve_std_absolute_array (i) = data.matches.(trial_name).overall.curve_std_absolute;
    straight_mean_signed_array (i) = data.matches.(trial_name).overall.straight_mean_signed;
    straight_std_signed_array (i) = data.matches.(trial_name).overall.straight_std_signed;
    curve_mean_signed_array (i) = data.matches.(trial_name).overall.curve_mean_signed;
    curve_std_signed_array (i) = data.matches.(trial_name).overall.curve_std_signed;
end

data.overall_errors.overall_mean_absolute = overall_mean_absolute_array;
data.overall_errors.overall_std_absolute = overall_std_absolute_array;
data.overall_errors.overall_mean_signed = overall_mean_signed_array;
data.overall_errors.overall_std_signed = overall_std_signed_array;
data.overall_errors.straight_mean_absolute = straight_mean_absolute_array;
data.overall_errors.straight_std_absolute = straight_std_absolute_array;
data.overall_errors.curve_mean_absolute = curve_mean_absolute_array;
data.overall_errors.curve_std_absolute = curve_std_absolute_array;
data.overall_errors.straight_mean_signed = straight_mean_signed_array;
data.overall_errors.straight_std_signed = straight_std_signed_array;
data.overall_errors.curve_mean_signed = curve_mean_signed_array;
data.overall_errors.curve_std_signed  = curve_std_signed_array;

% Take only the trials in which there are no balance perturbations and high accuracy
overall_mean_absolute_array = [];
overall_std_absolute_array = [];
overall_mean_signed_array = [];
overall_std_signed_array = [];
straight_mean_absolute_array = [];
straight_std_absolute_array = [];
curve_mean_absolute_array = [];
curve_std_absolute_array = [];
straight_mean_signed_array = [];
straight_std_signed_array = [];
curve_mean_signed_array = [];
curve_std_signed_array = [];
labels_balance0 = {}; %Variable for control

for i=1:numtrials
    trial_name = sprintf('Trial%04d', i);
    if str2double(data.prompt.balance(i,2)) == 0 && str2double(data.prompt.accuracy(i,2)) == 2 && (i>=3 && i<=29) %no perturbations and high accuracy
        overall_mean_absolute_array (end+1) = data.matches.(trial_name).overall.overall_mean_absolute;
        overall_std_absolute_array (end+1) = data.matches.(trial_name).overall.overall_std_absolute;
        overall_mean_signed_array (end+1) = data.matches.(trial_name).overall.overall_mean_signed;
        overall_std_signed_array (end+1) = data.matches.(trial_name).overall.overall_std_signed;
        straight_mean_absolute_array (end+1) = data.matches.(trial_name).overall.straight_mean_absolute;
        straight_std_absolute_array (end+1) = data.matches.(trial_name).overall.straight_std_absolute;
        curve_mean_absolute_array (end+1) = data.matches.(trial_name).overall.curve_mean_absolute;
        curve_std_absolute_array (end+1) = data.matches.(trial_name).overall.curve_std_absolute;
        straight_mean_signed_array (end+1) = data.matches.(trial_name).overall.straight_mean_signed;
        straight_std_signed_array (end+1) = data.matches.(trial_name).overall.straight_std_signed;
        curve_mean_signed_array (end+1) = data.matches.(trial_name).overall.curve_mean_signed;
        curve_std_signed_array (end+1) = data.matches.(trial_name).overall.curve_std_signed;
        labels_balance0{end+1} = trial_name;  
    end
end

data.overall_errors.balance0.overall_mean_absolute = mean(overall_mean_absolute_array);
data.overall_errors.balance0.overall_std_absolute = mean(overall_std_absolute_array);
data.overall_errors.balance0.overall_mean_signed = mean(overall_mean_signed_array);
data.overall_errors.balance0.overall_std_signed = mean(overall_std_signed_array);
data.overall_errors.balance0.straight_mean_absolute = mean(straight_mean_absolute_array);
data.overall_errors.balance0.straight_std_absolute = mean(straight_std_absolute_array);
data.overall_errors.balance0.curve_mean_absolute = mean(curve_mean_absolute_array);
data.overall_errors.balance0.curve_std_absolute = mean(curve_std_absolute_array);
data.overall_errors.balance0.straight_mean_signed = mean(straight_mean_signed_array);
data.overall_errors.balance0.straight_std_signed = mean(straight_std_signed_array);
data.overall_errors.balance0.curve_mean_signed = mean(curve_mean_signed_array);
data.overall_errors.balance0.curve_std_signed  = mean(curve_std_signed_array);

%% Store user decisions in a separate txt file
decisions = strings(1, 0);
if usedecisions == 0 && automatic_detection == 0
    for i=1:numtrials
        trialname = sprintf('Trial%04d', i);
        decisions = [decisions, data.steps.(trialname).decisions];
    end
    fileID = fopen(filepath + 'decisions.txt', 'a');
    for j = 1:length(decisions)
        fprintf(fileID, '%s', decisions(j));
    end
    fclose(fileID);
end

MeanLength_straights = data.meanlength_straights(1, :)';
MeanLength_straights_variability = data.meanlength_straights_variability(1, :)';
MeanLength_curves = data.meanlength_curves(1, :)';
MeanLength_curves_variability = data.meanlength_curves_variability(1, :)';

%% Write table T to CSV file
% Create table
T = table(Trial, Speed, Accuracy, Balance, MeanError_r_mm, MeanError_l_mm, MeanError_mm, MeanError_straight, MeanError_mm_median, MeanWidthStraights_mm, MeanWidthStraights_mm_median, ...
    StraightsWidthVariability_mm, MeanWidthCurves_mm, CurvesWidthVariability_mm, MeanLength_straights, MeanLength_straights_median, MeanLength_straights_variability, ...
    MeanLength_curves, MeanLength_curves_variability, AvgSpeed_m_per_s);

T.Properties.VariableNames = {'Trial Number', 'Walking Speed', 'Accuracy', 'Balance', 'Mean Error Right All Sections (mm)', 'Mean Error Left All Sections (mm)', ...
                              'Mean Error All Sections (mm)', 'Mean Error Straights', 'Median Error Straights', 'Mean Width Straights (mm)', 'Median Width Straights (mm)', 'Straights Width Variability (mm)', 'Mean Width Curves (mm)', ...
                              'Curves Width Variability (mm)', 'Mean Length Straights (mm)', 'Median Length Straights (mm)', 'Straights Length Variability (mm)', ...
                              'Mean Length Curves (mm)', 'Curves Length Variability (mm)', 'Average Speed (m/s)'};

% Specify the file path
filename = filepath + 'data_' + bmh + '.csv';
filename2 = filepath + 'data_' + bmh + '.xlsx';

if exist(filename, 'file') == 2 || exist(filename2, 'file') == 2
    % File exists
    choice = questdlg('The file already exists. Do you want to delete it and continue?', ...
                      'File Exists', ...
                      'Delete and Continue', 'Cancel', 'Cancel');
    if strcmp(choice, 'Delete and Continue')
        % Delete the existing file
        delete(filename);
        fprintf('Existing file deleted.\n');
    else
        fprintf('Operation cancelled.\n');
        return;  % Exit the script or function
    end
end

writetable(T, filename); %To import: name=readtable('C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\MATLAB scripts\Data\Trial1\data.csv'); To retrieve data: name{column, row}
writetable(T, filename2); 
disp(filename + " written")
save(filepath + 'data_' + bmh +'.mat', 'data');



