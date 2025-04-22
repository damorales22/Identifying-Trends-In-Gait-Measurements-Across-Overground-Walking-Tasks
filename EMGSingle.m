%% Segmenting EMG data by Gait Cycle
% This script loads a single trial, locates heel strikes, segments the EMG data
% based on these heel strike times, plots gait cycles on top of each other
% and then averages them
clc 
clear
close all

%% Variables
%filepath = "C:\Ability Lab\7_31_2024\";
%EMGpath = "C:\Ability Lab\7_31_2024\EMG\";

EMGpath =  "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\BMH10(September)\EMG\";
MCfile = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\BMH10(September)\MC\data_BMH10.mat";

trialnumber = 12;
label = "RSOL"; % choose the EMG channel to look at
delaytime = 1.5;
checktime = 0.2;
upperthreshold = 12;
maxsteptime = 3; % set maximum step length for segment

% bandpass filter inputs
f_low = 30;
f_high = 500;
n = 1; % second order (filtfilt doubles it to fourth order)

% lowpass filter inputs
f_lowpass = 6;
m = 2; % second order (filtfilt doubles it to fourth order)

%% Load Trial Data

EMGdata = struct();

% Look for files ending with '.mat'
matFiles = dir(fullfile(EMGpath, '*.mat'));

% Count the number of .mat files
numfiles = length(matFiles); 

if ismember('data.mat', {matFiles.name})
    numtrials = numfiles - 6;  % Subtract 1 extra for 'data.mat'
else
    numtrials = numfiles - 5;  % Subtract 5 for the MVCx files (max contraction)
end

for i=1:numtrials
    trialname = sprintf('Trial%04d', i);
    EMGdata.(trialname) = load(EMGpath + trialname + '.mat');
end

MVCfiles = dir(fullfile(EMGpath, 'MVC*.mat'));
for i = 1:length(MVCfiles)
    fileName = MVCfiles(i).name;
    fullFilePath = fullfile(EMGpath, fileName);
    EMGdata.(sprintf('MVC%d', i)) = load(fullFilePath);
end

MCdata = load(MCfile);
%aaaaa = load("C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\Trial1\MC\Trial0012")

trialname = sprintf('Trial%04d', trialnumber);
filename = EMGpath + trialname;
EMGname = EMGpath + trialname;
trial = load([num2str(filename) '.mat']);
data = load([num2str(EMGname) '.mat']);
fs = EMGdata.(trialname).Fs(1);
% framerate = EMGdata.(trialname).FrameRate; % get frame rate
% numframes = EMGdata.(trialname).Frames;
framerate = fs;
dim_numframes = size(EMGdata.(trialname).Data);
numframes = dim_numframes(1,2); %different for each trial
trialtime = numframes / framerate;

if startsWith(label, 'R')
    foot = 'RHEEL';
    side = "right";
else
    foot = 'LHEEL';
    side = "left";
end

%% Separates out markers - MoCap
% markers = trial.(trialname).Trajectories.Labeled.Labels;
% markerdata = struct();
% for i = 1:length(markers)
%    markerdata.(char(markers(i))) = trial.(trialname).Trajectories.Labeled.Data(i, :, :);
% end
% threshold = CalculateThreshold(trial, trialname);

%% Determine Foot Based on Label

% %% Locate Heel Strike Times
% % Extract and process heel coordinates
% coords = squeeze(markerdata.(foot)(1, 1:3, :))';
% 
% % Calculate velocity magnitude
% v_mag = vecnorm(diff(coords), 2, 2);
% 
% % Filter velocity data
% [b, a] = butter(m, f_lowpass / (framerate / 2), "low");
% v_mag_filt = filtfilt(b, a, v_mag);
% 
% % Identify heel strike windows
% time_window = round(checktime * framerate); % check 0.1 second intervals
% windows = identify_windows(v_mag_filt, threshold, time_window, upperthreshold);
% 
% % Clean up and remove overlapping windows
% checktime = 0.1; % 0.1 seconds
% windows = clean_windows(windows, checktime, framerate, delaytime);
% 
% % Convert to frames and seconds
% strikes.frames = round(windows(:, 1)); 
% strikes.secs = strikes.frames / framerate;
% 
% % Check with Plot
% figure();
% title("Locate Strikes Graph")
% hold on
% 
% % plot original data
% steps = 1/framerate;
% t = (3*steps/2:steps:trialtime-steps/2)';
% plot(t, v_mag)
% 
% % plot strike points
% scatter(windows(:, 1) / framerate, v_mag(windows(:, 1)))
% 
% xlim([0 10])
% legend(foot, 'window starts')
% hold off

%% Filter EMG Data
% Extract and preprocess EMG data
channel = find(contains(cellstr(EMGdata.(trialname).Channels), label)); %find the channel corresponding to the marker we wanna check (label)
%EMG = EMGdata.Data(channel, :);
specified_EMG = EMGdata.(trialname).Data(channel, :); %collect only the data from that sensor
specified_EMG = specified_EMG(~isnan(specified_EMG)); % remove unused rows

% Bandpass Filter
[b, a] = butter(n, [f_low f_high] / (fs / 2));
EMG_bp = filtfilt(b, a, specified_EMG);

% Rectify
EMG_abs = abs(EMG_bp);

% Lowpass Filter
[b, a] = butter(m, f_lowpass/(fs/2));
EMG_filtered = filtfilt(b, a, EMG_abs);

plot((1:50000)/fs,specified_EMG(1:50000))
figure
plot((1:50000)/fs,EMG_filtered(1:50000))
%% Segment EMG Data
% Segment EMG data based on strike times    
strikes_secs = MCdata.data.strikes.(trialname).(side).secs;
num_segments = length(strikes_secs) - 1;
maxlength = fs * maxsteptime;
segments = cell(num_segments, 1);

for i = 1:num_segments
    %start_idx = round(strikes_secs(i,1) * length(specified_EMG) / trialtime); %Get the index (from the second. Same result than with the next line. 
    start_idx = round(strikes_secs(i,1) * fs); %Get the index (from the second)
    %end_idx = round(strikes_secs((i+1),1) * length(specified_EMG) / trialtime); %i+1 and not the other column because we are only taking into consideration one leg and its cycle
    end_idx = round(strikes_secs((i+1),1) * fs); %i+1 and not the other column because we are only taking into consideration one leg and its cycle
    steptime = end_idx - start_idx; %how much time the step takes
    idx(i,1) = start_idx;
    idx(i,2) = end_idx;
    idx(i,3) = steptime;
    if end_idx <= length(EMG_filtered) && steptime < maxlength
        segments{i} = EMG_filtered(start_idx:end_idx); % Store segment in cell array
    else
        segments{i} = [];
    end 
end

% Filter out empty segments
segments = segments(~cellfun(@isempty, segments));
num_segments = length(segments);

% Plot the segments
figure;
hold on;
for i = 1:num_segments
    segment_length = length(segments{i});
    plot((0:segment_length-1) / segment_length * 100, segments{i}, 'Color', "#c0c0c0");
end

% Compute and plot average EMG amplitude
max_length = max(cellfun(@length, segments)); % find longest segment
avg = compute_average(segments);
plot((0:max_length-1) / max_length * 100, avg, 'Color', "#003366", 'LineWidth', 2); %Plot the average

xlabel('Gait Cycle (%)');
ylabel('EMG Amplitude (mV)');
title([trialname, ' ', label]);
hold off;

%% Helper Functions

function windows = identify_windows(v_mag_filt, threshold, time_window, upperthreshold)
    % Identify strike windows based on velocity magnitude and threshold
    num_samples = length(v_mag_filt);
    windows = NaN(num_samples, 2);
    i = 1;
    j = 1;

    while i < num_samples - time_window
        if all(v_mag_filt(i:i+time_window) < threshold)
            windows(j, 1) = i; % start time
            while v_mag_filt(i + time_window) < threshold && i < num_samples - time_window
                i = i + 1;
            end
            windows(j, 2) = i + time_window - 1; % end time
            j = j + 1;
        end
        i = i + 1;
    end

    windows = windows(~isnan(windows(:, 1)), :); % remove unused rows

    % remove windows where threshold is not exceeded since previous
    window_count = 1;
    for i = 1:length(windows) - 1
        % get segment
        segment = v_mag_filt(windows(window_count, 1):windows(window_count + 1, 1));
        if any(segment > upperthreshold)
            window_count = window_count + 1;
        else
            windows(window_count + 1, 1) = NaN;
            windows(window_count + 1, 2) = NaN;
            windows = windows(~isnan(windows(:, 1)), :); % remove row
        end
    end
end

function windows = clean_windows(windows, checktime, framerate, delaytime)
    % Clean up strike windows to remove overlaps and early windows
    check_samples = round(checktime * framerate);
    valid_indices = [true; diff(windows(:, 1)) > check_samples];
    windows = windows(valid_indices, :);
    windows = windows(windows(:, 1) >= delaytime * framerate, :);
end

function avg = compute_average(segments)
    % Compute average EMG amplitude across segments
    max_length = max(cellfun(@length, segments));
    steps = linspace(0, 100, max_length);
    num_segments = length(segments);
    avg = zeros(1, max_length);

    for i = 1:num_segments
        segment = segments{i};
        segment_length = length(segment);
        percent = linspace(0, 100, segment_length);
        interpolated = interp1(percent, segment, steps, 'linear', 'extrap');
        avg = avg + interpolated;
    end

    avg = avg / num_segments;
end