%% EMG Plotting by gait cycle across trials
% Takes trials and plots the EMG data for each muscle, segmented by gait
% cycle with an overlaid average
clc 
clear
close all

%% Variables
filepath = "C:\Ability Lab\6_24_2024\";
trialtitle = "Jo's Trials ";
numtrials = 9;

% Bandpass filter inputs
f_low = 30;
f_high = 500;
n = 1; % fourth order // bandpass filter and filtfilt both double the order

% lowpass filter inputs
f_lowpass = 6;
m = 2; % fourth order // filtfilt doubles the order

%% Get MoCo labels and EMG Sample Frequency
trialname = "Trial0001";
filename = filepath + trialname;
trial = load([num2str(filename) '.mat']); % load first trial
labels = trial.(trialname).Analog(4).Labels;
fs = trial.(trialname).Analog(4).Frequency;

%% Load Trials and Get Heel Strikes and Thresholds for all Trials
trials = struct();
for i = 1:numtrials
    trialname = sprintf('Trial%04d', i);
    filename = filepath + trialname;

    % load file
    trial = load([num2str(filename) '.mat']);

    % Add trial, threshold, and heel strikes to trials struct
    trials.(trialname).trial = trial;
    trials.(trialname).threshold = CalculateThreshold(trial, trialname);
    trials.(trialname).strikes = LocateHeelStrikes(trials.(trialname).trial, trialname, trials.(trialname).threshold);
end

%% Loop through labels 
for i = 1:length(labels) % iterate through muscles
    label = labels(1, i);
    side = "right"; % check which leg the EMG is on
    if ~startsWith(label , 'R')
        side = "left";
    end

    ax = gobjects(3, ceil(numtrials/3)); % initialize subplot array
    figurecount = 0; % initialize figure counter

    for j = 1:numtrials % plot trials
        if rem((j - 1), 9) == 0 % new figure every 9 trials
            figure()
            hold on
            figurecount = figurecount + 1;
        end
        trialname = sprintf('Trial%04d', j);
        trial = trials.(trialname).trial;
        strikes = trials.(trialname).strikes;

        % find subplot row & column
        row = mod(floor((j - 1) / 3), 3) + 1;
        column = rem(j - 1, 3) + 1;
        
        % plot 
        ax(row, column + (figurecount - 1)*3) = subplot(3, 3, j - (figurecount - 1)*9);
        EMGPlot(trial, trialname, strikes, fs, f_high, f_low, f_lowpass, n, m, label, side);
    end

    % Figure Formatting
    sgtitle(trialtitle + strrep(label, '_EMG 1', ''));  % Set title for entire figure

    % Standardize y-limits across subplots
    allYLim = cell2mat(arrayfun(@(x) ylim(x), ax, 'UniformOutput', false));
    maxYLim = max(max(allYLim(:, 2:2:numtrials*2/3)));
    minYLim = min(min(allYLim(:, 1:2:numtrials*2/3 - 1)));

    for j = 1:numel(ax)
        ylim(ax(j), [minYLim, maxYLim]);
    end
end


%% Plotting Function
function[] = EMGPlot(trial, trialname, strikes, fs, f_high, f_low, f_lowpass, n, m, label, side)
    %% Get frame info
    framerate = trial.(trialname).FrameRate;
    numframes = trial.(trialname).Frames;
    trialtime = numframes/framerate;

    %% Filter EMG data
    % isolate EMG data
    EMG = trial.(trialname).Analog(4).Data(find(trial.(trialname).Analog(4).Labels == string(label)), :); % extract EMG data for label
    EMG = EMG(~isnan(EMG(:))); % remove unused rows
    numsamples = trial.(trialname).Analog(4).NrOfSamples;
    
    % Bandpass Filter
    [b, a] = butter (n, [f_low f_high]/(fs/2));
    EMG_bp = filtfilt(b, a, EMG); 
    
    % Rectify
    EMG_abs = abs(EMG_bp);
    
    % Lowpass Filter
    [b, a] = butter(m, f_lowpass/(fs/2));
    EMG_filtered = filtfilt(b, a, EMG_abs);
    
    %% Segment EMG data
    num_segments = length(strikes.(side).secs) - 1; % One less to account for partial gait cycle at end
    segments = arrayfun(@(i) segmentEMGData(EMG_filtered, strikes.(side).secs(i), strikes.(side).secs(i+1), numsamples, trialtime), 1:num_segments, 'UniformOutput', false);

    % Filter out empty segments
    segments = segments(~cellfun(@isempty, segments));
    num_segments = length(segments);

    %% Plot
    hold on

    % Create steps (based on longest segment)
    max_length = max(cellfun(@length, segments)); % find longest segment
    steps = (0:max_length-1) / max_length * 100; 


    % Calculate average of all segment amplitudes
    avg_amplitude = mean(cellfun(@max, segments));
    threshold2 = avg_amplitude * 3;

    % Initialize array for average y data and valid segment counter
    avg = zeros(size(steps));
    valid_segments_count = 0; 
   

    % Plot the segments and calculate average y data
    for i = 1:num_segments
                
        % Check if segment should be included
        if any(segments{i} > threshold2)
            continue; % Skip this segment
        end
        
        % create stepping vector
        segment_length = length(segments{i});
        percent = (0:segment_length-1) / segment_length * 100; % calc percent steps
        
        % Plot the segment
        plot(percent, segments{i}, 'Color', "#c0c0c0");
    
        % Interpolate amplitude to match steps
        interpolated = interp1(percent, segments{i}, steps, 'linear', 'extrap');
    
        % Accumulate average y data
        avg = avg + interpolated; % add to average
        
        % Increment count of valid segments
        valid_segments_count = valid_segments_count + 1;
    end
    
    % Calculate average by dividing by the number of contributing segments
    if valid_segments_count > 0
        avg = avg ./ valid_segments_count;
    end

    % Plot the average y data
    plot(steps, avg, 'color', "#003366", 'LineWidth', 2);
    
    % Plot Appearance 
    xlabel('Gait Cycle (%)');
    ylabel('EMG Amplitude (mV)');
    axis("tight")
    title(trialname);
    hold off
end

%% Helper Function to Segment EMG Data
function segment = segmentEMGData(EMG_filtered, start_sec, end_sec, numsamples, trialtime)
    start_idx = round(start_sec * numsamples / trialtime); % Start index of segment
    end_idx = round(end_sec * numsamples / trialtime); % End index of segment
    if end_idx <= length(EMG_filtered)
        segment = EMG_filtered(start_idx:end_idx); % Store segment in cell array
    else
        segment = [];
    end
end