%% EMG Plotting by gait cycle across trials
% Takes trials and plots the EMG data for each muscle, segmented by gait
% cycle with an overlaid average
clc 
clear
close all

%% Variables
emgpath = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\Trial1\EMG\";
mocappath = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\Trial1\MC\";
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
trialname_emg = "Trial0001";
filename_emg = emgpath + trialname_emg;
trial_emg = load([num2str(filename_emg) '.mat']); % load first trial

%labels = trial.(trialname).Analog(4).Labels;
numChannels = size(trial_emg.Channels, 1); % Number of channels (rows)
labels = cell(numChannels, 1); % Cell array to store the first 3 letters
for i = 1:numChannels
    fullLabel = strtrim(trial_emg.Channels(i, :)); % Extract the full label from the row and remove trailing spaces
    labels{i} = fullLabel(1:3); % Extract the first three letters
end

%fs = trial.(trialname).Analog(4).Frequency;
fs = trial_emg.Fs(1,1);

%% Load Trials and Get Heel Strikes and Thresholds for all Trials
trials = struct();
for i = 1:numtrials
    trialname = sprintf('Trial%04d', i);
    %filename = filepath + trialname;
    filename = mocappath + trialname;

    % load file
    trial_mocap = load([num2str(filename) '.mat']);

    % Add trial, threshold, and heel strikes to trials struct
    trials.(trialname).trial = trial_mocap;
    trials.(trialname).threshold = CalculateThreshold(trial_mocap, trialname);
    trials.(trialname).strikes = LocateHeelStrikes(trials.(trialname).trial, trialname, trials.(trialname).threshold);
    trials.(trialname).emg = trial_emg;
end

%% Loop through labels 
for i = 1:length(labels) % iterate through muscles
    label = labels(i, 1);
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
        %trial = trials.(trialname).trial;
        trial = trials.(trialname);
        strikes = trials.(trialname).strikes;

        % find subplot row & column
        row = mod(floor((j - 1) / 3), 3) + 1;
        column = rem(j - 1, 3) + 1;
        
        % plot 
        ax(row, column + (figurecount - 1)*3) = subplot(3, 3, j - (figurecount - 1)*9);
        %EMGPlot(trial, trialname, strikes, fs, f_high, f_low, f_lowpass, n, m, label, side);
        EMGPlot(trial, trialname, strikes, fs, f_high, f_low, f_lowpass, n, m, i, side);
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
    break
end


%% Plotting Function
function[] = EMGPlot(trial, trialname, strikes, fs, f_high, f_low, f_lowpass, n, m, num_label, side)
    %% Get frame info
    % framerate = trial.(trialname).FrameRate;
    % numframes = trial.(trialname).Frames;
    framerate = trial.trial.(trialname).FrameRate;
    numframes = trial.trial.(trialname).Frames;    
    trialtime = numframes/framerate;

    %% Filter EMG data
    % isolate EMG data
    %EMG = trial.(trialname).Analog(4).Data(find(trial.(trialname).Analog(4).Labels == string(label)), :); % extract EMG data for label
    EMG = trial.emg.Data(num_label,:); % extract EMG data for label
    EMG = EMG(~isnan(EMG(:))); % remove unused rows
    %numsamples = trial.(trialname).Analog(4).NrOfSamples;
    num = size(trial.strikes.(side).frames); 
    numsamples = num(1,1)-1; %-1 because of the steps
    
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
        if length(segments{i}) < 2 % If it doesn't exist the second element
            %segments(i) = [];
            %segments = segments(~cellfun('isempty', segments));
            fprintf('Error at segments{%d} in %s\n', i, trialname);
            continue
        end
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