%% EMG Plotting by gait cycle across trials
% Takes trials and plots the EMG data for each muscle, segmented by gait
% cycle with an overlaid average
clc 
clear
close all

%% Variables
filepath = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\BMH10(September)\MC\";
EMGpath = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\BMH10(September)\EMG\";
trialtitle = "BMH05's Trials ";
numtrials = 27;
maxsteptime = 2; % set maximum step length for segment

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
EMGname = EMGpath + trialname;
trial = load([num2str(filename) '.mat']); % load first trial
EMGdata = load([num2str(EMGname) '.mat']);
labels = cellstr(EMGdata.Channels);
fs = EMGdata.Fs(1);

%% Load Trials and Get Heel Strikes and Thresholds for all Trials
trials = struct();
for i = 1:numtrials
    trialname = sprintf('Trial%04d', i);
    filename = filepath + trialname;
    EMGname = EMGpath + trialname;

    % load file
    trial = load([num2str(filename) '.mat']);

    % Add trial, threshold, and heel strikes to trials struct
    trials.(trialname).trial = trial;
    trials.(trialname).threshold = CalculateThreshold(trial, trialname);
    trials.(trialname).strikes = LocateHeelStrikes(trials.(trialname).trial, trialname, trials.(trialname).threshold);
    trials.(trialname).EMG = load([num2str(EMGname) '.mat']);
    disp([trialname ' loaded'])
end

%% Loop through labels 
for i = 1:length(labels) % iterate through muscles
    label = labels(i);
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
        EMGdata = trials.(trialname).EMG;

        % filter and segment EMG data
        [trials.(trialname).segments, trials.(trialname).average] = EMGProcess(EMGdata, strikes, fs, f_high, f_low, f_lowpass, n, m, i, side, maxsteptime);

        % find subplot row & column
        row = mod(floor((j - 1) / 3), 3) + 1;
        column = rem(j - 1, 3) + 1;

        % plot 
        ax(row, column + (figurecount - 1)*3) = subplot(3, 3, j - (figurecount - 1)*9);
        EMGPlot(trialname, trials.(trialname).segments);

        % Figure Formatting
        sgtitle(trialtitle + strrep(label, '_EMG 1', ''));  % Set title for entire figure
    end

    % Standardize y-limits across subplots
    ax = ax(isgraphics(ax));
    allYLim = cell2mat(arrayfun(@(x) ylim(x), ax, 'UniformOutput', false));
    maxYLim = max(max(allYLim(:, 2:2:end)));
    minYLim = min(min(allYLim(:, 1:2:end)));

    for j = 1:numel(ax)
        ylim(ax(j), [minYLim, maxYLim]);
    end
    disp(string(label) + " plotted")
end


%% EMG Processing Function
function[segments, avg] = EMGProcess(EMGData, strikes, fs, f_high, f_low, f_lowpass, n, m, channel, side, maxsteptime)
    %% Get frame info
    framerate = EMGData.Fs(1);
    numframes = length(EMGData.Data);
    trialtime = numframes/framerate;

    %% Filter EMG data
    % isolate EMG data
    EMG = EMGData.Data(channel, :); % extract EMG data for label
    EMG = EMG(~isnan(EMG(:))); % remove unused rows
    numsamples = length(EMG);
    
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
    maxlength = fs * maxsteptime;
    segments = arrayfun(@(i) segmentEMGData(EMG_filtered, strikes.(side).secs(i), strikes.(side).secs(i+1), numsamples, trialtime, maxlength), 1:num_segments, 'UniformOutput', false);

    % Filter out empty segments
    segments = segments(~cellfun(@isempty, segments));

    % create steps (based on longest segment)
    num_segments = length(segments);
    max_length = max(cellfun(@length, segments)); % find longest segment
    steps = (0:max_length-1) / max_length * 100;

    % Initialize array for average y data and valid segment counter
    avg = zeros(size(steps));

    % Interpolate segments
    for i = 1:num_segments
        % create stepping vector
        segment_length = length(segments{i});
        percent = (0:segment_length-1) / segment_length * 100; % calc percent steps
    
        % Interpolate amplitude to match steps
        segments{i} = interp1(percent, segments{i}, steps, 'linear', 'extrap');

        % Accumulate average y data
        avg = avg + segments{i}; % add to average
    end

    % Calculate average 
    avg = avg ./ num_segments;
end

%% Plotting Function
function[] = EMGPlot(trialname, segments)

    hold on
    % Create steps (based on longest segment)
    num_segments = length(segments);
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

%% Function to Segment EMG Data
function segment = segmentEMGData(EMG_filtered, start_sec, end_sec, numsamples, trialtime, maxlength)
    start_idx = round(start_sec * numsamples / trialtime); % Start index of segment
    end_idx = round(end_sec * numsamples / trialtime); % End index of segment
    steptime = end_idx - start_idx;
    if end_idx <= length(EMG_filtered) && steptime < maxlength
        segment = EMG_filtered(start_idx:end_idx); % Store segment in cell array
    else
        segment = [];
    end
end