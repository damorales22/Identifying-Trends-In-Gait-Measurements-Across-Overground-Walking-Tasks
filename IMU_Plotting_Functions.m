%% Plotting IMU data - Cleaned up, with functions
clc
clear
close all

%% Create Variables
filepath = "/Users/sarahli/Documents/Ability Lab/6_25_Data/"; % change
date = 'june_25';
numtrials = 9;

IMUcodes = ["742", "722", "511", "589", "568", "598", "668"];
IMUlocations = ["Right_Foot", "Left_Foot", "Right_Thigh", "Left_Thigh", "Right_Shank", "Left_Shank", "Pelvis"];

trialspeeds = [2 2 2 0 0 0 1 1 1]; % input trial speedprompts (0=slow, 1=medium, 2=fast)
trialaccuracies = [0 1 2 0 2 1 0 2 1]; % input trial accuracy prompts (0=none, 1=medium, 2=high)

trial = 1; % choose trial number for Figure 3

% crop_start = 700;
% crop_end = 1000;

% lowpass filter inputs
f_lowpass = 6;
m = 2; % fourth order // filtfilt doubles the order

totalfigures = 1; % count for graphs

%% Load IMU data from all trials
% Initializes struct for all IMU data
IMUdatastruct = struct();

% Initializes variables to determine whether a trial is missing
cutout = [];
numcutout= 1;

% Iterates through all IMUs for all trials
for i = 1:length(IMUcodes)
    sensor = IMUcodes(i);
    for j = 1:numtrials % iterate through trials

            % Check that trial hadn't previously cut out
            if ismember(j, cutout)
                continue
            end

            % Sets filename and path
            filename = [date '_trial_' num2str(j) '_00340' num2str(sensor) '.txt']; % change
            fullpath = join([filepath 'IMU Data/' filename], ""); % change
            
            % Try to open and process the file
            try
                fileID = fopen(fullpath, 'r');
                if fileID == -1
                    error('File not found'); % Trigger error if file not found
                end
           
            % Keeps track of which trials not found and skips it
            catch ME
                if ~ismember(j, cutout)
                    cutout(numcutout) = j;
                    numcutout = numcutout + 1;
                end
                continue
            end
    
            % Gets line containing frequecy
            fgetl(fileID); % Read first line + ignore
            secondLine = fgetl(fileID); % Read the second line
            
            % Gets frequency from line
            floats = regexp(secondLine, '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?', ...
                'match'); % Isolate numeric digits from the second line
            freq = str2double(floats); % Convert to num
            
            % Save freq/IMU/trial
            IMUdatastruct.(IMUlocations(i)).(['Trial' num2str(j)]).freq = freq;
            
            % Load IMU data
            opts = detectImportOptions(fullpath);
            opts.VariableNamingRule = 'preserve';
            IMUdatastruct.(IMUlocations(i)).(['Trial' num2str(j)]).data = readtable(fullpath, opts);
            
    end
end

disp(["trials where the file wasn't found: " cutout])

%% Load Trial Data to Get Heelstrikes
trialstruct = struct();
for i = 1:numtrials
    trialname = sprintf('Trial%04d', i);
    trialstruct.(trialname) = load([num2str(filepath) char(trialname) '.mat']); % load trial
    trialstruct.(trialname).strikes = LocateHeelStrikes(trialstruct.(trialname), trialname, ...
            CalculateThreshold(trialstruct.(trialname), trialname)); % get heel strikes
end

%% Create Labels for Graphs
labels = {};
for i = 1:9
    % Speeds
    if trialspeeds(i) == 0
        labels{i} = "Slow Speed";
    elseif trialspeeds(i) == 1
        labels{i} = "Medium Speed";
    elseif trialspeeds(i) == 2
        labels{i} = "Fast Speed";
    end
    % Accuracies
    if trialaccuracies(i) == 0
        labels{i} = labels{i} + " Low Acc.";
    elseif trialaccuracies(i) == 1
        labels{i} = labels{i} + " Med Acc.";
    elseif trialaccuracies(i) == 2
        labels{i} = labels{i} + " High Acc.";
    end
end

%% Run functions
% Changes these based on plotting preferances
type = "gyro";
selectsensors = ["Right_Foot", "Left_Thigh"];
selecttrials = [1 2 3];

totalfigures = PlotAllAverages(IMUcodes, IMUlocations, IMUdatastruct, ...
    totalfigures, numtrials, numcutout, type, trialstruct, m, f_lowpass, labels, cutout, selectsensors, selecttrials);

totalfigures = PlotAvgsAndGaits(IMUcodes, IMUlocations, IMUdatastruct, ...
    totalfigures, type, trialstruct, m, f_lowpass);

totalfigures = PlotFullTrial(IMUcodes, IMUlocations, IMUdatastruct, ...
    totalfigures, type, trialstruct, m, f_lowpass, labels, cutout);

%% Plot All Trial Averages Per Sensor
function[totalfigures] = PlotAllAverages(IMUcodes, IMUlocations, IMUdatastruct, ...
    totalfigures, numtrials, numcutout, type, trialstruct, m, f_lowpass, ...
    labels, cutout, selectsensors, selecttrials)
    % Checks accel or gyro
    if type == "accel"
        startdim = 2;
        enddim = 4;
        plot_title = " Acceleration";
        plot_y_label = 'Acceleration (m/s^2)';
    elseif type == "gyro"
        startdim = 8;
        enddim = 10;
        plot_title = " Gyroscope";
        plot_y_label = 'Angular Rate (rad/s)';
    end

    % iterate through IMU sensors
    for i = 1:length(IMUcodes)
        % Check whether desired trial
        if ~ismember(IMUlocations(i), selectsensors)
            continue
        end
        
        trials = fieldnames(IMUdatastruct.(IMUlocations(i)));
        figure(totalfigures)
        totalfigures = totalfigures + 1;
        hold on    
        % subplot(4, 2, i) % 2x4, one plot for each sensor
        % hold on
        % sensor = IMUcodes(i); 
    
        % determines which heelstrikes to use
        if contains(IMUlocations(i), "Right")
            side = "right";
        else
            side = "left";
        end
    
        % Initialize array to hold avgs
        avgs = cell(numtrials-numcutout, 1);
    
        for j = 1:length(trials) % iterate through trials
            % Check to see whether it's a selected trial
            trialnum = sscanf(trials{j}, 'Trial%d');
            if ~ismember(trialnum, selecttrials)
                continue
            end
            
            % trial variables
            trialname = sprintf('Trial%04d', trialnum);
            trialfr = trialstruct.(trialname).(trialname).FrameRate;
    
            % extract data from struct
            IMUfr = IMUdatastruct.(IMUlocations(i)).(trials{j}).freq;
            IMUdata = IMUdatastruct.(IMUlocations(i)).(trials{j}).data;
            % Column 1: Timestamp; Columns 2-4: Accelerometers X, Y, Z
    
            % demean data
            if width(IMUdata) < 21 % if IMU hasn't cut out
                continue
            end
            means = mean(IMUdata(:, 1:21)); % calculate mean values
    
            for k = startdim:enddim % iterate through dimensions (x, y, z) and subtract mean
               IMUdata(:, k) = IMUdata(:, k) - means(1, k);
            end
    
            % calculate vector magnitudes and convert to array
            IMUdata_2 = IMUdata(:, startdim:enddim).^2;
            ydata = table2array(sqrt(sum(IMUdata_2, 2)));
    
            % lowpass butterworth filter
            [b, a] = butter(m, f_lowpass/(IMUfr/2), "low");
            ydata = filtfilt(b, a, ydata);
    
            %% Segmentation
            % run helper function and store avg array in avgs
            avgs{j} = TrialAvgPlot(ydata, trialstruct.(trialname).strikes, side, IMUfr, trialfr); % need to interpolate strike times to match freq
        end
    
        % Create steps based on longest avg segment (fit same x axis)
        max_length = max(cellfun(@length, avgs)); % find longest segment
        steps = (0:max_length-1) / max_length * 100;    
    
        % create empty matrix for interpolated segments 
        interpolated = zeros(max_length, numtrials);
    
        % Interpolate and Plot Average Segments 
        for j = 1:numtrials-numcutout
    
            % Create stepping vector
            avg_length = length(avgs{j});
            percent = (0:avg_length-1) / avg_length * 100; % calc percent steps
    
            % Interpolate amplitude to match steps
            if isempty(avgs{j})
                continue
            end
    
            interpolated(:, j) = interp1(percent, avgs{j}, steps, 'linear', 'extrap');
    
            % plot segment
            plot(steps, interpolated(:, j)', "LineWidth", 1);            
        end
    
        % format plot
        xlim tight
        title(replace(num2str(IMUlocations(i)), "_", " ") + plot_title);  % Set title for entire figure
        xlabel ('Gait Cycle (%)')
        ylabel(plot_y_label)
        leftover_labels = labels;
        leftover_labels(cutout) = [];
        legend(leftover_labels)
        hold off
    end
end

%% Plot Averages Ontop of All Gait Cycles
function[totalfigures] = PlotAvgsAndGaits(IMUcodes, IMUlocations, IMUdatastruct, ...
    totalfigures, type, trialstruct, m, f_lowpass)
    % Checks accel or gyro
        if type == "accel"
            startdim = 2;
            enddim = 4;
            plot_title = " Acceleration";
            plot_y_label = 'Acceleration (m/s^2)';
        elseif type == "gyro"
            startdim = 8;
            enddim = 10;
            plot_title = " Gyroscope";
            plot_y_label = 'Angular Rate (rad/s)';
        end    
    
    % For subplots
    figure(totalfigures)
    totalfigures = totalfigures + 1;
    hold on 
    
    for i = 1:length(IMUcodes) % iterate through IMU sensors
        trials = fieldnames(IMUdatastruct.(IMUlocations(i)));
        
        % % For full figures
        % figure(totalfigures)
        % totalfigures = totalfigures + 1;
        % hold on 
        
        subplot(4, 2, i) % 2x4, one plot for each sensor
    
        % determines which heelstrikes to use
        if contains(IMUlocations(i), "Right")
            side = "right";
        else
            side = "left";
        end
        
        % iterate through trials
        for j = 1:1 %length(trials)
    
            % trial variables
            trialname = sprintf('Trial%04d', str2num(trials{j}(end)));
            trialfr = trialstruct.(trialname).(trialname).FrameRate;
    
            % extract data from struct
            IMUfr = IMUdatastruct.(IMUlocations(i)).(trials{j}).freq;
            IMUdata = IMUdatastruct.(IMUlocations(i)).(trials{j}).data;
            % Column 1: Timestamp; Columns 2-4: Accelerometers X, Y, Z
    
            % demean data
            if width(IMUdata) < 21 % if IMU hasn't cut out
                continue
            end
            means = mean(IMUdata(:, 1:21)); % calculate mean values
    
            for k = startdim:enddim % iterate through dimensions (x, y, z) and subtract mean
               IMUdata(:, k) = IMUdata(:, k) - means(1, k);
            end
    
            % calculate vector magnitudes and convert to array
            IMUdata_2 = IMUdata(:, startdim:enddim).^2;
            ydata = table2array(sqrt(sum(IMUdata_2, 2)));
    
            % lowpass butterworth filter
            [b, a] = butter(m, f_lowpass/(IMUfr/2), "low");
            ydata = filtfilt(b, a, ydata);
    
            % Segmentation + Plotting
            IMUAllAvg(ydata, trialstruct.(trialname).strikes, side, IMUfr, trialfr); % need to interpolate strike times to match freq
        end
    
        % format plot
        xlim tight
        title(replace(num2str(IMUlocations(i)), "_", " ") + plot_title);  % Set title for entire figure
        xlabel ('Gait Cycle (%)')
        ylabel(plot_y_label)
        hold off
    end
end

%% Plot Full Trial
function[totalfigures] = PlotFullTrial(IMUcodes, IMUlocations, IMUdatastruct, ...
    totalfigures, type, trialstruct, m, f_lowpass, labels, cutout)
    % Checks accel or gyro
    if type == "accel"
        startdim = 2;
        enddim = 4;
        plot_title = " Acceleration";
        plot_y_label = 'Acceleration (m/s^2)';
    elseif type == "gyro"
        startdim = 8;
        enddim = 10;
        plot_title = " Gyroscope";
        plot_y_label = 'Angular Rate (rad/s)';
    end    
    
    % iterate through IMU sensors    
    for i = length(IMUcodes)-1:length(IMUcodes)-1 
        trials = fieldnames(IMUdatastruct.(IMUlocations(i)));
        figure(totalfigures)
        totalfigures = totalfigures + 1;
        hold on  
        % subplot(4, 2, i) % 2x4, one plot for each sensor
        % hold on
        
        % iterate through trials
        for j = 1:1% length(trials) 
    
            % trial variables
            trialname = sprintf('Trial%04d', str2num(trials{j}(end)));
            trialfr = trialstruct.(trialname).(trialname).FrameRate;
    
            % extract data from struct
            IMUfr = IMUdatastruct.(IMUlocations(i)).(trials{j}).freq;
            IMUdata = IMUdatastruct.(IMUlocations(i)).(trials{j}).data;
            % Column 1: Timestamp; Columns 2-4: Accelerometers X, Y, Z
    
            % demean data
            if width(IMUdata) < 21 % if IMU hasn't cut out
                continue
            end
            means = mean(IMUdata(:, 1:21)); % calculate mean values
    
            for k = startdim:enddim % iterate through dimensions (x, y, z) and subtract mean
               IMUdata(:, k) = IMUdata(:, k) - means(1, k);
            end
    
            % calculate vector magnitudes and convert to array
            IMUdata_2 = IMUdata(:, startdim:enddim).^2;
            ydata = table2array(sqrt(sum(IMUdata_2, 2)));
    
            % lowpass butterworth filter
            [b, a] = butter(m, f_lowpass/(IMUfr/2), "low");
            ydata = filtfilt(b, a, ydata);
    
            % Create Time Vector
            time = 1:length(ydata);
    
            % plot accelerometer data
            plot(time, ydata)
    
        end
    
        % format plot
        xlim([1,720])
        % ylim([-20, 50])
        title(replace(num2str(IMUlocations(i)), "_", " ") + plot_title)
        xlabel('Frames')
        ylabel(plot_y_label)
        leftover_labels = labels;
        leftover_labels(cutout) = [];
        legend(leftover_labels)
        hold off
    end
end

%% Helper Function for Plotting All Averages
function[avg] = TrialAvgPlot(accel, strikes, side, IMUfr, trialfr)
    % Initialize segments to store segmented data
    num_segments = length(strikes.(side).secs) - 1; % one less than starts to account for partial gait cycle at end
    segments = cell(num_segments, 1); 

    % Segment the data based on time stamps
    for i = 1:num_segments
        start_idx = round(strikes.(side).frames(i)/trialfr*IMUfr); % start index
        end_idx = round(strikes.(side).frames(i + 1)/trialfr*IMUfr); % end index
        if end_idx <= length(accel)
            segments{i} = accel(start_idx:end_idx);
        end
    end

    hold on

    % Create steps (based on longest segment)
    max_length = max(cellfun(@length, segments)); % find longest segment
    steps = (0:max_length-1) / max_length * 100; 
    
    % Initialize a counter for valid segments
    valid_segments_count = 0;

    % Calculate average of all segment amplitudes
    total_amp = 0;  % Initialize an array to accumulate all segment amplitudes
    for i = 1:num_segments
        tempmax = max(segments{i});
        total_amp = total_amp + tempmax;
    end
    average_amplitude = total_amp/num_segments;  % Average amplitude of all segments
    
    % Define the threshold based on the average amplitude
    threshold2 = average_amplitude*3;

    % create empty matrix for interpolated segments 
    interpolated = zeros(max_length, num_segments);

    % Interpolate Segments
    for i = 1:num_segments

        % Create stepping vector
        segment_length = length(segments{i});
        percent = (0:segment_length-1) / segment_length * 100; % calc percent steps
        
        % Skip empty segments
        if isempty(segments{i})
            continue
        end

        % Interpolate amplitude to match steps
        interpolated(:, i) = interp1(percent, segments{i}, steps, 'linear', 'extrap');
                   
        % Increment count of valid segments
        valid_segments_count = valid_segments_count + 1;
    end

    % trim interpolated for unused segments
    interpolated = interpolated(:, 1:valid_segments_count);
    
    % Create average and std arrays
    avg = mean(interpolated, 2);
end

%% Helper Function - Plots all gait cycles with average
function[] = IMUAllAvg(accel, strikes, side, IMUfr, trialfr)
    % Segment Angle Data
    % Initialize segments to store segmented data
    num_segments = length(strikes.(side).secs) - 2; % one less than starts to account for partial gait cycle at end
    segments = cell(num_segments, 1); 
    
    
    
    % Segment the data based on time stamps
    for i = 1:num_segments
        start_idx = round(strikes.(side).frames(i)*IMUfr/trialfr); % start index
        end_idx = round(strikes.(side).frames(i + 1)*IMUfr/trialfr); % end index
        if end_idx <= length(accel)
            segments{i} = accel(start_idx:end_idx); 
        end
    end
    hold on

    % Create steps (based on longest segment)
    max_length = max(cellfun(@length, segments)); % find longest segment
    steps = (0:max_length-1) / max_length * 100; 

    % Initialize arrays for average y data
    avg = zeros(size(steps));
    
    % Initialize a counter for valid segments
    valid_segments_count = 0;
    
    % Calculate average of all segment amplitudes
    total_amp = 0;  % Initialize an array to accumulate all segment amplitudes
    for i = 1:num_segments
        tempmax = max(segments{i});
        total_amp = total_amp + tempmax;
    end
    average_amplitude = total_amp/num_segments;  % Average amplitude of all segments
    
    % Define the threshold based on the average amplitude
    threshold2 = average_amplitude*3;

    % Plot the segments and calculate average y data
    for i = 1:num_segments

        % Create stepping vector
        segment_length = length(segments{i});
        percent = (0:segment_length-1) / segment_length * 100; % calc percent steps
        
        % Plot the segment
        plot(percent, segments{i}, 'Color', "#c0c0c0");
          
        if isempty(segments{i})
            continue
        end

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
    ylabel('Joint Angle (rad)');
    axis("tight")
    hold off
end