%% Joint Angle Plotting by Gait Cycle
% This script takes trials and .mot files and plots joint angles,
% segmenting by and averaging across the gait cycle
clc
clear
close all

%% Variables
filepath = "C:\Ability Lab\6_25_2024\";
motfilepath = "C:\Ability Lab\6_25_2024\OpenSim Files\Sarah2\IK\";
trialtitle = "Sarah's Trials -";
numtrials = 9;
angles = ["ankle_angle" "hip_flexion" "knee_angle"];
joints = ["Ankle" "Hip" "Knee"];

% Choose which graphs to generate
separated = 1; % 0 = all gait cycles plotted, 1 = mean plotted w std shading

%% Load Trials and Get Heel Strikes and Thresholds for all Trials
trials = [];
for i = 1:numtrials

    % file and trial name
    trialname = sprintf('Trial%04d', i);
    filename = filepath + trialname;
    motfilename = motfilepath + trialname;

    % load files
    trial = load([num2str(filename) '.mat']);
    mot = importdata(motfilename + num2str('.mot'));

    % add trial, threshold, and heel strikes to trials struct
    trials.(trialname).trial = trial;
    trials.(trialname).mot = mot;
    trials.(trialname).threshold = CalculateThreshold(trial, trialname);
    trials.(trialname).strikes = LocateHeelStrikes(trials.(trialname).trial, trialname, trials.(trialname).threshold);
end

%% Create angle labels 
labels = [strcat(angles, '_r'), strcat(angles, '_l')];
jointlabels = strcat(' ', [joints, joints], ' Angle');
colors = colormap(turbo(numtrials)); % Define a set of colors

%% Plot Trials Separately (either with individual gait cycles or std)
% Loop through labels 
for i = 1:length(labels) % iterate through muscles
    label = string(labels(i));
    jointlabel = string(jointlabels(i));
    side = "right";
    if endsWith(label, 'l')
        side = "left";
    end

    ax = gobjects(3, ceil(numtrials/3)); % initialize subplot array
    figurecount = 0; % initialize figure counter

    for j = 1:numtrials % plot trials
        if rem((j - 1), 9) == 0 % new figure for every 9 trials
            figure()
            hold on
            figurecount = figurecount + 1;
        end

        trialname = sprintf('Trial%04d', j);
        strikes = trials.(trialname).strikes;
        mot = trials.(trialname).mot;
        
        % find subplot row & column
        row = mod(floor((j - 1) / 3), 3) + 1;
        column = rem(j - 1, 3) + 1;
        
        % plot 
        ax(row, column + (figurecount - 1)*3) = subplot(3, 3, j - (figurecount - 1) * 9);
        if separated == 0 
            AnglePlot(trialname, strikes, label, side, mot);
        else
            AvgPlot(trialname, strikes, label, side, mot);
        end
    end

    % Figure Formatting
    sgtitle(trialtitle + replace(side, {'r', 'l'}, {' R', ' L'}) + jointlabel);  % Set title for entire figure
    
    % Standardize y-limits across subplots
    allYLim = cell2mat(arrayfun(@(x) ylim(x), ax, 'UniformOutput', false));
    maxYLim = max(max(allYLim(:, 2:2:numtrials*2/3)));
    minYLim = min(min(allYLim(:, 1:2:numtrials*2/3 - 1)));

    for j = 1:numel(ax)
        ylim(ax(j), [minYLim, maxYLim]);
    end
end


%% Plot Trial Averages on Top of each other
for i = 1:length(labels) % iterate through muscles
    label = string(labels(1, i));
    jointlabel = string(jointlabels(1, i));
    side = "right";
    if endsWith(label, 'l')
        side = "left";
    end

    % Get averages for each trial
    avgs = cell(numtrials, 1);
    legendlabels = strings(numtrials, 1);

    for j = 1:numtrials 
        % isolate data
        trialname = sprintf('Trial%04d', j);
        strikes = trials.(trialname).strikes;
        mot = trials.(trialname).mot;
        
        % run helper function and store avg array in avgs
        avgs{j} = GetAvg(strikes, label, side, mot);

        % create plot label for trial
        legendlabels(j) = sprintf('Trial%04d', j); 
    end

    % Create steps based on longest avg segment (fit same x axis)
    max_length = max(cellfun(@length, avgs)); % find longest segment
    steps = (0:max_length-1) / max_length * 100;    
    
    % create figure
    fig = figure(i + length(labels));

    % create empty matrix for interpolated segments 
    interpolated = zeros(max_length, numtrials);

    % Interpolate and Plot Average Segments 
    hold on
    for j = 1:numtrials

        % Create stepping vector
        avg_length = length(avgs{j});
        percent = (0:avg_length-1) / avg_length * 100; % calc percent steps
       
        % Interpolate amplitude to match steps
        interpolated(:, j) = interp1(percent, avgs{j}, steps, 'linear', 'extrap');

        % plot segment
        plot(steps, interpolated(:, j)', "Color", colors(j, :), "LineWidth", 1);            
    end

    % Figure Formatting
    legend(legendlabels, "Location","best");
    sgtitle(trialtitle + replace(side, {'r', 'l'}, {' R', ' L'}) + jointlabel);  % title for entire figure
    hold off
end

%% Plotting Function - Plots all gait cycles with average
function[] = AnglePlot(trialname, strikes, label, side, mot)
    % Get Angle Data
    data = mot.data(:, find(mot.colheaders == label)) * (180/pi);

    %% Segment Angle Data
    % Initialize segments to store segmented data
    num_segments = length(strikes.(side).secs) - 1; % one less than starts to account for partial gait cycle at end
    segments = cell(num_segments, 1); 
    
    
    % Segment the data based on time stamps
    for i = 1:num_segments
        start_idx = strikes.(side).frames(i); % start index
        end_idx = strikes.(side).frames(i + 1); % end index
        segments{i} = data(start_idx:end_idx); 
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
        
        % Check if segment should be included 
        if any(segments{i} > threshold2)
            continue; % Skip this segment
        end
        
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
    ylabel('Joint Angle (rad)');
    axis("tight")
    title(trialname);
    hold off
end


%% Plotting Function - Plots average gait cycle with standard deviation shading
function[] = AvgPlot(trialname, strikes, label, side, mot)
    % Get Angle Data
    data = mot.data(:, find(mot.colheaders == label)) * (180/pi);

    %% Segment Angle Data
    % Initialize segments to store segmented data
    num_segments = length(strikes.(side).secs) - 1; % one less than starts to account for partial gait cycle at end
    segments = cell(num_segments, 1); 
    
    % Segment the data based on time stamps
    for i = 1:num_segments
        start_idx = strikes.(side).frames(i); % start index 
        end_idx = strikes.(side).frames(i + 1); % end index 
        segments{i} = data(start_idx:end_idx); 
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
        
        % Check if segment should be included 
        if any(segments{i} > threshold2)
            continue; % Skip this segment
        end
        
        % Interpolate amplitude to match steps
        interpolated(:, i) = interp1(percent, segments{i}, steps, 'linear', 'extrap');
                   
        % Increment count of valid segments
        valid_segments_count = valid_segments_count + 1;
    end

    % trim interpolated for unused segments
    interpolated = interpolated(:, 1:valid_segments_count);
    
    % Create average and std arrays
    avg = mean(interpolated, 2)';
    deviation = std(interpolated, 0, 2)';

    % Plot shaded Standard Deviation
    fill([steps, flip(steps)], [avg+deviation, flip(avg-deviation)], [0.753 0.753 0.753], 'EdgeColor','none')

    % Plot Average
    plot(steps, avg, 'color', "#003366", 'LineWidth', 2);
    
    % Plot Appearance 
    xlabel('Gait Cycle (%)');
    ylabel('Joint Angle (rad)');
    axis("tight")
    title(trialname);
    hold off
end

%% Function - Returns Average for inmputted trial
function[avg] = GetAvg(strikes, label, side, mot)
    % Get Angle Data
    data = mot.data(:, find(mot.colheaders == label)) * (180/pi);

    %% Segment Angle Data
    % Initialize segments to store segmented data
    num_segments = length(strikes.(side).secs) - 1; % one less than starts to account for partial gait cycle at end
    segments = cell(num_segments, 1); 
    
    % Segment the data based on time stamps
    for i = 1:num_segments
        start_idx = strikes.(side).frames(i); % start index
        end_idx = strikes.(side).frames(i + 1); % end index
        segments{i} = data(start_idx:end_idx);
    end

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
        
        % Check if segment should be included 
        if any(segments{i} > threshold2)
            continue; % Skip this segment
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