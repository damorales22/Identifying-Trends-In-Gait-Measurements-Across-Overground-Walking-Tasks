%% New LocateSteps
% Uses z position of the marker to calculate the positions that must be used to calculate the steps for each foot.

function[steps] = LocateSteps(trial, trialname, counter, usedecisions, automatic_detection, filepath)

   %% Parameters for filtering and thresholds

    fc = 0.1;                  % Cut-off frequency for the Butterworth filter

    min_value_for_max = 75;    % Minimum value for the maxs
    max_value_flag_max = 120;  % Upper threshold to flag it as a possible max => step
    min_separation_maxs = 50;  % Minimum horizontal separation for maxima; to keep detecting trips
    prominence_threshold = 50; % Maxima prominence. Higher value to be more strict with the detected maxima 
    diff_max_min = 80;         % Maximum vertical difference between the maxima and the absolute minima 30 frames back

    max_value_for_min = 30;    % Maximum value for the mins
    min_separation_mins = 50;  % Minimum horizontal separation for minima
    
    % get markerdata
    try
        markers = trial.(trialname).Trajectories.Labeled.Labels;
    catch %if the internal name of that trial is changed
        folder_names = fieldnames(trial);  % Get all field names in 'trial'
        trialname = folder_names{1};  % Assuming the only folder is the first one
        markers = trial.(trialname).Trajectories.Labeled.Labels;
    end
    markerdata = struct();
    for i = 1:length(markers)
       markerdata.(char(markers(i))) = trial.(num2str(trialname)).Trajectories.Labeled.Data(i, :, :);
    end

    decisions = strings(0, 1);

    if usedecisions == 1
        filepath = fullfile(filepath, 'decisions.txt');
        fileID = fopen(filepath, 'r');
        if fileID == -1
            error('No .txt file found. Introduce the decisions manually');
        end
        decisionsArray = textscan(fileID, '%s');  % Read as strings
        fclose(fileID);
        decisionsString = strjoin(decisionsArray{1});
        decisionsCells = cell(1, length(decisionsString)); % Initialize the cell array to store each letter

        for i = 1:length(decisionsString) % Store each letter in a separate cell in the array
            decisionsArray{i} = decisionsString(i);
        end
    end

    %% ASSUMPTION: the position when each foot (MTH5 or LMH5) is on the ground is the middle distance between the peaks of the heel markers when representing their z coordinate
    % Since the don't really care about the timeframe but about the position, this method is okay.

    % Extract right and left foot z-coordinates
    right_position = squeeze(markerdata.RHEEL(1,3,:)); % ONLY z coordinate; BEFORE: MTH5 
    left_position = squeeze(markerdata.LHEEL(1,3,:));  % ONLY z coordinate
    
    % Filter parameters
    [b, a] = butter(4, fc, 'low');  % 4th order Butterworth filter
    
    % Filter the right and left foot signals
    right_filtered = filtfilt(b, a, right_position);
    left_filtered = filtfilt(b, a, left_position);
    
    % Find absolute minimum for shifting both signals to start from zero
    min_right = min(right_filtered);
    min_left = min(left_filtered);
    
    % Shift both signals
    right_shifted = right_filtered - min_right;
    left_shifted = left_filtered - min_left;
    
    % Calculate absolute maximums
    max_right = max(right_shifted);
    max_left = max(left_shifted);
    
    % Find local maxima and minima for both shifted signals
    [local_maxima_right, max_locs_right] = findpeaks(right_shifted);
    [local_minima_right, min_locs_right] = findpeaks(-right_shifted); % Invert to find minima
    local_minima_right = -local_minima_right; % Return minima to original values
    
    [local_maxima_left, max_locs_left] = findpeaks(left_shifted);
    [local_minima_left, min_locs_left] = findpeaks(-left_shifted); % Invert to find minima
    local_minima_left = -local_minima_left; % Return minima to original values
    
    % Filter maxima above the threshold
    high_max_locs_right = max_locs_right(local_maxima_right > min_value_for_max);
    high_maxima_right = local_maxima_right(local_maxima_right > min_value_for_max);
    high_max_locs_left = max_locs_left(local_maxima_left > min_value_for_max);
    high_maxima_left = local_maxima_left(local_maxima_left > min_value_for_max);
    
    check_high_max_locs_right = max_locs_right(local_maxima_right <= min_value_for_max);
    check_high_maxima_right = local_maxima_right(local_maxima_right <= min_value_for_max);
    check_high_max_locs_left = max_locs_left(local_maxima_left <= min_value_for_max);
    check_high_maxima_left = local_maxima_left(local_maxima_left <= min_value_for_max);
    
    check_right_maxs(:,1) = check_high_max_locs_right;
    check_right_maxs(:,2) = check_high_maxima_right;
    check_left_maxs(:,1) = check_high_max_locs_left; 
    check_left_maxs(:,2) = check_high_maxima_left;
    
    % Store right maxima
    right_maxs(:,1) = high_max_locs_right; % X locations of right maxima
    right_maxs(:,2) = high_maxima_right;   % Values of right maxima
    
    % Filter minima below the threshold and ensure a separation of at least 100 units
    right_mins_filtered = min_locs_right(local_minima_right < max_value_for_min);
    right_minima_filtered = local_minima_right(local_minima_right < max_value_for_min);
    
    % Enforce separation of at least min_separation_mins for right foot minima
    valid_right_mins = [];
    valid_right_min_values = [];
    for i = 1:length(right_mins_filtered)
        if isempty(valid_right_mins) || abs(right_mins_filtered(i) - valid_right_mins(end)) >= min_separation_mins
            valid_right_mins = [valid_right_mins; right_mins_filtered(i)];
            valid_right_min_values = [valid_right_min_values; right_minima_filtered(i)];
        end
    end
    right_mins(:,1) = valid_right_mins;
    right_mins(:,2) = valid_right_min_values;
    
    % Store left maxima
    left_maxs(:,1) = high_max_locs_left; % X locations of left maxima
    left_maxs(:,2) = high_maxima_left;   % Values of left maxima
    
    % Filter minima below the threshold and ensure a separation of at least 100 units
    left_mins_filtered = min_locs_left(local_minima_left < max_value_for_min);
    left_minima_filtered = local_minima_left(local_minima_left < max_value_for_min);
    
    % Enforce separation of at least min_separation_mins for left foot minima
    valid_left_mins = [];
    valid_left_min_values = [];
    for i = 1:length(left_mins_filtered)
        if isempty(valid_left_mins) || abs(left_mins_filtered(i) - valid_left_mins(end)) >= min_separation_mins
            valid_left_mins = [valid_left_mins; left_mins_filtered(i)];
            valid_left_min_values = [valid_left_min_values; left_minima_filtered(i)];
        end
    end
    
    left_mins(:,1) = valid_left_mins;
    left_mins(:,2) = valid_left_min_values;
    
    %% Ensure maxima are separated by at least min_separation_maxs for both feet max height regarding a window and check prominence
    
    valid_right_maxs = [];
    valid_right_max_values = [];
    
    % Loop through each maximum in right_maxs
    for i = 1:length(right_maxs)
        % Define the window around the current maximum position
        window_start = max(1, right_maxs(i, 1) - 30);
        window_end = right_maxs(i, 1);
        
        % Find the absolute maximum within this window
        [local_max_value, local_max_index] = max(right_shifted(window_start:window_end));
        % Find the absolute minimum within this window
        local_min_value = min(right_shifted(window_start:window_end));
        
        diff = local_max_value - local_min_value;
       
        if automatic_detection == 1 && usedecisions == 0  % Check AUTOMATICALLY if the difference between max and min is below threshold
            if diff < diff_max_min
                continue; % Skip this maximum if condition is met
            end
        else % use the user instructions or preloaded prompts
            if abs(diff - diff_max_min) < 10 %see if there is only a +-10 difference, so it could be a step  
                if usedecisions == 1
                    if string(decisionsArray{counter}) == "y"
                        counter = counter + 1;
                        continue; % Skip this maximum if user decides to remove it
                    elseif string(decisionsArray{counter}) == "n"
                        counter = counter + 1;
                    end
                else % do it manually
                    figure;
                    hold on;
                    plot(right_shifted, 'b', 'DisplayName', 'Right Foot');
                    plot(left_shifted, 'Color', [1 0.5 0], 'DisplayName', 'Left Foot');
            
                    plot(right_maxs(i, 1), right_maxs(i, 2), 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Inspected Max');
            
                    title(sprintf('Inspecting Maximum at frame %d - %s', right_maxs(i, 1), trialname));
                    xlim([window_start-200 window_end+200])
                    xlabel('Index');
                    ylabel('Amplitude');
                    legend('show');
                    hold off
                    
                    % Prompt user for removal confirmation
                    if diff < diff_max_min
                        prompt = 'Recommended action: Y. Do you want to remove this maximum? (y/n): ';
                    else 
                        prompt = 'Recommended action: N. Do you want to remove this maximum? (y/n): ';
                    end
                    user_input = input(prompt, 's');
                    decisions(end+1) = string(user_input);
                    
                    if strcmpi(user_input, 'y')
                        continue; % Skip this maximum if user decides to remove it
                    end
                end  
            elseif diff < diff_max_min  % Check if the difference between max and min is below threshold (and didn't fulfill last if) 
                continue; % Skip this maximum if condition is met
            end
        end
  
        % Adjust index to be absolute within the entire signal
        local_max_index = local_max_index + window_start - 1;
        
        % Check if this local maximum matches the current maximum and meets conditions
        if (local_max_index == right_maxs(i, 1)) && (isempty(valid_right_maxs) || abs(right_maxs(i, 1) - valid_right_maxs(end, 1)) >= min_separation_maxs) && ...
                (right_maxs(i, 2) >= prominence_threshold)
       
            valid_right_maxs = [valid_right_maxs; right_maxs(i, :)];
            valid_right_max_values = [valid_right_max_values; right_maxs(i, 2)];
        end
    end
    
    right_maxs = valid_right_maxs;
    
    valid_left_maxs = [];
    valid_left_max_values = [];
    
    % Loop through each maximum in left_maxs
    for i = 1:length(left_maxs)
        % Define the window around the current maximum position
        window_start = max(1, left_maxs(i, 1) - 30);
        window_end = left_maxs(i, 1);
        
        % Find the absolute maximum within this window
        [local_max_value, local_max_index] = max(left_shifted(window_start:window_end));
        % Find the absolute minimum within this window
        local_min_value = min(left_shifted(window_start:window_end));
        
        diff = local_max_value - local_min_value;
        
        if automatic_detection == 1 && usedecisions == 0
            if diff < diff_max_min
                continue; % Skip this maximum if condition is met
            end
        else
            if abs(diff - diff_max_min) < 10 %see if there is only a +-10 difference, so it could be a step. Ask the user  
                if usedecisions == 1
                    if string(decisionsArray{counter}) == "y"
                        counter = counter + 1;
                        continue; % Skip this maximum if user decides to remove it
                    elseif string(decisionsArray{counter}) == "n"
                        counter = counter + 1;
                    end
                else
                    figure;
                    hold on;
                    plot(right_shifted, 'b', 'DisplayName', 'Right Foot');
                    plot(left_shifted, 'Color', [1 0.5 0], 'DisplayName', 'Left Foot');
            
                    plot(left_maxs(i, 1), left_maxs(i, 2), 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Inspected Max');
            
                    title(sprintf('Inspecting Maximum at frame %d - %s', left_maxs(i, 1), trialname));
                    xlim([window_start-200 window_end+200])
                    xlabel('Index');
                    ylabel('Amplitude');
                    legend('show');
                    hold off;
                    % Prompt user for removal confirmation
                    if diff < diff_max_min
                        prompt = 'Recommended action: Y. Do you want to remove this maximum? (y/n): ';
                    else 
                        prompt = 'Recommended action: N. Do you want to remove this maximum? (y/n): ';
                    end
                    user_input = input(prompt, 's');
                    decisions(end+1) = string(user_input);
                    
                    if strcmpi(user_input, 'y')
                        continue; % Skip this maximum if user decides to remove it
                    end
                end
            elseif diff < diff_max_min  % Check if the difference between max and min is below threshold (and didn't fulfill last if) 
                continue; % Skip this maximum if condition is met
            end
        end

        % Adjust index to be absolute within the entire signal
        local_max_index = local_max_index + window_start - 1;
        
        % Check if this local maximum matches the current maximum and meets conditions
        if (local_max_index == left_maxs(i, 1)) && (isempty(valid_left_maxs) || abs(left_maxs(i, 1) - valid_left_maxs(end, 1)) >= min_separation_maxs) && ...
           (left_maxs(i, 2) >= prominence_threshold)
       
            valid_left_maxs = [valid_left_maxs; left_maxs(i, :)];
            valid_left_max_values = [valid_left_max_values; left_maxs(i, 2)];
        end
    end
    
    left_maxs = valid_left_maxs;
    
    
    %% Cell array with max and min
    
    % Convert right_mins and right_maxs to cell arrays
    right_mins_cell = num2cell(right_mins);
    right_maxs_cell = num2cell(right_maxs);
    
    % Combine right minima and maxima
    right_combined = [right_mins_cell; right_maxs_cell]; % Concatenate minima and maxima
    
    % Create a column for type ("min" or "max")
    type_column = [repmat({"min"}, size(right_mins_cell, 1), 1); repmat({"max"}, size(right_maxs_cell, 1), 1)];
    
    % Combine with the type column
    right_combined_with_type = [right_combined, type_column];
    
    % Sort based on the frame (first column)
    right_combined = sortrows(right_combined_with_type, 1);

    for i=1:size(right_combined(:,1))
        if strcmp(right_combined{i, 1}, NaN)
            continue
        end
        %% Eliminate certain minima 
        if strcmp(right_combined{i, 3}, "max")
            n = 1;
            max_size = size(right_combined(:,1));
    
            while n+i < max_size(1,1) %it will stop in the break always, unless we go out of the matrix. create the window in between two consecutive maxs (between i (current max) and n (next max))
                if strcmp(right_combined{n+i, 3}, "min")
                    n = n+1;
                elseif strcmp(right_combined{n+i, 3}, "max")
                    break
                end
            end
            if i + n <= max_size
                min_value = min(cell2mat(right_combined(i:i+n, 2)));
                for j=i+1:i+n-1 %+-1 to avoid the maxima
                    if right_combined{j, 2} == min_value
                        continue
                    else
                        right_combined{j, 1} = NaN;
                        right_combined{j, 2} = NaN;
                        right_combined{j, 3} = NaN;
                    end
                end
            end
        end
    end
    
    %% Clean
    % Create a logical index for rows that do not contain NaN
    rows_to_keep = ~any(cellfun(@(x) iscell(x) || (isnumeric(x) && isnan(x)), right_combined), 2);
    
    % Keep only the rows that do not contain NaN
    right_combined = right_combined(rows_to_keep, :);
    
    close all

    %% Plot of the maxima and minima
    figure;
    hold on;
    
    plot(right_shifted, 'b', 'DisplayName', 'Right Foot');
    plot(left_shifted, 'Color', [1 0.5 0], 'DisplayName', 'Left Foot');
    
    % Plot maxima and minima for right foot
    plot(right_maxs(:,1), right_maxs(:,2), 'bo', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Right Foot Maxima');
    plot(right_mins(:,1), right_mins(:,2), 'bs', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Right Foot Minima');
    
    % Plot maxima and minima for left foot
    plot(left_maxs(:,1), left_maxs(:,2), 'o', 'Color', [1 0.5 0], 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Left Foot Maxima');
    plot(left_mins(:,1), left_mins(:,2), 's', 'Color', [1 0.5 0], 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Left Foot Minima');
    
    
    xlim([0 300]);
    xlabel('Time (frames)');
    ylabel('Z-coordinate (position)');
    legend('show');
    title('Right and Left Foot Maxima and Minima');
    hold off;
    
    %% Second plot for crosses between consecutive maxima
    figure;
    hold on;
    
    % Plot the shifted right and left foot signals
    h1 = plot(right_shifted, 'b'); % Right
    h2 = plot(left_shifted, 'Color', [1 0.5 0]); % Left
    
    % Plot the maxima points for right and left feet
    h5 = plot(right_maxs(:,1), right_maxs(:,2), 'bo', 'MarkerSize', 8, 'LineWidth', 1.5); % Blue circles for right foot maxima
    h6 = plot(left_maxs(:,1), left_maxs(:,2), 'o', 'Color', [1 0.5 0], 'MarkerSize', 8, 'LineWidth', 1.5); % Orange circles for left foot maxima
    
    % Plot crosses at the midpoint between consecutive maxima for right foot
    for i = 1:length(right_maxs)-1
        x_mid_right = (right_maxs(i,1) + right_maxs(i+1,1)) / 2; % Midpoint on x-axis
        x_mid_right_rounded(i,1) = round(x_mid_right); % Round to nearest integer
        y_mid_right = right_shifted(x_mid_right_rounded(i,1)); % Corresponding y-value
        h3 = plot(x_mid_right_rounded(i,1), y_mid_right, 'bx', 'MarkerSize', 10, 'LineWidth', 2); % Blue cross for right step
    end
    
    % Plot crosses at the midpoint between consecutive maxima for left foot
    for i = 1:length(left_maxs)-1
        x_mid_left = (left_maxs(i,1) + left_maxs(i+1,1)) / 2; % Midpoint on x-axis
        x_mid_left_rounded(i,1) = round(x_mid_left); % Round to nearest integer
        y_mid_left = left_shifted(x_mid_left_rounded(i,1)); % Corresponding y-value
        h4 = plot(x_mid_left_rounded(i,1), y_mid_left, 'x', 'Color', [1 0.5 0], 'MarkerSize', 10, 'LineWidth', 2); % Orange cross for left step
    end
    
    % Set x-axis limits
    xlim([0 300]);
    xlabel('Time (frames)');
    ylabel('Z-coordinate (position)');
    legend([h1 h2 h3 h4 h5 h6], {'Right Foot (RHEEL)', 'Left Foot (LHEEL)', 'Right Step', 'Left Step', 'Right Foot Maxima', 'Left Foot Maxima'});
    title('Right and Left Foot Steps');
    hold off;
    
    %% Find Step Times and store them with the user decisions

    fs = trial.(trialname).FrameRate;
    steps = struct();

    steps.right.frames = x_mid_right_rounded; % step frame
    steps.right.secs = steps.right.frames ./ fs; % step times indexed by time in sec

    steps.left.frames = x_mid_left_rounded; % step times indexed by frame
    steps.left.secs = steps.left.frames ./ fs; % step times indexed by time in sec

    steps.counter = counter;

    if ~isfield(steps, 'decisions')
        steps.decisions = [];  % Initialize as an empty array
    end
    if ~isempty(decisions)
        steps.decisions = [steps.decisions,decisions];
    end
end