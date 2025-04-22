%% Walking Speed Function
% Calculates the average wlaking speed during a trial, taking the loaded
% trial and markers as inputs and outputting the average walking speed

function [speed, stride_time] = WalkingSpeed(trial, markerdata, trialname, steps)
    % isolate location of pelvis
    %coords_RASIS = squeeze(markerdata.RASIS(1,1:2,:))'; % x, y in columns
    %coords_LASIS = squeeze(markerdata.LASIS(1,1:2,:))';
    try
        coords_RASIS = permute(markerdata.RASIS(1,1:2,:), [3, 2, 1]); %rearrange dimensions (eliminate the first dimension) and only take x,y
    catch
        coords_RASIS = permute(markerdata.RPSIS(1,1:2,:), [3, 2, 1]); %rearrange dimensions (eliminate the first dimension) and only take x,y
    end
    try
        coords_LASIS = permute(markerdata.LASIS(1,1:2,:), [3, 2, 1]);
    catch
        coords_LASIS = permute(markerdata.LPSIS(1,1:2,:), [3, 2, 1]);
    end

    coords(:, 1) = (coords_RASIS(:, 1) + coords_LASIS(:, 1))./2; %calculate mid point for x
    coords(:, 2) = (coords_RASIS(:, 2) + coords_LASIS(:, 2))./2; %calculate mid point for y
    coords(any(isnan(coords), 2), :) = [];

    if isempty(coords)
        coords_RASIS = permute(markerdata.RPSIS(1,1:2,:), [3, 2, 1]);
        coords_LASIS = permute(markerdata.LPSIS(1,1:2,:), [3, 2, 1]);
        clear coords
        coords(:, 1) = (coords_RASIS(:, 1) + coords_LASIS(:, 1))./2; %calculate mid point for x
        coords(:, 2) = (coords_RASIS(:, 2) + coords_LASIS(:, 2))./2; %calculate mid point for y
        coords(any(isnan(coords), 2), :) = [];
    end

    % Trim the first 10 and last 10 seconds
    framerate = trial.(trialname).FrameRate;
    numframes = trial.(trialname).Frames;
    trim_time = 10;
    try
        coords_trimmed = coords(trim_time*framerate:numframes - trim_time*framerate, :);
    catch
        coords_trimmed = coords;
    end

    % Find distance travelled by the pelvis
    dist = 0;
    dist_byframe_trimmed = 0;
    for i = 1:length(coords_trimmed)-1
        disttemp = sqrt((coords_trimmed(i, 1)-coords_trimmed(i+1, 1))^2 + (coords_trimmed(i, 2)-coords_trimmed(i+1, 2))^2); %Pythagoras th (sqrt(x^2+y^2))
        if ~isnan(disttemp)
            dist = dist + disttemp;
            dist_byframe_trimmed(i) = dist;
        end
    end
    
    %% Calculate walking speed and var (std) - Method 1
    incrementalDistance = [0, diff(dist_byframe_trimmed)]; % Distance per frame (first frame has no previous frame)
    timePerFrame = 1 / framerate; % Time duration of each frame (seconds)
    speed_overtime = (incrementalDistance / timePerFrame) / 1000; % Speed in meters per second (m/s)

    trialtime = numframes/framerate - 20; % account for trimming (10+10 secs)
    avg = dist/trialtime;
    speed.avg = avg/1000; % convert from mm/s to m/s
    %speed.std = std(speed_overtime(2:end));
    
    %% Calculate walking speed and var (std) - Method 2: for each step (additional calculation of the stride time)
    % Identify first feet moved to calculate stride times
    if steps.right.frames(1) > steps.left.frames(1)
        first_feet = "left";
        second_feet = "right";
    else
        first_feet = "right";
        second_feet = "left";
    end

    steps_ordered = [];
    i = 1; % Index for first feet
    j = 1; % Index for second feet
    
    while i <= length(steps.(first_feet).frames) && j <= length(steps.(second_feet).frames)
        
        max_i = i; % Store the initial index of first feet
        while i <= length(steps.(first_feet).frames) && steps.(first_feet).frames(i) < steps.(second_feet).frames(j)
            max_i = i;
            i = i + 1;
        end
        
        if steps.(first_feet).frames(max_i) < steps.(second_feet).frames(j)
            steps_ordered = [steps_ordered; steps.(first_feet).frames(max_i), steps.(second_feet).frames(j)]; %First column: when the step starts. Second column: when the step ends
            j = j + 1;
        else
            j = j + 1;
        end
    end

    % Find distance travelled by the pelvis
    dist = 0;
    dist_byframe = 0;
    for i = 1:length(coords)-1
        disttemp = sqrt((coords(i, 1)-coords(i+1, 1))^2 + (coords(i, 2)-coords(i+1, 2))^2); %Pythagoras th (sqrt(x^2+y^2))
        dist = dist + disttemp;
        dist_byframe(i) = dist;
    end

    % Walking speed std: calculation of the speed for each step, and then do the std of all
    steps_concatenated = reshape(steps_ordered.', 1, []);
    try
        for i=1:size(steps_concatenated,2)-1
            distance_step = (dist_byframe(steps_concatenated(i+1))-dist_byframe(steps_concatenated(i)))/1000; % Distance covered in that step (in meters)
            time_step = (steps_concatenated(i+1) - steps_concatenated(i))/framerate;
            speed_step(i) = distance_step/time_step;
        end
        speed_step = speed_step(~isnan(speed_step));
    
        speed_step = speed_step(6:end-5); % Remove the first and last 5 steps
        %speed.avg = mean(speed_step); % I think it is better to calculate the speed in the first way, but the result is the almost the same (+- 0.01 changes)
        speed.std = std(speed_step);
    catch
        speed.std = NaN;
    end

    % Stride times calulation
    stride_times = [diff(steps_ordered(:, 1)); diff(steps_ordered(:, 2))]; % frames
    stride_time_avg = mean(stride_times)/framerate;
    stride_time_std = std(stride_times)/framerate;
    stride_time_median = median(stride_times)/framerate;

    stride_time.avg = stride_time_avg;
    stride_time.std = stride_time_std;
    stride_time.median = stride_time_median;

end
