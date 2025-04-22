%% Calculate Threshold Function
% Used alongside the heel strike function, this function takes a trial and calculates the ideal threshold based on the trial's walking speed via the
% WalkingSpeed function. Important for the EMG data, as it must be very precise. The threshold is the max magnitude of the velocity of the foot in all 
% x, y and z directions put together, and anything above the threshold is considered the foot moving. Below is considered the foot on the ground. 
% Trial and error was done to determine which threshold is best for speeds slow medium and fast, like inputting values, running the function, and seeing 
% which value worked best to get the exact point of heel strike, and then once she had 6 values she plotted it and just found the best fit line of
% those 6 points. The threshold is the max velocity that we still consider as “still”

% function[threshold] = CalculateThreshold(trial, trialname)
%     % get markerdata
%     markers = trial.(trialname).Trajectories.Labeled.Labels;
%     markerdata = struct();
%     for i = 1:length(markers)
%        markerdata.(char(markers(i))) = trial.(num2str(trialname)).Trajectories.Labeled.Data(i, :, :);
%     end
%     markerdata.FrameRate = trial.(trialname).FrameRate;
%     markerdata.Frames = trial.(trialname).Frames;
% 
%     % get average speed
%     speed = WalkingSpeed(trial, markerdata, trialname);
% 
%     % calculate threshold
%     threshold = 10.95*speed.avg - 6.56;
% end

function[threshold] = CalculateThreshold(trial, trialname)
    % get markerdata
    markers = trial.(trialname).Trajectories.Labeled.Labels;
    markerdata = struct();
    for i = 1:length(markers)
       markerdata.(char(markers(i))) = trial.(num2str(trialname)).Trajectories.Labeled.Data(i, :, :);
    end

    %% ASSUMPTION: the moment and position when each foot (MTH5 or LMH5) is (approximately) the moment when the z coordinate of each marker cross each other (observed in the MoCap representation) 
    % Since the don't really care about the timeframe but about the position, this method is okay.

    right_position = squeeze(markerdata.RMTH5(1,3,:)); % ONLY z coordinate
    left_position = squeeze(markerdata.LMTH5(1,3,:)); % ONLY z coordinate

    diff_signal = left_position - right_position;
    crossing_indices = find(diff_signal(1:length(diff_signal)-1) .* diff_signal(2:length(diff_signal)) < 0);

    % Pre-allocate arrays for storing results
    right_step = [];
    left_step = [];
    crossing_points = [];
    
    for i = 1:length(crossing_indices)
        idx = crossing_indices(i);
        
        % Check the slope at the crossing point for each signal
        left_slope = left_position(idx+1) - left_position(idx);
        right_slope = right_position(idx+1) - right_position(idx);
        
        % Store the crossing point and the foot that is decreasing
        crossing_points = [crossing_points; idx]; % or replace idx with the actual x-value if needed
        if left_slope < 0
            left_step = [left_step; idx];
        elseif right_slope < 0
            right_step = [right_step; idx];
        end
    end

    figure
    hold on
    plot(500:1000, right_position(500:1000))
    plot(500:1000, left_position(500:1000))
    legend('Right foot', 'Left foot')
    hold off

    % markerdata.FrameRate = trial.(trialname).FrameRate;
    % markerdata.Frames = trial.(trialname).Frames;
    % 
    % % get average speed
    speed = WalkingSpeed(trial, markerdata, trialname);
    % 
    % % calculate threshold
    threshold = 10.95*speed.avg - 6.56;
end