%% This Function takes a loaded trial and marker data as an input and returns the frames that correspond to heel strikes for each foot. LocateHeelStrikes gives us 
% a much more accurate timestamp but the end of the step isn’t considered so it’s definetly still useful- for segmentation each gait cycle the heel strike has to be 
% super accurate for emg processing

function[strikes] = LocateHeelStrikes(trial, trialname, threshold)
    %% Variables
    delaytime = 1.5; % time in seconds before first located heel strike
    checktime = 0.2;
    upperthreshold = 12; % threshold that must be exceeded before next window can be counted (mm/frame)

    % lowpass filter inputs
    f_low = 6;
    m = 2; % fourth order // filtfilt doubles the order
    
    %% Separates out markers
    markers = trial.(trialname).Trajectories.Labeled.Labels;
    markerdata = struct();
    for i = 1:length(markers)
       markerdata.(char(markers(i))) = trial.(num2str(trialname)).Trajectories.Labeled.Data(i, :, :);
    end
    
    %% Locate Heel Strike Times
    % look for times where foot is stationary
    
    % extract all coordinates from heels
    coords_r = squeeze(markerdata.RHEEL(1,1:3,:))'; % x, y, z in columns
    coords_l = squeeze(markerdata.LHEEL(1,1:3,:))'; 
    
    % calculate velocity
    v_r = (coords_r(2:length(coords_r), :) - coords_r(1:length(coords_r)-1, :));
    v_mag_r = sqrt((v_r(:,1)).^2 + (v_r(:,2)).^2 + (v_r(:,3)).^2); %Magnitude
    
    v_l = (coords_l(2:length(coords_l), :) - coords_l(1:length(coords_l)-1, :));
    v_mag_l = sqrt((v_l(:,1)).^2 + (v_l(:,2)).^2 + (v_l(:,3)).^2);
    
    % filter data
    framerate = trial.(trialname).FrameRate; % get frame rate
    [b, a] = butter(m, f_low/(framerate/2), "low"); %Low-pass Butterworth filter
    try
        v_mag_r_filt = filtfilt(b, a, v_mag_r); %Data is filtered twice—once. This ensures that there is no phase shift. Important for preserving the timing of the velocity signals.
        v_mag_l_filt = filtfilt(b, a, v_mag_l);
    catch
        error("Error filtering. There are probably missing frames. Check v_mag_r and v_mag_l, and qtm files")
    end
    
    % create time window
    time_window = 0.05*framerate; % check 0.1 second intervals
    
    % find window starts and ends for foot
    windows_r = NaN(length(v_mag_r_filt),2); % first column for start times, second for end times
    i = 1; % indexing through frames
    j = 1; % count identified windows
    
    while i < length(v_mag_r_filt) - time_window
        if ((v_mag_r_filt(i) < threshold) && (v_mag_r_filt(i+time_window) < threshold)) %Threshold is not longer defined. However, you could adjust a line to the crossing points
            windows_r(j,1) = i; % set start time
            while v_mag_r_filt(i+time_window) < threshold && (i < length(v_mag_r_filt) - time_window)
                i = i+1; % index through frames until threshold exceeded
            end
            i = i + time_window - 1; % go to end of window
            windows_r(j,2) = i; % set end time
            j = j+1; % next window
        else 
            i = i+1; % next frame
        end
    end
    
    % find window starts and ends for left foot
    windows_l = NaN(length(v_mag_l_filt),2); % first column for start times, second for end times
    i = 1; % indexing through frames
    j = 1; % count identified windows
    
    while i < length(v_mag_l_filt) - time_window
        if ((v_mag_l_filt(i) < threshold) && (v_mag_l_filt(i+time_window) < threshold))
            windows_l(j,1) = i; % set start time
            while (v_mag_l_filt(i+time_window) < threshold) && (i < length(v_mag_l_filt) - time_window)
                i = i+1; % index through frames until threshold exceeded
            end
            i = i + time_window - 1; % go to end of window
            windows_l(j,2) = i; % set end time
            j = j+1; % next window
        else 
            i = i+1; % next frame
        end
    end
    
    % check that windows aren't back to back
    window_count = 2;
    for i=2:length(windows_r(:,1)) - 1
        if window_count >= length(windows_r(:,1))
            break
        elseif windows_r(window_count, 1) - windows_r(window_count - 1, 2) < checktime * framerate
            windows_r(window_count,1) = NaN;
            windows_r(window_count,2) = NaN;
            windows_r = windows_r(~isnan(windows_r(:, 1)), :); % remove row
        else
            window_count = window_count + 1;
        end
    end
    
    window_count = 2;
    for i=2:length(windows_l(:,1)) - 1
        if window_count >= length(windows_l(:,1))
            break
        elseif windows_l(window_count,1) - windows_l(window_count-1,2) < checktime * framerate
            windows_l(window_count,1) = NaN;
            windows_l(window_count,2) = NaN;
            windows_l = windows_l(~isnan(windows_l(:, 1)), :); % remove row
        else
            window_count = window_count +1;
        end
    end
    
    % remove windows before delay
    windows_r = windows_r(windows_r(:, 1) >= delaytime * framerate, :); 
    windows_l = windows_l(windows_l(:, 1) >= delaytime * framerate, :);

    % remove windows where threshold is not exceeded since previous
    window_count = 1;
    for i = 1:length(windows_r) - 1
        % get segment
        segment = v_mag_r(windows_r(window_count, 1):windows_r(window_count + 1, 1));
        if any(segment > upperthreshold)
            window_count = window_count + 1;
        else
            windows_r(window_count + 1, 1) = NaN;
            windows_r(window_count + 1, 2) = NaN;
            windows_r = windows_r(~isnan(windows_r(:, 1)), :); % remove row
        end
    end

    window_count = 1;
    for i = 1:length(windows_l) - 1
        % get segment
        segment = v_mag_l(windows_l(window_count, 1):windows_l(window_count + 1, 1));
        if any(segment > upperthreshold)
            window_count = window_count + 1;
        else
            windows_l(window_count + 1, 1) = NaN;
            windows_l(window_count + 1, 2) = NaN;
            windows_l = windows_l(~isnan(windows_l(:, 1)), :); % remove row
        end
    end

    
    % Find Heel Strike Times
    strikes = struct();
    %strikes.right.frames = round(windows_r(:, 1)); % strike times indexed by frame; only starting point
    strikes.right.frames = round(windows_r); %starting and ending point
    strikes.right.secs = strikes.right.frames ./ framerate; % strike times indexed by time in sec
    
    %strikes.left.frames = round(windows_l(:, 1)); % step times indexed by frame; only starting point
    strikes.left.frames = round(windows_l); %starting and ending point
    strikes.left.secs = strikes.left.frames ./ framerate; % step times indexed by time in sec
    
    %% Check with Plot
    % numframes = trial.(trialname).Frames;
    % trialtime = numframes/framerate;
    % steps = 1/framerate;
    % t = 3*steps/2:steps:trialtime - steps/2;
    % 
    % % get start and end point values
    % s_points_r = v_mag_r(windows_r(:, 1));
    % e_points_r = v_mag_r(windows_r(:, 2));
    % 
    % s_points_l = v_mag_l(windows_l(:, 1));
    % e_points_l = v_mag_l(windows_l(:, 2));
    % 
    % % plot 
    % figure();
    % title("Find Steps Graph")
    % hold on
    % 
    % % plot(t, v_mag_l_filt)
    % % plot(t, v_mag_r_filt)
    % plot(t,v_mag_r) %Visualize the points where the heel strike (or stationary period) begins and ends for both feet
    % plot(t,v_mag_l)
    % 
    % scatter(windows_r(:, 1)/framerate, s_points_r)
    % scatter(windows_r(:, 2)/framerate, e_points_r)
    % 
    % scatter(windows_l(:, 1)/framerate, s_points_l)
    % scatter(windows_l(:, 2)/framerate, e_points_l)
    % 
    % xlim([0 10])
    % legend('right', 'left')
    % hold off
end