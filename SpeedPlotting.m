%% Plotting Speed over time for trials
clear 
clc
close all

%% create variables
numtrials = 9; 
delaytime = 1.5; % ignore troughs before this time

% lowpass filter inputs
f_lowpass = 6;
m = 2; % fourth order // filtfilt doubles the order

% Define a set of colors
colors = colormap(turbo(numtrials));

%% Get Data
figure()
hold on
for i = 1:numtrials
    trialname = sprintf('Trial%04d', i);
    
    % Loads Trial Data
    trial = load(['C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\Trial1\MC\' num2str(trialname) '.mat']); % load trial
    framerate = trial.(trialname).FrameRate;
    numframes = trial.(trialname).Frames;

    % Separates out markers
    markers = trial.(trialname).Trajectories.Labeled.Labels;
    markerdata = struct();
    for j = 1:length(markers)
       markerdata.(char(markers(j))) = trial.(num2str(trialname)).Trajectories.Labeled.Data(j, :, :);
    end

    % Track Pelvis
    coords_RASIS = squeeze(markerdata.RASIS(1,1:2,:))'; % x, y in columns
    coords_LASIS = squeeze(markerdata.LASIS(1,1:2,:))';

    coords(:, 1) = (coords_RASIS(:, 1) + coords_LASIS(:, 1))./2;
    coords(:, 2) = (coords_RASIS(:, 2) + coords_LASIS(:, 2))./2;

   
    % calculate velocity
    v = (coords(1:length(coords)-1, :) - coords(2:length(coords), :));
    v_mag = sqrt((v(:,1)).^2 + (v(:,2)).^2).*framerate;
    v_mag = v_mag./1000;

    % filter data
    [b, a] = butter(m, f_lowpass/(framerate/2), "low");
    v_mag_filt = filtfilt(b, a, v_mag);
    
    % smooth data
    v_mag_smooth = smooth(v_mag, 5000);
    
    % plot
    labels = cell(1, numtrials);
    for j = 1:numtrials % create labels
        labels{j} = ['Trial ' num2str(j)];
    end

    step = 1/framerate;
    trialtime = numframes/framerate;
    t = 3*step/2:step:trialtime - step/2;
    plot(t, v_mag_smooth, "Color", colors(i, :))

end

xlim([10 trialtime-10])
legend(labels, "Location", "southeast")
ylabel("Speed (m/s)")
xlabel("Time (s)")
title("Speed vs. Time", 'FontSize', 20)
hold off