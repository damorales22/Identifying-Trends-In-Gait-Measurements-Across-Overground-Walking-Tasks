%% Runs Heel Strikes Function
clear
clc
close all

% Variables
filepath = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\BMH10(September)\MC\";
EMGpath = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\BMH10(September)\EMG\";
trialname = "Trial0004";

% Bandpass filter inputs
f_low = 30;
f_high = 500;
n = 1; % fourth order // bandpass filter and filtfilt both double the order

% lowpass filter inputs
f_lowpass = 6;
m = 2; % fourth order // filtfilt doubles the order

% Load Mocap data
filename = filepath + trialname;
trial = load([num2str(filename) '.mat']);

% run heel strikes
threshold = CalculateThreshold(trial, trialname);
strikes = LocateHeelStrikes(trial, trialname, threshold);

% Load EMG data
EMGname = EMGpath + trialname;
EMGdata = load([num2str(EMGname) '.mat']);
fs = EMGdata.Fs(1); %sampling frequency


for i = 1:12 %12 channels
    data = EMGdata.Data(i, :); %data from that channel
    
    % Bandpass Filter
    [b, a] = butter (n, [f_low f_high]/(fs/2));
    EMG_bp = filtfilt(b, a, data); 
    
    % Rectify
    EMG_abs = abs(EMG_bp);
    
    % Lowpass Filter
    [b, a] = butter(m, f_lowpass/(fs/2));
    EMG_filtered = filtfilt(b, a, EMG_abs);
    clf;

    %figure()
    hold on
    plot(EMG_filtered)
    if i <=6 %Plot the points; right
        scatter(strikes.right.secs.*fs, EMG_filtered(round(strikes.right.secs.*fs))) %scatter to represent as circles
    else %left
        scatter(strikes.left.secs.*fs, EMG_filtered(round(strikes.left.secs.*fs)))
    end
    hold off
    title("Channel " + i)
    axis("tight")
    
    % figure
    % plot(data)
    % axis("tight")
    % 
    % figure
    % plot(EMG_bp)
    % axis("tight")
    % 
    % figure
    % plot(EMG_abs)
    % axis("tight")

end

