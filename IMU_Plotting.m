%% Clear Previous
clear 
clc
close all

%% Create Variables 
IMUcodes = ["742", "722", "511", "589", "668", "598", "568"];
IMUlocations = ["Right Foot", "Left Foot", "Right Thigh", "Left Thigh", "Pelvis", "Left Shank", "Right Shank"];
numtrials = 2;
trials = ["jumping jacks", "walking"]; % trial types
trial = 1; % choose trial number for Figure 3

crop_start = 700;
crop_end = 1000;

%% Plot Accelerations
figure(1)
for i = 1:length(IMUcodes) % iterate through IMU sensors
    subplot(4, 2, i) % 2x4, one plot for each sensor
    hold on
    sensor = IMUcodes(i); 
    for j = 1:numtrials % iterate through trials

        % load data
        filename = ['june_test_' num2str(j) '_00340' num2str(sensor) '.txt'];
        IMUdata = readtable(['C:\Ability Lab\IMU Processing\imudata_june\' filename]);
        % Column 1: Timestamp; Columns 2-4: Accelerometers X, Y, Z
        
        % demean data
        means = mean(IMUdata); % calculate mean values

        for k = 2:4 % iterate through dimensions (x, y, z) and subtract mean
           IMUdata(:, k) = IMUdata(:, k) - means(1, k);
        end

        % calculate vector magnitudes and convert to array
        IMUdata_2 = IMUdata(:, 2:4).^2;
        accel = table2array(sqrt(sum(IMUdata_2, 2)));
        
        % Create Time Vector
        time = 1:length(accel);

        % plot accelerometer data
        plot(time, accel)

    end

    % format plot
    xlim tight
    % ylim([-20, 50])
    title(num2str(IMUlocations(i)))
    xlabel ('Time')
    ylabel('Acceleration (m/s^2)')
    legend([num2str(trials(1))], [num2str(trials(2))])
    hold off
end

%% Plot Zoomed Trials

figure(2)
for i = 1:length(IMUcodes) % iterate through IMU sensors
    subplot(4, 2, i) % 2x4, one plot for each sensor
    hold on
    sensor = IMUcodes(i); 
    for j = 1:numtrials % iterate through trials
        
         % load data
        filename = ['june_test_' num2str(j) '_00340' num2str(sensor) '.txt'];
        IMUdata = readtable(['C:\Ability Lab\IMU Processing\imudata_june\' filename]);
        % Column 1: Timestamp; Columns 2-4: Accelerometers X, Y, Z
        
        % demean data
        means = mean(IMUdata); % calculate mean values

        for k = 2:4 % iterate through dimensions (x, y, z) and subtract mean
           IMUdata(:, k) = IMUdata(:, k) - means(1, k);
        end

        % calculate vector magnitudes and convert to array
        IMUdata_2 = IMUdata(:, 2:4).^2;
        accel = table2array(sqrt(sum(IMUdata_2, 2)));
        
        % Create Time Vector
        time = 1:length(accel);

        % plot accelerometer data
        plot(time, accel)
    end

    % format plot
    xlim([crop_start, crop_end])
    % ylim([-20, 50])
    title(num2str(IMUlocations(i)))
    xlabel ('Time')
    ylabel('Acceleration (m/s^2)')
    legend([num2str(trials(1))], [num2str(trials(2))])
    hold off
end


%% Plot Single Zoomed Trial

figure(3)
for i = 1:length(IMUcodes) % iterate through IMU sensors
    subplot(4, 2, i) % 2x4, one plot for each sensor
    hold on
    sensor = IMUcodes(i); 
    for j = trial
        
         % load data
        filename = ['june_test_' num2str(j) '_00340' num2str(sensor) '.txt'];
        IMUdata = readtable(['C:\Ability Lab\IMU Processing\imudata_june\' filename]);
        % Column 1: Timestamp; Columns 2-4: Accelerometers X, Y, Z
        
        % demean data
        means = mean(IMUdata); % calculate mean values

        for k = 2:4 % iterate through dimensions (x, y, z) and subtract mean
           IMUdata(:, k) = IMUdata(:, k) - means(1, k);
        end

        % calculate vector magnitudes and convert to array
        IMUdata_2 = IMUdata(:, 2:4).^2;
        accel = table2array(sqrt(sum(IMUdata_2, 2)));
        
        % Create Time Vector
        time = 1:length(accel);

        % plot accelerometer data
        plot(time, accel)
    end

    % format plot
    xlim([crop_start, crop_end])
    % ylim([-20, 50])
    title(num2str(IMUlocations(i)))
    xlabel ('Time')
    ylabel('Acceleration (m/s^2)')
    legend([num2str(trials(j))])
    hold off
end

%% Plot Gyroscope 
figure(4)
for i = 1:length(IMUcodes) % iterate through IMU sensors
    subplot(4, 2, i) % 2x4, one plot for each sensor
    hold on
    sensor = IMUcodes(i); 
    for j = 1:numtrials % iterate through trials

        % load data
        filename = ['june_test_' num2str(j) '_00340' num2str(sensor) '.txt'];
        IMUdata = readtable(['C:\Ability Lab\IMU Processing\imudata_june\' filename]);
        % Column 1: Timestamp; Columns 8-10: Gyroscope X, Y, Z
        
        % demean data
        means = mean(IMUdata); % calculate mean values

        for k = 8:10 % iterate through dimensions (x, y, z) and subtract mean
           IMUdata(:, k) = IMUdata(:, k) - means(1, k);
        end

        % calculate vector magnitudes and convert to array
        IMUdata_2 = IMUdata(:, 8:10).^2;
        gyro = table2array(sqrt(sum(IMUdata_2, 2)));
        
        % Create Time Vector
        time = 1:length(gyro);

        % plot accelerometer data
        plot(time, gyro)

    end

    % format plot
    xlim tight
    % ylim([-20, 50])
    title(num2str(IMUlocations(i)))
    xlabel ('Time')
    ylabel('Angular Rate (rad/s)')
    legend([num2str(trials(1))], [num2str(trials(2))])
    hold off
end

%% Plot Zoomed Gyroscope

figure(5)
for i = 1:length(IMUcodes) % iterate through IMU sensors
    subplot(4, 2, i) % 2x4, one plot for each sensor
    hold on
    sensor = IMUcodes(i); 
    for j = 1:numtrials % iterate through trials
        
         % load data
        filename = ['june_test_' num2str(j) '_00340' num2str(sensor) '.txt'];
        IMUdata = readtable(['C:\Ability Lab\IMU Processing\imudata_june\' filename]);
        % Column 1: Timestamp; Columns 2-4: Accelerometers X, Y, Z
        
        % demean data
        means = mean(IMUdata); % calculate mean values

        for k = 8:10 % iterate through dimensions (x, y, z) and subtract mean
           IMUdata(:, k) = IMUdata(:, k) - means(1, k);
        end

        % calculate vector magnitudes and convert to array
        IMUdata_2 = IMUdata(:, 8:10).^2;
        gyro = table2array(sqrt(sum(IMUdata_2, 2)));
        
        % Create Time Vector
        time = 1:length(gyro);

        % plot accelerometer data
        plot(time, gyro)
    end

    % format plot
    xlim([crop_start, crop_end])
    % ylim([-20, 50])
    title(num2str(IMUlocations(i)))
    xlabel ('Time')
    ylabel('Angular Rate (rad/s)')
    legend([num2str(trials(1))], [num2str(trials(2))])
    hold off
end