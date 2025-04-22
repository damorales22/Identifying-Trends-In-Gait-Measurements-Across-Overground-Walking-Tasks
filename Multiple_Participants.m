%% Take the MoCap data.mat dataset of each participant and analyze
clear 
clc
close all

%% User parameters
path = 'C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\'; 

n = 3; %Column 3 LONGITUDINAL ERROR (distance to the infinite line), Column 4 SIGNED LONGITUDINAL ERROR, Column 5 TOTAL ERROR (distance to the red tape), Column 6 LATERAL ERROR; (of matches.side.match)
section = 2; %1 for ALL, 2 for only STRAIGHT, 3 for only CURVES
foot = 1; %1 for BOTH FEET analysis (combined), 2 for only RIGHT foot, 3 for only LEFT foot

plot1 = 0; %Errors grouped by condition and overall bar plot   
    subtract_all_1 = 0;    %0 do nothing, 1 to subtract the error mean across all trials to the error bars. ONLY applied when using the ABOSLUTE or SIGNED error and ONLY when grouped BY SPEED.
    subtract_speeds_1 = 0; %0 do nothing, 1 to subtract the error mean according to each speed to the error bars. ONLY applied when using the ABOSLUTE or SIGNED error and ONLY when grouped BY SPEED.

plot2 = 0; %Bar Plot for Accuracy Conditions (3; any) and a given Balance and Speed
    speed_input_2 = 2;   % 0=slow, 1=medium, 2=fast, 3=any
    balance_input_2 = 0; % 0=none, 1=medium, 2=high, 3=any; (disturbances)

    % Substracting the means do not make sense in this plot, as the means are JUST for a=2
    subtract_all_2 = 0;    %0 do nothing, 1 to subtract the error mean across all trials to the error bars. ONLY applied when using the ABOSLUTE or SIGNED error and ONLY when grouped by speed.
    subtract_speeds_2 = 0; %0 do nothing, 1 to subtract the error mean according to each speed to the error bars. ONLY applied when using the ABOSLUTE or SIGNED error and ONLY when grouped by speed.

plot4 = 0; %Customized multiple plots of the given conditions
    speed_input = 3;    % 0=slow, 1=medium, 2=fast, 3=any
    accuracy_input = 0; % 0=none, 1=medium, 2=high, 3=any
    balance_input = 0;  % 0=none, 1=medium, 2=high, 3=any; (disturbances)

    subtract_all_4 = 0;    %0 do nothing, 1 to subtract the error mean ACROSS ALL TRIALS to the error bars. ONLY applied when using the ABOSLUTE or SIGNED error
    subtract_speeds_4 = 0; %0 do nothing, 1 to subtract the error mean according to EACH SPEED to the error bars. ONLY applied when using the ABOSLUTE or SIGNED error

plot5 = 0; % Plot the metabolic data across all participants
    normalize = 1; %Normalize to the weigth of each participant

plot6 = 0; % Plot all conditions
   plot_number = 4; % 1: speed vs metabolics. 2: speed vs step width. 3: speed vs step length. 4: speed vs step error.
         % 5: metabolics vs step width. 6: metabolics vs step length.  7: metabolics vs step error. 8: step width vs step length. 
         % 9: step width vs step error. 10: step length vs step error. 11: speed vs width variability. 12: metabolics vs width variability. 
         % 13: step error vs width variability. 14: step error vs head angle. 15: step error vs head angle (just colors)
         % 16: step error vs head angle (shape and colors). 17: step error vs trunk angle. 18: step error vs trunk angle (shape and colors). 
         % 19: speed vs stride time. 20: trunk angle vs step width. 21: trunk angle vs step width var. 22: stride time vs step length var. 
         % 23: stride time var vs step length var. 24: stride time var vs step error 25: step width var vs step length. 
         % 26: step width var vs step length var. 27: step width var vs step width 28: step width var vs metabolics.
         % 29: step length vs metabolics. 30: step length var vs metabolics. 31: step width vs metabolics. 32: slv vs speed
         % 33: EE vs walking speed var. 34: Trunk angle vs walking speed var. 35: Trunk angle vs EE. 
         % 36: SLV vs step width. 37: Stride time vs step width. 38: Stride time var vs step width. 39: Head angle vs slv. 
         % 40: Trunk angle vs slv. 41: Stride time vs stride time var. 42: Trunk angle vs head angle. 43: Head angle vs step error var.
         % 44: step error vs step error var. 45: head angle vs step width
                   
   color_shape_option_1_1 = 0; % 0 - all for THAT plot number. 1 - Color: accuracy. Shape: balance. 2 - Color: balance. Shape: speed. 3 - Color: speed. Shape: accuracy              
   color_shape_option_1_2 = 0; % 0 - off. 1 - Run ALL plots for that color_shape_option_1_1
   show_cluster_analysis1 = 0; % Show clustering plots
   normalize_metabolics1 = 1;
   show_line = 1; % Plot regression line
       show_formula1 = 1;
   subtract_all_6 = 0;    % Only when using the step error
   subtract_speeds_6 = 0; % Only when using the step error

plot7 = 0; % Participants' performance and propmts 
   shaded_plot = 0;
   plot_eachparticipant = 1; % Points in the bar plots representing each participant
   subtract_all_7 = 0;    % Only when using the step error. Not used anymore
   subtract_speeds_7 = 0; % Only when using the step error. A baseline for each speed (3)

plot8 = 0; % Compare the baselines; the only output is the txt

plot9 = 1; % Survey responses
   color_shape_option_2 = 1; % 1 - Color: accuracy. Shape: balance. 2 - Color: balance. Shape: speed. 3 - Color: speed. Shape: accuracy
   show_cluster_analysis2 = 0;
   show_formula2 = 1;
   normalize_metabolics2 = 1;

plot10 = 0; % Pairwise correlation matrix
   option_correlations = 1; % Choose what correlations to show. 1: all. 2: close to +1. 3: close to -1. 4: close to 0
   show_all_plots = 0; % 0 to only show the main one, 1 to show also the 9 plots (one for each condition). The excel file is generated anyway

%% Load data and define constants

subfolders = dir(path);
data_names = {};
imported_data = struct();

% The structure of the data should be path\BMHxx\MC\ and path\BMHxx\Metabolics.csv 
for i = 1:length(subfolders)
    folder_name = fullfile(subfolders(i).name, 'MC');  % Construct folder name correctly

    if startsWith(folder_name, 'BMH') 

        metabolics_file_path = fullfile(path, subfolders(i).name, "Metabolics.csv");
        if exist(metabolics_file_path, 'file')
            loaded_metabolics = readtable(metabolics_file_path);  % Use readtable if it's a CSV file
            bmh = subfolders(i).name(1:5);  % Get the correct BMH name
            metabolics.(bmh) = loaded_metabolics;  % Store metabolic data
        end

        if exist(fullfile(path, folder_name, 'data.mat'), 'file') || exist(fullfile(path, folder_name, 'data_' + string(bmh) + '.mat'), 'file')
            try
                data_path = fullfile(path, folder_name, 'data.mat');
                loaded_data = load(data_path);
            catch
                data_path = fullfile(path, folder_name, 'data_' + string(bmh) + '.mat');
                loaded_data = load(data_path);
            end
            data_names{end+1} = bmh; 
            imported_data.(bmh) = loaded_data.data; 
        end
    end
end

colors = [
    0.1 0.6 0.1; % Green for Slow
    0.1 0.1 0.8; % Blue for Medium
    0.8 0.1 0.1  % Red for Fast
    0.5 0.5 0.5  % Gray for Unknown
];

colors_participants = [ % Colors for each participant
        0.2, 0.7, 0.8;  % Blueish
        0.3, 0.7, 0.4;  % Greenish
        0.6, 0.3, 0.7;  % Purple
        0.9, 0.4, 0.3;  % Yellow-orange
        0.1, 0.4, 0.6;  % Dark teal
        0.8, 0.6, 0.1;  % Brown-orange
        0.4, 0.6, 0.2;  % Olive-green
        0.2, 0.3, 0.7;  % Dark blue
        0.7, 0.2, 0.4;  % Reddish-pink
        0.1, 0.6, 0.5;  % Cyan-like
        0.5, 0.4, 0.7;  % Muted violet
        0.6, 0.2, 0.5;  % Magenta-like
    ];
colors_participants = min(colors_participants * 1.4, 1); % Increase intensity, limit to 1

accuracy_categories = ["Low", "Medium", "High"]; 
balance_categories = ["None", "Medium", "High"];
speed_categories = ["Slow", "Medium", "Fast"];

%% Labels the different speeds with balance = 0, and accuracy = 2; Step error and baselines

baseline_errors = struct();

for i=1:length(data_names) % go through all participants
    bmh = string(data_names{i});
    [speed_s0a2b0, speed_s1a2b0, speed_s2a2b0] = BaselineErrors(imported_data,bmh,n,section);
    baseline_errors.(bmh).speed_s0a2b0 = speed_s0a2b0;
    baseline_errors.(bmh).speed_s1a2b0 = speed_s1a2b0;
    baseline_errors.(bmh).speed_s2a2b0 = speed_s2a2b0;
end

%% Error classification

error = [];

for i=1:size(data_names,2)

    bmh = data_names{i};
    numtrials(i) = numel(fieldnames(imported_data.(bmh).matches));
    %dataset_name = eval(data_names{i});

   for k=1:numtrials(i)  % iterate through trials and add data to the struct
        
        trialname = sprintf('Trial%04d', k);
    
        if trialname == "Trial0004"
            continue
        end
    
        % Extract error values for the current trial (distance to the line)
        if section ==1
            errorData_l = imported_data.(bmh).matches.(trialname).left.match(:,n);   
            errorData_r = imported_data.(bmh).matches.(trialname).right.match(:,n);
        elseif section==2
            errorData_l = imported_data.(bmh).matches.(trialname).left.match_straight(:,n);   
            errorData_r = imported_data.(bmh).matches.(trialname).right.match_straight(:,n);
        elseif section==3
            errorData_l = imported_data.(bmh).matches.(trialname).left.match_curve(:,n);   
            errorData_r = imported_data.(bmh).matches.(trialname).right.match_curve(:,n);
        else
            error("section must be a value between 1 and 3")
        end 
    
        speed = imported_data.(bmh).prompt.speed(k,1);
        accuracy = imported_data.(bmh).prompt.accuracy(k,1);
        balance = imported_data.(bmh).prompt.balance(k,1);
       
        % Combine left and right errors and add to struct
        if foot==1
            errorData = [errorData_l; errorData_r];   %To use from BOTH feet
        elseif foot==2
            errorData = errorData_r;                  %To use only RIGHT foot
        elseif foot==3
            errorData = errorData_l;                  %To use only LEFT foot
        else
            error("foot must be a number between 1 and 3")
        end
        
        name = trialname + "_" + string(bmh); 

        error.byspeed.(speed).(name) = errorData;
        error.byaccuracy.(accuracy).(name) = errorData;
        error.bybalance.(balance).(name) = errorData;
    end
end
close all %to not display the steperror graph

%% Sort data by label
participants = "multiple";
[data_labels, data_ordered_bylabel, surveys] = Sort_bylabel(imported_data, metabolics, baseline_errors, participants);

%% Create BoxPlots
if plot1 == 1
    BoxPlot_error(imported_data,error,n,section,foot,colors,subtract_speeds_1,subtract_all_1)
end

%% Filtered Bar Plot for Accuracy Condition (Balance == 0 & Trials 3-29) == Accuracy=*bars*, Balance=0, Speed=3
if plot2 == 1
    Accuracy_bars_error(imported_data,error,n,section,foot,colors,subtract_speeds_2,subtract_all_2,accuracy_categories,speed_input_2, balance_input_2,speed_s0a2b0,speed_s1a2b0,speed_s2a2b0)
end

%% Customized multiple plots
if plot4 == 1
    Multiple_customized(imported_data,error,n,section,foot,speed_input,accuracy_input,balance_input,subtract_speeds_4,subtract_all_4)
end

%% Metabolic data for this participant across all trials
if plot5 == 1
    Metabolics(imported_data,error,colors,metabolics,baseline_errors,normalize)
end

%% Plot all conditions according to different parameters
if plot6 == 1
    plot_number_max = 43;
    if color_shape_option_1_1 == 0 && color_shape_option_1_2 == 1
        error("Choose what combination to plot with color_shape_option_1_1")
    end
    if color_shape_option_1_2 == 1
        for i=1:plot_number_max
            plot_number = i;
            All_conditions(data_labels,plot_number,subtract_all_6,subtract_speeds_6,baseline_errors,show_cluster_analysis1,participants,show_line,show_formula1,color_shape_option_1_1,normalize_metabolics1)
        end
    elseif color_shape_option_1_1 == 0
        for i=1:3
            All_conditions(data_labels,plot_number,subtract_all_6,subtract_speeds_6,baseline_errors,show_cluster_analysis1,participants,show_line,show_formula1,i,normalize_metabolics1)
        end
    else
        All_conditions(data_labels,plot_number,subtract_all_6,subtract_speeds_6,baseline_errors,show_cluster_analysis1,participants,show_line,show_formula1,color_shape_option_1_1,normalize_metabolics1)
    end
end

%% Participants' performance, answers and prompts
if plot7 == 1
    Performance_prompts(data_labels,imported_data,baseline_errors,colors,subtract_all_7,subtract_speeds_7,participants,shaded_plot,plot_eachparticipant,colors_participants,metabolics)
end

%% Compare baselines
if plot8 == 1
    Compare_baselines(imported_data,error,baseline_errors)
end

%% Survey responses
if plot9 == 1
    Survey_responses(data_labels, color_shape_option_2, participants, show_cluster_analysis2, show_formula2,normalize_metabolics2)
end

%% Pairwise correlation matrix
if plot10 == 1
    Correlation_matrix(data_ordered_bylabel, participants, option_correlations,show_all_plots)
end
