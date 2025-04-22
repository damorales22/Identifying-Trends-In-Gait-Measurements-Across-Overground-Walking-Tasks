clear 
clc
close all

%% User variables and inputs

filepath = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\BMH10\";
%filepath = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\BMH02\";
%filepath = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\BMH08\";
%filepath = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\BMH20\";
%filepath = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\BMH06\";
%filepath = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\BMH19\";

n = 3; %Column 3 LONGITUDINAL ERROR (distance to the infinite line), Column 4 SIGNED LONGITUDINAL ERROR, Column 5 TOTAL ERROR (distance to the red tape), Column 6 LATERAL ERROR; (of matches.side.match)
section = 2; %1 for ALL, 2 for only STRAIGHT, 3 for only CURVES
foot = 1; %1 for BOTH FEET analysis (combined), 2 for only RIGHT foot, 3 for only LEFT foot

plot1 = 0; %Errors grouped by condition and overall bar plot   
    subtract_all_1 = 0;    %0 do nothing, 1 to subtract the error mean across all trials to the error bars. ONLY applied when using the SIGNED error.
    subtract_speeds_1 = 0; %0 do nothing, 1 to subtract the error mean according to each speed to the error bars. ONLY applied when using the SIGNED error.

plot2 = 0; %Bar Plot for Accuracy Conditions (3; any) and a given Balance and Speed
    speed_input_2 = 0;   % 0=slow, 1=medium, 2=fast, 3=any
    balance_input_2 = 0; % 0=none, 1=medium, 2=high, 3=any; (disturbances)

    % Substracting the means do not make sense in this plot, as the means are JUST for a=2
    subtract_all_2 = 0;    %0 do nothing, 1 to subtract the error mean across all trials to the error bars. ONLY applied when using the ABOSLUTE or SIGNED error and ONLY when grouped by speed.
    subtract_speeds_2 = 0; %0 do nothing, 1 to subtract the error mean according to each speed to the error bars. ONLY applied when using the ABOSLUTE or SIGNED error and ONLY when grouped by speed.

plot3 = 0; %Step error for every single step for just a specific trial
    trial_analyze = "Trial0024";
    subtract_all_3 = 0;    %0 do nothing, 1 to subtract the error mean ACROSS ALL TRIALS to the error bars. ONLY applied when using the SIGNED error
    subtract_speeds_3 = 0; %0 do nothing, 1 to subtract the error mean according to EACH SPEED to the error bars. ONLY applied when using the SIGNED error

plot4 = 1; %Customized multiple plots of the given conditions
    speed_input = 3;    % 0=slow, 1=medium, 2=fast, 3=any
    accuracy_input = 2; % 0=none, 1=medium, 2=high, 3=any
    balance_input = 0;  % 0=none, 1=medium, 2=high, 3=any; (disturbances)

    subtract_all_4 = 0;    %0 do nothing, 1 to subtract the error mean ACROSS ALL TRIALS to the error bars. ONLY applied when using the SIGNED error
    subtract_speeds_4 = 0; %0 do nothing, 1 to subtract the error mean according to EACH SPEED to the error bars. ONLY applied when using the SIGNED error

plot5 = 0; % Plot the metabolic data for this participant
    normalize = 0; %Normalize to the weigth of each participant

plot6 = 0; % Plot all conditions
    plot_number = 18; % 1: speed vs metabolics.      2: speed vs step width.        3: speed vs step length.        4: speed vs step error.
                     % 5: metabolics vs step width. 6: metabolics vs step length.  7: metabolics vs step error.    8: step width vs step length. 
                     % 9: step width vs step error. 10: step length vs step error. 11: speed vs width variability. 12: metabolics vs width variability. 
                     % 13: step error vs width variability. 14: step error vs head angle. 15: step error vs head angle (just colors)
                     % 16: step error vs head angle (shape and colors). 17: step error vs trunk angle. 18: step error vs trunk angle (shape and colors). 
                     % 19: speed vs stride time. 20: trunk angle vs step width. 21: trunk angle vs step width var. 22: stride time vs step length var. 
                     % 23: stride time var vs step length var. 24: stride time var vs step error. 25: step width var vs step length. 
                     % 26: step width var vs step length var. 27: step width var vs step width 28: step width var vs metabolics.
                     % 29: step length vs metabolics. 30: step length var vs metabolics. 31: step width vs metabolics.

    color_shape_option_1 = 1; % 1 - Color: accuracy. Shape: balance. 2 - Color: balance. Shape: speed. 3 - Color: speed. Shape: accuracy
    show_cluster_analysis1 = 0; % Show clustering plots
    show_line = 1; % Plot regression line
        show_formula1 = 0;
    normalize_metabolics1 = 0;
    subtract_all_6 = 0;    % Only when using the step error
    subtract_speeds_6 = 0; % Only when using the step error

plot7 = 0; % Participants' performance and propmts 
    shaded_plot = 0;
    plot_eachparticipant = 1; % Points in the bar plots representing each participant
    subtract_all_7 = 0;    % Only when using the step error
    subtract_speeds_7 = 0; % Only when using the step error. A baseline for each speed (3)

plot8 = 0; % Compare the baselines; the only output is the txt

plot9 = 0; % Survey responses
    color_shape_option_2 = 1; % 1 - Color: accuracy. Shape: balance. 2 - Color: balance. Shape: speed. 3 - Color: speed. Shape: accuracy
    show_cluster_analysis2 = 0;
    normalize_metabolics1 = 0;
    show_formula2 = 0;

plot10 = 0; % Pairwise correlation matrix
    option_correlations = 1; % Choose what correlations to show. 1: all. 2: close to +1. 3: close to -1. 4: close to 0
    show_all_plots = 0; % 0 to only show the main one, 1 to show also the 9 plots (one for each condition). The excel file is generated anyway

%% Other variables

bmh = regexp(filepath + "MC\", 'BMH\d{2}', 'match', 'once');
try
    imported = load(filepath + "MC\" + "data.mat");
catch 
    imported = load(filepath + "MC\" + "data_" + bmh + ".mat");
end
imported_data.(bmh) = imported.data;
numtrials = size(imported_data.(bmh).walkingspeed,2); %Number of trials

metabolics.(bmh) = load(filepath + "Metabolics.csv");

colors = [
    0.1 0.6 0.1; % Green for Slow
    0.1 0.1 0.8; % Blue for Medium
    0.8 0.1 0.1  % Red for Fast
    0.5 0.5 0.5  % Gray for Unknown
];

accuracy_categories = ["Low", "Medium", "High"]; 
balance_categories = ["None", "Medium", "High"];
speed_categories = ["Slow", "Medium", "Fast"];

%% Labels the different speeds with balance = 0, and accuracy = 2

baseline_errors = struct();

[speed_s0a2b0, speed_s1a2b0, speed_s2a2b0] = BaselineErrors(imported_data,bmh,n,section);
baseline_errors.(bmh).speed_s0a2b0 = speed_s0a2b0;
baseline_errors.(bmh).speed_s1a2b0 = speed_s1a2b0;
baseline_errors.(bmh).speed_s2a2b0 = speed_s2a2b0;

clear imported;

%% Labels the different speeds with balance = 0, and accuracy = 2
for i = 3:29

    trialname = sprintf('Trial%04d', i);
    
    speed_val = str2double(imported_data.(bmh).prompt.speed(i, 2));
    accuracy_val = str2double(imported_data.(bmh).prompt.accuracy(i, 2));
    balance_val = str2double(imported_data.(bmh).prompt.balance(i, 2));

    if speed_val == 0 && accuracy_val == 2 && balance_val == 0
        label_s0a2b0 = trialname;
        if n == 3 % longitudinal - absolute
            if section == 1 
                speed_s0a2b0 = imported_data.(bmh).matches.(trialname).overall.overall_mean_absolute; %All
            elseif section == 2
                speed_s0a2b0 = imported_data.(bmh).matches.(trialname).overall.straight_mean_absolute; %Straight
            elseif section == 3
                speed_s0a2b0 = imported_data.(bmh).matches.(trialname).overall.curve_mean_absolute; %Curves
            end
        elseif n == 4 % longitudinal - signed
            if section == 1
                speed_s0a2b0 = imported_data.(bmh).matches.(trialname).overall.overall_mean_signed;
            elseif section == 2
                speed_s0a2b0 = imported_data.(bmh).matches.(trialname).overall.straight_mean_signed;
            elseif section == 3
                speed_s0a2b0 = imported_data.(bmh).matches.(trialname).overall.curve_mean_signed;
            end
        end
    elseif speed_val == 1 && accuracy_val == 2 && balance_val == 0
        label_s1a2b0 = trialname;
        if n == 3 % longitudinal - absolute
            if section == 1
                speed_s1a2b0 = imported_data.(bmh).matches.(trialname).overall.overall_mean_absolute;
            elseif section == 2
                speed_s1a2b0 = imported_data.(bmh).matches.(trialname).overall.straight_mean_absolute;
            elseif section == 3
                speed_s1a2b0 = imported_data.(bmh).matches.(trialname).overall.curve_mean_absolute;
            end
        elseif n == 4 % longitudinal - signed
            if section == 1
                speed_s1a2b0 = imported_data.(bmh).matches.(trialname).overall.overall_mean_signed;
            elseif section == 2
                speed_s1a2b0 = imported_data.(bmh).matches.(trialname).overall.straight_mean_signed;
            elseif section == 3
                speed_s1a2b0 = imported_data.(bmh).matches.(trialname).overall.curve_mean_signed;
            end
        end
    elseif speed_val == 2 && accuracy_val == 2 && balance_val == 0
        label_s2a2b0 = trialname;
        if n == 3 % longitudinal - absolute
            if section == 1
                speed_s2a2b0 = imported_data.(bmh).matches.(trialname).overall.overall_mean_absolute;
            elseif section == 2
                speed_s2a2b0 = imported_data.(bmh).matches.(trialname).overall.straight_mean_absolute;
            elseif section == 3
                speed_s2a2b0 = imported_data.(bmh).matches.(trialname).overall.curve_mean_absolute;
            end
        elseif n == 4 % longitudinal - signed
            if section == 1
                speed_s2a2b0 = imported_data.(bmh).matches.(trialname).overall.overall_mean_signed;
            elseif section == 2
                speed_s2a2b0 = imported_data.(bmh).matches.(trialname).overall.straight_mean_signed;
            elseif section == 3
                speed_s2a2b0 = imported_data.(bmh).matches.(trialname).overall.curve_mean_signed;
            end
        end
    end
end

%% Get Data

error = [];
error.byspeed = [];
error.byaccuracy = [];
error.bybalance = [];

k=1;

% iterate through trials and add data to the struct
while k <= numtrials
    
    trialname = sprintf('Trial%04d', k);

    if trialname == "Trial0004"
        k = k+1;
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

    % Ensure the fields exists in the struct
    if ~isfield(error.byspeed, speed)
        error.byspeed.(speed) = struct();
    end
    if ~isfield(error.byaccuracy, accuracy)
        error.byaccuracy.(accuracy) = struct();
    end
    if ~isfield(error.bybalance, balance)
        error.bybalance.(balance) = struct();
    end

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
    
    error.byspeed.(speed).(trialname) = errorData;
    error.byaccuracy.(accuracy).(trialname) = errorData;
    error.bybalance.(balance).(trialname) = errorData;

    k=k+1;
end
close all %to not display the steperror graph

%% Sort data by label
participants = "single";
[data_labels, data_ordered_bylabel, surveys] = Sort_bylabel(imported_data, metabolics, baseline_errors, participants, bmh);

%% Create BoxPlots
if plot1 == 1
    BoxPlot_error(imported_data,error,n,section,foot,colors,subtract_speeds_1,subtract_all_1, bmh)
end

%% Filtered Bar Plot for Accuracy Condition (Balance == 0 & Trials 3-29) == Accuracy=*bars*, Balance=0, Speed=3
if plot2 == 1
    Accuracy_bars_error(imported_data,error,n,section,foot,colors,subtract_speeds_2,subtract_all_2,accuracy_categories,speed_input_2, balance_input_2,speed_s0a2b0,speed_s1a2b0,speed_s2a2b0, bmh)
end

%% Customized multiple plots
if plot4 == 1
    Multiple_customized(imported_data,error,n,section,foot,speed_input,accuracy_input,balance_input,subtract_speeds_4,subtract_all_4, bmh)
end

%% Metabolic data for this participant across all trials
if plot5 == 1
    Metabolics(imported_data,error,colors,metabolics,baseline_errors,normalize, bmh)
end

%% Plot all conditions according to different parameters
if plot6 == 1
    All_conditions(data_labels,plot_number,subtract_all_6,subtract_speeds_6,baseline_errors,show_cluster_analysis1,participants,show_line,show_formula1,color_shape_option_1, bmh)
end

%% Participants' performance, answers and prompts
if plot7 == 1
    Performance_prompts(data_labels,imported_data,baseline_errors,colors,subtract_all_7,subtract_speeds_7,participants,shaded_plot,plot_eachparticipant,colors_participants,bmh)
end

%% Compare baselines
if plot8 == 1
    Compare_baselines(imported_data,error,baseline_errors, bmh)
end

%% Survey responses
if plot9 == 1
    Survey_responses(data_labels, color_shape_option_2, participants, show_cluster_analysis2, show_formula2, bmh)
end

%% Pairwise correlation matrix
if plot10 == 1
    Correlation_matrix(data_ordered_bylabel, participants, option_correlations,show_all_plots, bmh)
end