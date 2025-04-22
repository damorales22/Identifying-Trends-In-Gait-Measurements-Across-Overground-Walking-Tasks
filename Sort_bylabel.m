function [data_labels, data_ordered_bylabel, surveys] = Sort_bylabel(imported_data, metabolic_data, baseline_errors, participants,bmh)

    path = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\";
    
    bmh_list = string(fieldnames(imported_data));

    filename = 'C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\Correlations.xlsx';
    filename2 = 'C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\Correlations_with_surveys.xlsx';
    file_path = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\Survey and Prompts.xlsx";

    %% Surveys
    surveys = struct();

    % Transform difficulties to numbers
    target1 = "Extremly easy";
    target2 = "Somewhat easy";
    target3 = "Neither easy nor difficult";
    target4 = "Somewhat difficult";
    target5 = "Extremely difficult";

    for i = 1:length(bmh_list)
        bmh = bmh_list{i};
        range = 'B4:J30';
        data = readtable(file_path, 'Sheet', bmh, 'Range', range);
        
        surveys.balance_struct.not_normalized.(bmh) = data{:, 1};
        surveys.foot_placement_struct.not_normalized.(bmh) = data{:, 2};
        surveys.walking_speed_struct.not_normalized.(bmh) = data{:, 3};
        surveys.metabolics_struct.not_normalized.(bmh) = data{:, 4};
        
        surveys.balance_struct.normalized.(bmh) = data{:, 5};
        surveys.foot_placement_struct.normalized.(bmh) = data{:, 6};
        surveys.walking_speed_struct.normalized.(bmh) = data{:, 7};
        surveys.metabolics_struct.normalized.(bmh) = data{:, 8};

        %surveys.difficulty_struct.(bmh) = data{:, 9};
        difficulties = data{:, 9};
        
        %for j=1:length(surveys.difficulty_struct.(bmh))
        for j=1:length(difficulties)
            text_difficulty = string(data{j,9});

            if strcmp(text_difficulty, target1) || editDistance(text_difficulty, target1) <= 2 %if there is a 2 letters or less difference between the two texts
                numeric_val = 1;
            elseif strcmp(text_difficulty, target2) || editDistance(text_difficulty, target2) <= 2
                numeric_val = 2;
            elseif strcmp(text_difficulty, target3) || editDistance(text_difficulty, target3) <= 2
                numeric_val = 3;
            elseif strcmp(text_difficulty, target4) || editDistance(text_difficulty, target4) <= 2
                numeric_val = 4;
            elseif strcmp(text_difficulty, target5) || editDistance(text_difficulty, target5) <= 2
                numeric_val = 5;
            end
            surveys.difficulty_struct.(bmh)(j,1) = numeric_val;
        end
    end
    
    %% Sort by label
    data_labels.s0a0b0 = [];
    data_labels.s0a0b1 = [];
    data_labels.s0a0b2 = [];
    data_labels.s0a1b0 = [];
    data_labels.s0a1b1 = [];
    data_labels.s0a1b2 = [];
    data_labels.s0a2b0 = [];
    data_labels.s0a2b1 = [];
    data_labels.s0a2b2 = [];
    
    data_labels.s1a0b0 = [];
    data_labels.s1a0b1 = [];
    data_labels.s1a0b2 = [];
    data_labels.s1a1b0 = [];
    data_labels.s1a1b1 = [];
    data_labels.s1a1b2 = [];
    data_labels.s1a2b0 = [];
    data_labels.s1a2b1 = [];
    data_labels.s1a2b2 = [];
    
    data_labels.s2a0b0 = [];
    data_labels.s2a0b1 = [];
    data_labels.s2a0b2 = [];
    data_labels.s2a1b0 = [];
    data_labels.s2a1b1 = [];
    data_labels.s2a1b2 = [];
    data_labels.s2a2b0 = [];
    data_labels.s2a2b1 = [];
    data_labels.s2a2b2 = [];
    
    metabolics.normalized = [];
    metabolics.not_normalized = [];
    walking_speed = [];
    step_width = [];
    step_length = [];
    step_error = [];
    survey_balance = [];
    survey_foot_placement = [];
    survey_walking_speed = [];
    survey_metabolics = [];
    survey_difficulty = [];
    age = [];
    height = [];
    weight = [];
    
    if participants == "single"
        metabolics.not_normalized = metabolic_data.(bmh); % W
        metabolics.normalized = metabolic_data.(bmh)/imported_data.(bmh).participant_data.weight; % W/kg
        walking_speed = imported_data.(bmh).walkingspeed; % m/s
        walking_speed_std = imported_data.(bmh).walkingspeed_std; % m/s
        step_width = imported_data.(bmh).meanwidthstraights/10; % cm
        step_width_var = imported_data.(bmh).straightsvariability/10; % cm  
        step_length = imported_data.(bmh).meanlength_straights/10; % cm
        step_length_var = imported_data.(bmh).meanlength_straights_variability/10; % cm
        stride_time = imported_data.(bmh).stride_time; % s
        stride_time_var = imported_data.(bmh).stride_std; % s
        step_error = imported_data.(bmh).overall_errors.straight_mean_absolute/10; % cm
        step_error_std = imported_data.(bmh).overall_errors.straight_std_absolute/10; % cm
        step_error_signed = imported_data.(bmh).overall_errors.straight_mean_signed/10; % cm
        step_error_signed_std = imported_data.(bmh).overall_errors.straight_std_signed/10; % cm
        
        survey_balance = surveys.balance_struct.normalized.(bmh)'; 
        survey_foot_placement = surveys.foot_placement_struct.normalized.(bmh)'; 
        survey_walking_speed = surveys.walking_speed_struct.normalized.(bmh)'; 
        survey_metabolics = surveys.metabolics_struct.normalized.(bmh)'; 
        survey_difficulty = surveys.difficulty_struct.(bmh)'; 

        age = imported_data.(bmh).participant_data.age;  
        height = imported_data.(bmh).participant_data.height;  
        weight = imported_data.(bmh).participant_data.weight;  
     
        mean_head_angle = nan(35, 1);
        mean_trunk_angle = nan(35, 1);
        mean_stj_headset_angle = nan(35, 1);
        mean_dist_gaze = nan(35, 1);
        for j=3:29 % For each bmh
            trialname = sprintf('Trial%04d', j);
            mean_head_angle(j,1) = real(mean(imported_data.(bmh).head_angles.(trialname)(240:end-240), 'omitnan')); % Do not take into account the first and last 2 seconds
            mean_trunk_angle(j,1) = imported_data.(bmh).trunk_angles.(trialname); % Do not take into account the first and last 2 seconds and remove NaN
            mean_stj_headset_angle(j,1) = mean(imported_data.(bmh).stj_headset_angles.(trialname)(240:end-240), 'omitnan'); % Do not take into account the first and last 2 seconds       
            mean_dist_gaze(j,1) = imported_data.(bmh).dist_gaze.(trialname); % Already in cm
        end
        head_angle(i,:) = mean_head_angle;
        trunk_angle(i,:) = mean_trunk_angle;
        stj_headset_angle(i,:) = mean_stj_headset_angle;
        dist_gaze(i,:) = mean_dist_gaze;
    
    elseif participants == "multiple"
        for i=1:size(fieldnames(imported_data),1) % number of bmhs
            names = string(fieldnames(imported_data));
            bmh = names(i);
            metabolics.not_normalized(i,:) = table2array(metabolic_data.(bmh)); % W
            metabolics.normalized(i,:) = table2array(metabolic_data.(bmh))/imported_data.(bmh).participant_data.weight; % W/kg
            walking_speed(i,:) = imported_data.(bmh).walkingspeed; % m/s
            walking_speed_std(i,:) = imported_data.(bmh).walkingspeed_std; % m/s
            step_width(i,:) = imported_data.(bmh).meanwidthstraights/10; % cm
            step_width_var(i,:) = imported_data.(bmh).straightsvariability/10; % cm
            step_length(i,:) = imported_data.(bmh).meanlength_straights/10; % cm
            step_length_var(i,:) = imported_data.(bmh).meanlength_straights_variability/10; % cm
            stride_time(i,:) = imported_data.(bmh).stride_time; % s
            stride_time_var(i,:) = imported_data.(bmh).stride_std; % a
            step_error(i,:) = imported_data.(bmh).overall_errors.straight_mean_absolute/10; % cm
            step_error_std(i,:) = imported_data.(bmh).overall_errors.straight_std_absolute/10; % cm
            step_error_signed(i,:) = imported_data.(bmh).overall_errors.straight_mean_signed/10; % cm
            step_error_signed_std(i,:) = imported_data.(bmh).overall_errors.straight_std_signed/10; % cm
    
            survey_balance(i,:) = surveys.balance_struct.normalized.(bmh)'; 
            survey_foot_placement(i,:) = surveys.foot_placement_struct.normalized.(bmh)'; 
            survey_walking_speed(i,:) = surveys.walking_speed_struct.normalized.(bmh)'; 
            survey_metabolics(i,:) = surveys.metabolics_struct.normalized.(bmh)'; 
            survey_difficulty(i,:) = surveys.difficulty_struct.(bmh); 

            age(i,:) = imported_data.(bmh).participant_data.age;  
            height(i,:) = imported_data.(bmh).participant_data.height;  
            weight(i,:) = imported_data.(bmh).participant_data.weight;  
    
            mean_head_angle = nan(35, 1);
            mean_trunk_angle = nan(35, 1);
            mean_stj_headset_angle = nan(35, 1);
            mean_dist_gaze = nan(35, 1);
            for j=3:29 % For all trials of each bmh
                trialname = sprintf('Trial%04d', j);
                mean_head_angle(j,1) = real(mean(imported_data.(bmh).head_angles.(trialname)(240:end-240), 'omitnan')); % Do not take into account the first and last 2 seconds and remove NaN
                mean_trunk_angle(j,1) = imported_data.(bmh).trunk_angles.(trialname); % Do not take into account the first and last 2 seconds and remove NaN
                mean_stj_headset_angle(j,1) = mean(imported_data.(bmh).stj_headset_angles.(trialname)(240:end-240), 'omitnan'); % Do not take into account the first and last 2 seconds                  
                mean_dist_gaze(j,1) = imported_data.(bmh).dist_gaze.(trialname); % Already in cm
            end
            head_angle(i,:) = mean_head_angle;
            trunk_angle(i,:) = mean_trunk_angle;
            stj_headset_angle(i,:) = mean_stj_headset_angle; 
            dist_gaze(i,:) = mean_dist_gaze;
        end
    end
    
    list_labels = fieldnames(data_labels); %labels stored in the struct
    bmh_list = string(fieldnames(baseline_errors));
    
    for i=1:size(fieldnames(imported_data),1) %through bmhs
        names = string(fieldnames(imported_data));
        bmh = names(i);
        numtrials = size(imported_data.(bmh).labels,2);
        labels = imported_data.(bmh).labels;  % Format: 'sxaybz'
    
        for j = 3:29 %through trials
            current_label = labels{j}; % from the ordered label array
            trialname = sprintf('Trial%04d', j); 
            subtracted = 0;
            for k=1:size(list_labels,1)
                label_struct = list_labels{k};
                if strcmp(current_label, label_struct)
                    data_labels.(current_label).walking_speed(i,:) = walking_speed(i,j);
                    data_labels.(current_label).walking_speed_std(i,:) = walking_speed_std(i,j);
                    data_labels.(current_label).step_width(i,:) = step_width(i,j);
                    data_labels.(current_label).step_width_var(i,:) = step_width_var(i,j);
                    data_labels.(current_label).step_length(i,:) = step_length(i,j);
                    data_labels.(current_label).step_length_var(i,:) = step_length_var(i,j);
                    data_labels.(current_label).stride_time(i,:) = stride_time(i,j);
                    data_labels.(current_label).stride_time_var(i,:) = stride_time_var(i,j);
                    data_labels.(current_label).step_error(i,:) = step_error(i,j);
                    data_labels.(current_label).step_error_std(i,:) = step_error_std(i,j);
                    data_labels.(current_label).step_error_signed(i,:) = step_error_signed(i,j);
                    data_labels.(current_label).step_error_signed_std(i,:) = step_error_signed_std(i,j);
                    data_labels.(current_label).head_angle(i,:) = head_angle(i,j);
                    data_labels.(current_label).trunk_angle(i,:) = trunk_angle(i,j);
                    data_labels.(current_label).metabolics.not_normalized(i,:) = metabolics.not_normalized(i,j);
                    data_labels.(current_label).metabolics.normalized(i,:) = metabolics.normalized(i,j);
                    data_labels.(current_label).stj_headset_angle(i,:) = stj_headset_angle(i,j);
                    data_labels.(current_label).dist_gaze(i,:) = dist_gaze(i,j);
                    data_labels.(current_label).age(i,:) = age(i,1);
                    data_labels.(current_label).height(i,:) = height(i,1);
                    data_labels.(current_label).weight(i,:) = weight(i,1);
    
                    if subtracted == 0
                        j = j-2;
                        subtracted = 1; 
                    end
                    data_labels.(current_label).survey_balance(i,:) = survey_balance(i,j);
                    data_labels.(current_label).survey_foot_placement(i,:) = survey_foot_placement(i,j);
                    data_labels.(current_label).survey_walking_speed(i,:) = survey_walking_speed(i,j);
                    data_labels.(current_label).survey_metabolics(i,:) = survey_metabolics(i,j);
                    data_labels.(current_label).survey_difficulty(i,:) = survey_difficulty(i,j);
                end
            end
        end
    end
    %% Store each variable sorted by label number (s0a0b0,s0a0b1...)
    
    data_ordered_bylabel = struct();
    
    for i=1:size(list_labels,1)
        current_label = list_labels{i};
    
        data_ordered_bylabel.walking_speed(i,:) = data_labels.(current_label).walking_speed;
        data_ordered_bylabel.walking_speed_std(i,:) = data_labels.(current_label).walking_speed_std;
        data_ordered_bylabel.step_width(i,:) = data_labels.(current_label).step_width;
        data_ordered_bylabel.step_width_var(i,:) = data_labels.(current_label).step_width_var;
        data_ordered_bylabel.step_length(i,:) = data_labels.(current_label).step_length;
        data_ordered_bylabel.step_length_var(i,:) = data_labels.(current_label).step_length_var;
        data_ordered_bylabel.stride_time(i,:) = data_labels.(current_label).stride_time;
        data_ordered_bylabel.stride_time_var(i,:) = data_labels.(current_label).stride_time_var;
        data_ordered_bylabel.step_error(i,:) = data_labels.(current_label).step_error;
        data_ordered_bylabel.step_error_std(i,:) = data_labels.(current_label).step_error_std;
        data_ordered_bylabel.step_error_signed(i,:) = data_labels.(current_label).step_error_signed;
        data_ordered_bylabel.step_error_signed_std(i,:) = data_labels.(current_label).step_error_signed_std;
        data_ordered_bylabel.head_angle(i,:) = data_labels.(current_label).head_angle;
        data_ordered_bylabel.trunk_angle(i,:) = data_labels.(current_label).trunk_angle;
        data_ordered_bylabel.metabolics.not_normalized(i,:) = data_labels.(current_label).metabolics.not_normalized;
        data_ordered_bylabel.metabolics.normalized(i,:) = data_labels.(current_label).metabolics.normalized;
        data_ordered_bylabel.dist_gaze(i,:) = data_labels.(current_label).dist_gaze;
    
        % Surveys
        data_ordered_bylabel.surveys_balance(i,:) = data_labels.(current_label).survey_balance;
        data_ordered_bylabel.surveys_foot_placement(i,:) = data_labels.(current_label).survey_foot_placement; 
        data_ordered_bylabel.surveys_walking_speed(i,:) = data_labels.(current_label).survey_walking_speed; 
        data_ordered_bylabel.surveys_metabolics(i,:) = data_labels.(current_label).survey_metabolics; 
        data_ordered_bylabel.surveys_difficulty(i,:) = data_labels.(current_label).survey_difficulty; 

        data_ordered_bylabel.age(i,:) = data_labels.(current_label).age;
        data_ordered_bylabel.height(i,:) = data_labels.(current_label).height;
        data_ordered_bylabel.weight(i,:) = data_labels.(current_label).weight;
    
    end
end
