%function [] = Correlation_matrix(imported_data,error,baseline_errors,metabolic_data,bmh)
function [] = Correlation_matrix(data_ordered_bylabel, participants, option_correlations, show_all_plots, bmh)

    if participants == "multiple"
        bmh_text = "Across all Participants";
    else
        bmh_text = string(bmh);
    end

    filename = 'C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\Correlations.xlsx';
    filename2 = 'C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\Correlations_with_surveys.xlsx';

    %% Analysis

    positions_s0 = [1,2,3,4,5,6,7,8,9]; % Positions of the s0 labels in the arrays
    positions_s1 = [10,11,12,13,14,15,16,17,18];
    positions_s2 = [19,20,21,22,23,24,25,26,27];
    positions_a0 = [1,2,3,10,11,12,19,20,21];
    positions_a1 = [4,5,6,13,14,15,22,23,24];
    positions_a2 = [7,8,9,16,17,18,25,26,27];
    positions_b0 = [1,4,7,10,13,16,19,22,25];
    positions_b1 = [2,5,8,11,14,17,20,23,26];
    positions_b2 = [3,6,9,12,15,18,21,24,27];

    walking_speed_byprompt = struct();
    walking_speed_std_byprompt = struct();
    metabolics_byprompt = struct();
    step_width_byprompt = struct();
    step_width_var_byprompt = struct();
    step_length_byprompt = struct();
    step_length_var_byprompt = struct();
    stride_time_byprompt = struct();
    stride_time_var_byprompt = struct();
    step_error_byprompt = struct();
    step_error_std_byprompt = struct();
    head_angle_byprompt = struct();
    trunk_angle_byprompt = struct();

    survey_balance_byprompt = struct();
    survey_foot_placement_byprompt = struct();     
    survey_walking_speed_byprompt = struct();
    survey_metabolics_byprompt = struct();
    survey_difficulty_byprompt = struct();


    %% Correlations (all) with survey responses 

    % Concatenate all points of the interest labels 
    % walking_speed_vector = data_ordered_bylabel.walking_speed(:);
    % walking_speed_std_vector = data_ordered_bylabel.walking_speed_std(:);
    % metabolics_vector = data_ordered_bylabel.metabolics.normalized(:);
    % step_width_vector = data_ordered_bylabel.step_width(:);
    % step_width_var_vector = data_ordered_bylabel.step_width_var(:);
    % step_length_vector = data_ordered_bylabel.step_length(:);
    % step_length_var_vector = data_ordered_bylabel.step_length_var(:);
    % stride_time_vector = data_ordered_bylabel.stride_time(:);
    % stride_time_var_vector = data_ordered_bylabel.stride_time_var(:);
    % step_error_vector = data_ordered_bylabel.step_error(:);
    % step_error_var_vector = data_ordered_bylabel.step_error_std(:);
    % head_angle_vector = data_ordered_bylabel.head_angle(:);
    % trunk_angle_vector = data_ordered_bylabel.trunk_angle(:);
    % dist_gaze_vector = data_ordered_bylabel.dist_gaze(:);

    % surveys_balance_vector = data_ordered_bylabel.surveys_balance(:);
    % surveys_foot_placement_vector = data_ordered_bylabel.surveys_foot_placement(:);
    % surveys_walking_speed_vector = data_ordered_bylabel.surveys_walking_speed(:);
    % surveys_metabolics_vector = data_ordered_bylabel.surveys_metabolics(:);
    % surveys_difficulty_vector = data_ordered_bylabel.surveys_difficulty(:);

    % Average to eliminate noise
    walking_speed_vector = nanmean(data_ordered_bylabel.walking_speed, 2);
    walking_speed_std_vector = nanmean(data_ordered_bylabel.walking_speed_std, 2);
    metabolics_vector = nanmean(data_ordered_bylabel.metabolics.normalized, 2);
    step_width_vector = nanmean(data_ordered_bylabel.step_width, 2);
    step_width_var_vector = nanmean(data_ordered_bylabel.step_width_var, 2);
    step_length_vector = nanmean(data_ordered_bylabel.step_length, 2);
    step_length_var_vector = nanmean(data_ordered_bylabel.step_length_var, 2);
    stride_time_vector = nanmean(data_ordered_bylabel.stride_time, 2);
    stride_time_var_vector = nanmean(data_ordered_bylabel.stride_time_var, 2);
    step_error_vector = nanmean(data_ordered_bylabel.step_error, 2);
    step_error_var_vector = nanmean(data_ordered_bylabel.step_error_std, 2);
    head_angle_vector = nanmean(data_ordered_bylabel.head_angle, 2);
    trunk_angle_vector = nanmean(data_ordered_bylabel.trunk_angle, 2);
    dist_gaze_vector = nanmean(data_ordered_bylabel.dist_gaze, 2);

    surveys_balance_vector = mean(data_ordered_bylabel.surveys_balance, 2);
    surveys_foot_placement_vector = mean(data_ordered_bylabel.surveys_foot_placement, 2);
    surveys_walking_speed_vector = mean(data_ordered_bylabel.surveys_walking_speed, 2);
    surveys_metabolics_vector = mean(data_ordered_bylabel.surveys_metabolics, 2);
    surveys_difficulty_vector = mean(data_ordered_bylabel.surveys_difficulty, 2);


    % Prompts
    walking_speed_byprompt.s0 = data_ordered_bylabel.walking_speed(positions_s0, :);
    walking_speed_byprompt.s1 = data_ordered_bylabel.walking_speed(positions_s1, :);
    walking_speed_byprompt.s2 = data_ordered_bylabel.walking_speed(positions_s2, :);
    walking_speed_byprompt.a0 = data_ordered_bylabel.walking_speed(positions_a0, :);
    walking_speed_byprompt.a1 = data_ordered_bylabel.walking_speed(positions_a1, :);
    walking_speed_byprompt.a2 = data_ordered_bylabel.walking_speed(positions_a2, :);
    walking_speed_byprompt.b0 = data_ordered_bylabel.walking_speed(positions_b0, :);
    walking_speed_byprompt.b1 = data_ordered_bylabel.walking_speed(positions_b1, :);
    walking_speed_byprompt.b2 = data_ordered_bylabel.walking_speed(positions_b2, :);

    walking_speed_std_byprompt.s0 = data_ordered_bylabel.walking_speed_std(positions_s0, :);
    walking_speed_std_byprompt.s0(any(isnan(walking_speed_std_byprompt.s0), 2), :) = [];  
    walking_speed_std_byprompt.s1 = data_ordered_bylabel.walking_speed_std(positions_s1, :);
    walking_speed_std_byprompt.s1(any(isnan(walking_speed_std_byprompt.s1), 2), :) = [];
    walking_speed_std_byprompt.s2 = data_ordered_bylabel.walking_speed_std(positions_s2, :);
    walking_speed_std_byprompt.s2(any(isnan(walking_speed_std_byprompt.s2), 2), :) = [];
    walking_speed_std_byprompt.a0 = data_ordered_bylabel.walking_speed_std(positions_a0, :);
    walking_speed_std_byprompt.a0(any(isnan(walking_speed_std_byprompt.a0), 2), :) = [];
    walking_speed_std_byprompt.a1 = data_ordered_bylabel.walking_speed_std(positions_a1, :);
    walking_speed_std_byprompt.a1(any(isnan(walking_speed_std_byprompt.a1), 2), :) = [];
    walking_speed_std_byprompt.a2 = data_ordered_bylabel.walking_speed_std(positions_a2, :);
    walking_speed_std_byprompt.a2(any(isnan(walking_speed_std_byprompt.a2), 2), :) = [];
    walking_speed_std_byprompt.b0 = data_ordered_bylabel.walking_speed_std(positions_b0, :);
    walking_speed_std_byprompt.b0(any(isnan(walking_speed_std_byprompt.b0), 2), :) = [];
    walking_speed_std_byprompt.b1 = data_ordered_bylabel.walking_speed_std(positions_b1, :);
    walking_speed_std_byprompt.b1(any(isnan(walking_speed_std_byprompt.b1), 2), :) = [];
    walking_speed_std_byprompt.b2 = data_ordered_bylabel.walking_speed_std(positions_b2, :);
    walking_speed_std_byprompt.b2(any(isnan(walking_speed_std_byprompt.b2), 2), :) = [];

    metabolics_byprompt.s0 = data_ordered_bylabel.metabolics.normalized(positions_s0, :);
    metabolics_byprompt.s1 = data_ordered_bylabel.metabolics.normalized(positions_s1, :);
    metabolics_byprompt.s2 = data_ordered_bylabel.metabolics.normalized(positions_s2, :);
    metabolics_byprompt.a0 = data_ordered_bylabel.metabolics.normalized(positions_a0, :);
    metabolics_byprompt.a1 = data_ordered_bylabel.metabolics.normalized(positions_a1, :);
    metabolics_byprompt.a2 = data_ordered_bylabel.metabolics.normalized(positions_a2, :);
    metabolics_byprompt.b0 = data_ordered_bylabel.metabolics.normalized(positions_b0, :);
    metabolics_byprompt.b1 = data_ordered_bylabel.metabolics.normalized(positions_b1, :);
    metabolics_byprompt.b2 = data_ordered_bylabel.metabolics.normalized(positions_b2, :);

    step_width_byprompt.s0 = data_ordered_bylabel.step_width(positions_s0, :);
    step_width_byprompt.s1 = data_ordered_bylabel.step_width(positions_s1, :);
    step_width_byprompt.s2 = data_ordered_bylabel.step_width(positions_s2, :);
    step_width_byprompt.a0 = data_ordered_bylabel.step_width(positions_a0, :);
    step_width_byprompt.a1 = data_ordered_bylabel.step_width(positions_a1, :);
    step_width_byprompt.a2 = data_ordered_bylabel.step_width(positions_a2, :);
    step_width_byprompt.b0 = data_ordered_bylabel.step_width(positions_b0, :);
    step_width_byprompt.b1 = data_ordered_bylabel.step_width(positions_b1, :);
    step_width_byprompt.b2 = data_ordered_bylabel.step_width(positions_b2, :);

    step_width_var_byprompt.s0 = data_ordered_bylabel.step_width_var(positions_s0, :);
    step_width_var_byprompt.s1 = data_ordered_bylabel.step_width_var(positions_s1, :);
    step_width_var_byprompt.s2 = data_ordered_bylabel.step_width_var(positions_s2, :);
    step_width_var_byprompt.a0 = data_ordered_bylabel.step_width_var(positions_a0, :);
    step_width_var_byprompt.a1 = data_ordered_bylabel.step_width_var(positions_a1, :);
    step_width_var_byprompt.a2 = data_ordered_bylabel.step_width_var(positions_a2, :);
    step_width_var_byprompt.b0 = data_ordered_bylabel.step_width_var(positions_b0, :);
    step_width_var_byprompt.b1 = data_ordered_bylabel.step_width_var(positions_b1, :);
    step_width_var_byprompt.b2 = data_ordered_bylabel.step_width_var(positions_b2, :);

    step_length_byprompt.s0 = data_ordered_bylabel.step_length(positions_s0, :);
    step_length_byprompt.s1 = data_ordered_bylabel.step_length(positions_s1, :);
    step_length_byprompt.s2 = data_ordered_bylabel.step_length(positions_s2, :);
    step_length_byprompt.a0 = data_ordered_bylabel.step_length(positions_a0, :);
    step_length_byprompt.a1 = data_ordered_bylabel.step_length(positions_a1, :);
    step_length_byprompt.a2 = data_ordered_bylabel.step_length(positions_a2, :);
    step_length_byprompt.b0 = data_ordered_bylabel.step_length(positions_b0, :);
    step_length_byprompt.b1 = data_ordered_bylabel.step_length(positions_b1, :);
    step_length_byprompt.b2 = data_ordered_bylabel.step_length(positions_b2, :);

    step_length_var_byprompt.s0 = data_ordered_bylabel.step_length_var(positions_s0, :);
    step_length_var_byprompt.s1 = data_ordered_bylabel.step_length_var(positions_s1, :);
    step_length_var_byprompt.s2 = data_ordered_bylabel.step_length_var(positions_s2, :);
    step_length_var_byprompt.a0 = data_ordered_bylabel.step_length_var(positions_a0, :);
    step_length_var_byprompt.a1 = data_ordered_bylabel.step_length_var(positions_a1, :);
    step_length_var_byprompt.a2 = data_ordered_bylabel.step_length_var(positions_a2, :);
    step_length_var_byprompt.b0 = data_ordered_bylabel.step_length_var(positions_b0, :);
    step_length_var_byprompt.b1 = data_ordered_bylabel.step_length_var(positions_b1, :);
    step_length_var_byprompt.b2 = data_ordered_bylabel.step_length_var(positions_b2, :);

    stride_time_byprompt.s0 = data_ordered_bylabel.stride_time(positions_s0, :);
    stride_time_byprompt.s1 = data_ordered_bylabel.stride_time(positions_s1, :);
    stride_time_byprompt.s2 = data_ordered_bylabel.stride_time(positions_s2, :);
    stride_time_byprompt.a0 = data_ordered_bylabel.stride_time(positions_a0, :);
    stride_time_byprompt.a1 = data_ordered_bylabel.stride_time(positions_a1, :);
    stride_time_byprompt.a2 = data_ordered_bylabel.stride_time(positions_a2, :);
    stride_time_byprompt.b0 = data_ordered_bylabel.stride_time(positions_b0, :);
    stride_time_byprompt.b1 = data_ordered_bylabel.stride_time(positions_b1, :);
    stride_time_byprompt.b2 = data_ordered_bylabel.stride_time(positions_b2, :);

    stride_time_var_byprompt.s0 = data_ordered_bylabel.stride_time_var(positions_s0, :);
    stride_time_var_byprompt.s1 = data_ordered_bylabel.stride_time_var(positions_s1, :);
    stride_time_var_byprompt.s2 = data_ordered_bylabel.stride_time_var(positions_s2, :);
    stride_time_var_byprompt.a0 = data_ordered_bylabel.stride_time_var(positions_a0, :);
    stride_time_var_byprompt.a1 = data_ordered_bylabel.stride_time_var(positions_a1, :);
    stride_time_var_byprompt.a2 = data_ordered_bylabel.stride_time_var(positions_a2, :);
    stride_time_var_byprompt.b0 = data_ordered_bylabel.stride_time_var(positions_b0, :);
    stride_time_var_byprompt.b1 = data_ordered_bylabel.stride_time_var(positions_b1, :);
    stride_time_var_byprompt.b2 = data_ordered_bylabel.stride_time_var(positions_b2, :);

    step_error_byprompt.s0 = data_ordered_bylabel.step_error(positions_s0, :);
    step_error_byprompt.s1 = data_ordered_bylabel.step_error(positions_s1, :);
    step_error_byprompt.s2 = data_ordered_bylabel.step_error(positions_s2, :);
    step_error_byprompt.a0 = data_ordered_bylabel.step_error(positions_a0, :);
    step_error_byprompt.a1 = data_ordered_bylabel.step_error(positions_a1, :);
    step_error_byprompt.a2 = data_ordered_bylabel.step_error(positions_a2, :);
    step_error_byprompt.b0 = data_ordered_bylabel.step_error(positions_b0, :);
    step_error_byprompt.b1 = data_ordered_bylabel.step_error(positions_b1, :);
    step_error_byprompt.b2 = data_ordered_bylabel.step_error(positions_b2, :);

    step_error_var_byprompt.s0 = data_ordered_bylabel.step_error_std(positions_s0, :);
    step_error_var_byprompt.s1 = data_ordered_bylabel.step_error_std(positions_s1, :);
    step_error_var_byprompt.s2 = data_ordered_bylabel.step_error_std(positions_s2, :);
    step_error_var_byprompt.a0 = data_ordered_bylabel.step_error_std(positions_a0, :);
    step_error_var_byprompt.a1 = data_ordered_bylabel.step_error_std(positions_a1, :);
    step_error_var_byprompt.a2 = data_ordered_bylabel.step_error_std(positions_a2, :);
    step_error_var_byprompt.b0 = data_ordered_bylabel.step_error_std(positions_b0, :);
    step_error_var_byprompt.b1 = data_ordered_bylabel.step_error_std(positions_b1, :);
    step_error_var_byprompt.b2 = data_ordered_bylabel.step_error_std(positions_b2, :);

    head_angle_byprompt.s0 = data_ordered_bylabel.head_angle(positions_s0, :);
    head_angle_byprompt.s1 = data_ordered_bylabel.head_angle(positions_s1, :);
    head_angle_byprompt.s2 = data_ordered_bylabel.head_angle(positions_s2, :);
    head_angle_byprompt.a0 = data_ordered_bylabel.head_angle(positions_a0, :);
    head_angle_byprompt.a1 = data_ordered_bylabel.head_angle(positions_a1, :);
    head_angle_byprompt.a2 = data_ordered_bylabel.head_angle(positions_a2, :);
    head_angle_byprompt.b0 = data_ordered_bylabel.head_angle(positions_b0, :);
    head_angle_byprompt.b1 = data_ordered_bylabel.head_angle(positions_b1, :);
    head_angle_byprompt.b2 = data_ordered_bylabel.head_angle(positions_b2, :);

    trunk_angle_byprompt.s0 = data_ordered_bylabel.trunk_angle(positions_s0, :);
    trunk_angle_byprompt.s1 = data_ordered_bylabel.trunk_angle(positions_s1, :);
    trunk_angle_byprompt.s2 = data_ordered_bylabel.trunk_angle(positions_s2, :);
    trunk_angle_byprompt.a0 = data_ordered_bylabel.trunk_angle(positions_a0, :);
    trunk_angle_byprompt.a1 = data_ordered_bylabel.trunk_angle(positions_a1, :);
    trunk_angle_byprompt.a2 = data_ordered_bylabel.trunk_angle(positions_a2, :);
    trunk_angle_byprompt.b0 = data_ordered_bylabel.trunk_angle(positions_b0, :);
    trunk_angle_byprompt.b1 = data_ordered_bylabel.trunk_angle(positions_b1, :);
    trunk_angle_byprompt.b2 = data_ordered_bylabel.trunk_angle(positions_b2, :);

    dist_gaze_byprompt.s0 = data_ordered_bylabel.dist_gaze(positions_s0, :);
    dist_gaze_byprompt.s1 = data_ordered_bylabel.dist_gaze(positions_s1, :);
    dist_gaze_byprompt.s2 = data_ordered_bylabel.dist_gaze(positions_s2, :);
    dist_gaze_byprompt.a0 = data_ordered_bylabel.dist_gaze(positions_a0, :);
    dist_gaze_byprompt.a1 = data_ordered_bylabel.dist_gaze(positions_a1, :);
    dist_gaze_byprompt.a2 = data_ordered_bylabel.dist_gaze(positions_a2, :);
    dist_gaze_byprompt.b0 = data_ordered_bylabel.dist_gaze(positions_b0, :);
    dist_gaze_byprompt.b1 = data_ordered_bylabel.dist_gaze(positions_b1, :);
    dist_gaze_byprompt.b2 = data_ordered_bylabel.dist_gaze(positions_b2, :);

    surveys_balance_byprompt.s0 = data_ordered_bylabel.surveys_balance(positions_s0, :);
    surveys_balance_byprompt.s1 = data_ordered_bylabel.surveys_balance(positions_s1, :);
    surveys_balance_byprompt.s2 = data_ordered_bylabel.surveys_balance(positions_s2, :);
    surveys_balance_byprompt.a0 = data_ordered_bylabel.surveys_balance(positions_a0, :);
    surveys_balance_byprompt.a1 = data_ordered_bylabel.surveys_balance(positions_a1, :);
    surveys_balance_byprompt.a2 = data_ordered_bylabel.surveys_balance(positions_a2, :);
    surveys_balance_byprompt.b0 = data_ordered_bylabel.surveys_balance(positions_b0, :);
    surveys_balance_byprompt.b1 = data_ordered_bylabel.surveys_balance(positions_b1, :);
    surveys_balance_byprompt.b2 = data_ordered_bylabel.surveys_balance(positions_b2, :);

    surveys_foot_placement_byprompt.s0 = data_ordered_bylabel.surveys_foot_placement(positions_s0, :);
    surveys_foot_placement_byprompt.s1 = data_ordered_bylabel.surveys_foot_placement(positions_s1, :);
    surveys_foot_placement_byprompt.s2 = data_ordered_bylabel.surveys_foot_placement(positions_s2, :);
    surveys_foot_placement_byprompt.a0 = data_ordered_bylabel.surveys_foot_placement(positions_a0, :);
    surveys_foot_placement_byprompt.a1 = data_ordered_bylabel.surveys_foot_placement(positions_a1, :);
    surveys_foot_placement_byprompt.a2 = data_ordered_bylabel.surveys_foot_placement(positions_a2, :);
    surveys_foot_placement_byprompt.b0 = data_ordered_bylabel.surveys_foot_placement(positions_b0, :);
    surveys_foot_placement_byprompt.b1 = data_ordered_bylabel.surveys_foot_placement(positions_b1, :);
    surveys_foot_placement_byprompt.b2 = data_ordered_bylabel.surveys_foot_placement(positions_b2, :);

    surveys_walking_speed_byprompt.s0 = data_ordered_bylabel.surveys_walking_speed(positions_s0, :);
    surveys_walking_speed_byprompt.s1 = data_ordered_bylabel.surveys_walking_speed(positions_s1, :);
    surveys_walking_speed_byprompt.s2 = data_ordered_bylabel.surveys_walking_speed(positions_s2, :);
    surveys_walking_speed_byprompt.a0 = data_ordered_bylabel.surveys_walking_speed(positions_a0, :);
    surveys_walking_speed_byprompt.a1 = data_ordered_bylabel.surveys_walking_speed(positions_a1, :);
    surveys_walking_speed_byprompt.a2 = data_ordered_bylabel.surveys_walking_speed(positions_a2, :);
    surveys_walking_speed_byprompt.b0 = data_ordered_bylabel.surveys_walking_speed(positions_b0, :);
    surveys_walking_speed_byprompt.b1 = data_ordered_bylabel.surveys_walking_speed(positions_b1, :);
    surveys_walking_speed_byprompt.b2 = data_ordered_bylabel.surveys_walking_speed(positions_b2, :);

    surveys_metabolics_byprompt.s0 = data_ordered_bylabel.surveys_metabolics(positions_s0, :);
    surveys_metabolics_byprompt.s1 = data_ordered_bylabel.surveys_metabolics(positions_s1, :);
    surveys_metabolics_byprompt.s2 = data_ordered_bylabel.surveys_metabolics(positions_s2, :);
    surveys_metabolics_byprompt.a0 = data_ordered_bylabel.surveys_metabolics(positions_a0, :);
    surveys_metabolics_byprompt.a1 = data_ordered_bylabel.surveys_metabolics(positions_a1, :);
    surveys_metabolics_byprompt.a2 = data_ordered_bylabel.surveys_metabolics(positions_a2, :);
    surveys_metabolics_byprompt.b0 = data_ordered_bylabel.surveys_metabolics(positions_b0, :);
    surveys_metabolics_byprompt.b1 = data_ordered_bylabel.surveys_metabolics(positions_b1, :);
    surveys_metabolics_byprompt.b2 = data_ordered_bylabel.surveys_metabolics(positions_b2, :);

    surveys_difficulty_byprompt.s0 = data_ordered_bylabel.surveys_difficulty(positions_s0, :);
    surveys_difficulty_byprompt.s1 = data_ordered_bylabel.surveys_difficulty(positions_s1, :);
    surveys_difficulty_byprompt.s2 = data_ordered_bylabel.surveys_difficulty(positions_s2, :);
    surveys_difficulty_byprompt.a0 = data_ordered_bylabel.surveys_difficulty(positions_a0, :);
    surveys_difficulty_byprompt.a1 = data_ordered_bylabel.surveys_difficulty(positions_a1, :);
    surveys_difficulty_byprompt.a2 = data_ordered_bylabel.surveys_difficulty(positions_a2, :);
    surveys_difficulty_byprompt.b0 = data_ordered_bylabel.surveys_difficulty(positions_b0, :);
    surveys_difficulty_byprompt.b1 = data_ordered_bylabel.surveys_difficulty(positions_b1, :);
    surveys_difficulty_byprompt.b2 = data_ordered_bylabel.surveys_difficulty(positions_b2, :);

   % Calculate the correlations
    correlation_walkingspeed_walkingspeed_std = corr(walking_speed_vector, walking_speed_std_vector);
    correlation_walkingspeed_metabolics = corr(walking_speed_vector, metabolics_vector);
    correlation_walkingspeed_step_width = corr(walking_speed_vector, step_width_vector);
    correlation_walkingspeed_step_width_var = corr(walking_speed_vector, step_width_var_vector);
    correlation_walkingspeed_step_length = corr(walking_speed_vector, step_length_vector);
    correlation_walkingspeed_step_length_var = corr(walking_speed_vector, step_length_var_vector);
    correlation_walkingspeed_stride_time = corr(walking_speed_vector, stride_time_vector);
    correlation_walkingspeed_stride_time_var = corr(walking_speed_vector, stride_time_var_vector);
    correlation_walkingspeed_step_error = corr(walking_speed_vector, step_error_vector);
    correlation_walkingspeed_step_error_var = corr(walking_speed_vector, step_error_var_vector);
    correlation_walkingspeed_head_angle = corr(walking_speed_vector, head_angle_vector);
    correlation_walkingspeed_trunk_angle = corr(walking_speed_vector, trunk_angle_vector);
    correlation_walkingspeed_dist_gaze = corr(walking_speed_vector, dist_gaze_vector);
    correlation_walkingspeed_surveys_balance = corr(walking_speed_vector, surveys_balance_vector);
    correlation_walkingspeed_surveys_foot_placement = corr(walking_speed_vector, surveys_foot_placement_vector);
    correlation_walkingspeed_surveys_walking_speed = corr(walking_speed_vector, surveys_walking_speed_vector);
    correlation_walkingspeed_surveys_metabolics = corr(walking_speed_vector, surveys_metabolics_vector);
    correlation_walkingspeed_surveys_difficulty = corr(walking_speed_vector, surveys_difficulty_vector);
    correlation_walkingspeed_pspeed = Correlation_prompts(walking_speed_byprompt, "s0", "s1", "s2");
    correlation_walkingspeed_paccuracy = Correlation_prompts(walking_speed_byprompt, "a0", "a1", "a2");
    correlation_walkingspeed_pbalance = Correlation_prompts(walking_speed_byprompt, "b0", "b1", "b2");

    correlation_walkingspeed_std_metabolics = corr(walking_speed_std_vector, metabolics_vector);
    correlation_walkingspeed_std_step_width = corr(walking_speed_std_vector, step_width_vector);
    correlation_walkingspeed_std_step_width_var = corr(walking_speed_std_vector, step_width_var_vector);
    correlation_walkingspeed_std_step_length = corr(walking_speed_std_vector, step_length_vector);
    correlation_walkingspeed_std_step_length_var = corr(walking_speed_std_vector, step_length_var_vector);
    correlation_walkingspeed_std_stride_time = corr(walking_speed_std_vector, stride_time_vector);
    correlation_walkingspeed_std_stride_time_var = corr(walking_speed_std_vector, stride_time_var_vector);
    correlation_walkingspeed_std_step_error = corr(walking_speed_std_vector, step_error_vector);
    correlation_walkingspeed_std_step_error_var = corr(walking_speed_std_vector, step_error_var_vector);
    correlation_walkingspeed_std_head_angle = corr(walking_speed_std_vector, head_angle_vector);
    correlation_walkingspeed_std_trunk_angle = corr(walking_speed_std_vector, trunk_angle_vector);
    correlation_walkingspeed_std_dist_gaze = corr(walking_speed_std_vector, dist_gaze_vector);
    correlation_walkingspeed_std_surveys_balance = corr(walking_speed_std_vector, surveys_balance_vector);
    correlation_walkingspeed_std_surveys_foot_placement = corr(walking_speed_std_vector, surveys_foot_placement_vector);
    correlation_walkingspeed_std_surveys_walking_speed = corr(walking_speed_std_vector, surveys_walking_speed_vector);
    correlation_walkingspeed_std_surveys_metabolics = corr(walking_speed_std_vector, surveys_metabolics_vector);
    correlation_walkingspeed_std_surveys_difficulty = corr(walking_speed_std_vector, surveys_difficulty_vector);
    correlation_walkingspeed_std_pspeed = Correlation_prompts(walking_speed_std_byprompt, "s0", "s1", "s2");
    correlation_walkingspeed_std_paccuracy = Correlation_prompts(walking_speed_std_byprompt, "a0", "a1", "a2");
    correlation_walkingspeed_std_pbalance = Correlation_prompts(walking_speed_std_byprompt, "b0", "b1", "b2");

    correlation_metabolics_step_width = corr(metabolics_vector,step_width_vector);
    correlation_metabolics_step_width_var = corr(metabolics_vector, step_width_var_vector);
    correlation_metabolics_step_length = corr(metabolics_vector, step_length_vector);
    correlation_metabolics_step_length_var = corr(metabolics_vector, step_length_var_vector);
    correlation_metabolics_stride_time = corr(metabolics_vector, stride_time_vector);
    correlation_metabolics_stride_time_var = corr(metabolics_vector, stride_time_var_vector);
    correlation_metabolics_step_error = corr(metabolics_vector, step_error_vector);
    correlation_metabolics_step_error_var = corr(metabolics_vector, step_error_var_vector);
    correlation_metabolics_head_angle = corr(metabolics_vector, head_angle_vector);
    correlation_metabolics_trunk_angle = corr(metabolics_vector, trunk_angle_vector);
    correlation_metabolics_dist_gaze = corr(metabolics_vector, dist_gaze_vector);
    correlation_metabolics_surveys_balance = corr(metabolics_vector, surveys_balance_vector);
    correlation_metabolics_surveys_foot_placement = corr(metabolics_vector, surveys_foot_placement_vector);
    correlation_metabolics_surveys_walking_speed = corr(metabolics_vector, surveys_walking_speed_vector);
    correlation_metabolics_surveys_metabolics = corr(metabolics_vector, surveys_metabolics_vector);
    correlation_metabolics_surveys_difficulty = corr(metabolics_vector, surveys_difficulty_vector);
    correlation_metabolics_pspeed = Correlation_prompts(metabolics_byprompt, "s0", "s1", "s2");
    correlation_metabolics_paccuracy = Correlation_prompts(metabolics_byprompt, "a0", "a1", "a2");
    correlation_metabolics_pbalance = Correlation_prompts(metabolics_byprompt, "b0", "b1", "b2");

    correlation_step_width_step_width_var = corr(step_width_vector, step_width_var_vector);
    correlation_step_width_step_length = corr(step_width_vector, step_length_vector);
    correlation_step_width_step_length_var = corr(step_width_vector, step_length_var_vector);
    correlation_step_width_stride_time = corr(step_width_vector, stride_time_vector);
    correlation_step_width_stride_time_var = corr(step_width_vector, stride_time_var_vector);
    correlation_step_width_step_error = corr(step_width_vector, step_error_vector);
    correlation_step_width_step_error_var = corr(step_width_vector, step_error_var_vector);
    correlation_step_width_head_angle = corr(step_width_vector, head_angle_vector);
    correlation_step_width_trunk_angle = corr(step_width_vector, trunk_angle_vector);
    correlation_step_width_dist_gaze = corr(step_width_vector, dist_gaze_vector);
    correlation_step_width_surveys_balance = corr(step_width_vector, surveys_balance_vector);
    correlation_step_width_surveys_foot_placement = corr(step_width_vector, surveys_foot_placement_vector);
    correlation_step_width_surveys_walking_speed = corr(step_width_vector, surveys_walking_speed_vector);
    correlation_step_width_surveys_metabolics = corr(step_width_vector, surveys_metabolics_vector);
    correlation_step_width_surveys_difficulty = corr(step_width_vector, surveys_difficulty_vector);
    correlation_step_width_pspeed = Correlation_prompts(step_width_byprompt, "s0", "s1", "s2");
    correlation_step_width_paccuracy = Correlation_prompts(step_width_byprompt, "a0", "a1", "a2");
    correlation_step_width_pbalance = Correlation_prompts(step_width_byprompt, "b0", "b1", "b2");

    correlation_step_width_var_step_length = corr(step_width_var_vector, step_length_vector);
    correlation_step_width_var_step_length_var = corr(step_width_var_vector, step_length_var_vector);
    correlation_step_width_var_stride_time = corr(step_width_var_vector, stride_time_vector);
    correlation_step_width_var_stride_time_var = corr(step_width_var_vector, stride_time_var_vector);
    correlation_step_width_var_step_error = corr(step_width_var_vector, step_error_vector);
    correlation_step_width_var_step_error_var = corr(step_width_var_vector, step_error_var_vector);
    correlation_step_width_var_head_angle = corr(step_width_var_vector, head_angle_vector);
    correlation_step_width_var_trunk_angle = corr(step_width_var_vector, trunk_angle_vector);
    correlation_step_width_var_dist_gaze = corr(step_width_var_vector, dist_gaze_vector);
    correlation_step_width_var_surveys_balance = corr(step_width_var_vector, surveys_balance_vector);
    correlation_step_width_var_surveys_foot_placement = corr(step_width_var_vector, surveys_foot_placement_vector);
    correlation_step_width_var_surveys_walking_speed = corr(step_width_var_vector, surveys_walking_speed_vector);
    correlation_step_width_var_surveys_metabolics = corr(step_width_var_vector, surveys_metabolics_vector);
    correlation_step_width_var_surveys_difficulty = corr(step_width_var_vector, surveys_difficulty_vector);
    correlation_step_width_var_pspeed = Correlation_prompts(step_width_var_byprompt, "s0", "s1", "s2");
    correlation_step_width_var_paccuracy = Correlation_prompts(step_width_var_byprompt, "a0", "a1", "a2");
    correlation_step_width_var_pbalance = Correlation_prompts(step_width_var_byprompt, "b0", "b1", "b2");

    correlation_step_length_step_length_var = corr(step_length_vector, step_length_var_vector);
    correlation_step_length_stride_time = corr(step_length_vector, stride_time_vector);
    correlation_step_length_stride_time_var = corr(step_length_vector, stride_time_var_vector);
    correlation_step_length_step_error = corr(step_length_vector, step_error_vector);
    correlation_step_length_step_error_var = corr(step_length_vector, step_error_var_vector);
    correlation_step_length_head_angle = corr(step_length_vector, head_angle_vector);
    correlation_step_length_trunk_angle = corr(step_length_vector, trunk_angle_vector);
    correlation_step_length_dist_gaze = corr(step_length_vector, dist_gaze_vector);
    correlation_step_length_surveys_balance = corr(step_length_vector, surveys_balance_vector);
    correlation_step_length_surveys_foot_placement = corr(step_length_vector, surveys_foot_placement_vector);
    correlation_step_length_surveys_walking_speed = corr(step_length_vector, surveys_walking_speed_vector);
    correlation_step_length_surveys_metabolics = corr(step_length_vector, surveys_metabolics_vector);
    correlation_step_length_surveys_difficulty = corr(step_length_vector, surveys_difficulty_vector);
    correlation_step_length_pspeed = Correlation_prompts(step_length_byprompt, "s0", "s1", "s2");
    correlation_step_length_paccuracy = Correlation_prompts(step_length_byprompt, "a0", "a1", "a2");
    correlation_step_length_pbalance = Correlation_prompts(step_length_byprompt, "b0", "b1", "b2");

    correlation_step_length_var_stride_time = corr(step_length_var_vector, stride_time_vector);
    correlation_step_length_var_stride_time_var = corr(step_length_var_vector, stride_time_var_vector);
    correlation_step_length_var_step_error = corr(step_length_var_vector, step_error_vector);
    correlation_step_length_var_step_error_var = corr(step_length_var_vector, step_error_var_vector);
    correlation_step_length_var_head_angle = corr(step_length_var_vector, head_angle_vector);
    correlation_step_length_var_trunk_angle = corr(step_length_var_vector, trunk_angle_vector);
    correlation_step_length_var_dist_gaze = corr(step_length_var_vector, dist_gaze_vector);
    correlation_step_length_var_surveys_balance = corr(step_length_var_vector, surveys_balance_vector);
    correlation_step_length_var_surveys_foot_placement = corr(step_length_var_vector, surveys_foot_placement_vector);
    correlation_step_length_var_surveys_walking_speed = corr(step_length_var_vector, surveys_walking_speed_vector);
    correlation_step_length_var_surveys_metabolics = corr(step_length_var_vector, surveys_metabolics_vector);
    correlation_step_length_var_surveys_difficulty = corr(step_length_var_vector, surveys_difficulty_vector);
    correlation_step_length_var_pspeed = Correlation_prompts(step_length_var_byprompt, "s0", "s1", "s2");
    correlation_step_length_var_paccuracy = Correlation_prompts(step_length_var_byprompt, "a0", "a1", "a2");
    correlation_step_length_var_pbalance = Correlation_prompts(step_length_var_byprompt, "b0", "b1", "b2");

    correlation_stride_time_stride_time_var = corr(stride_time_vector, stride_time_var_vector);
    correlation_stride_time_step_error = corr(stride_time_vector, step_error_vector);
    correlation_stride_time_step_error_var = corr(stride_time_vector, step_error_var_vector);
    correlation_stride_time_head_angle = corr(stride_time_vector, head_angle_vector);
    correlation_stride_time_trunk_angle = corr(stride_time_vector, trunk_angle_vector);
    correlation_stride_time_dist_gaze = corr(stride_time_vector, dist_gaze_vector);
    correlation_stride_time_surveys_balance = corr(stride_time_vector, surveys_balance_vector);
    correlation_stride_time_surveys_foot_placement = corr(stride_time_vector, surveys_foot_placement_vector);
    correlation_stride_time_surveys_walking_speed = corr(stride_time_vector, surveys_walking_speed_vector);
    correlation_stride_time_surveys_metabolics = corr(stride_time_vector, surveys_metabolics_vector);
    correlation_stride_time_surveys_difficulty = corr(stride_time_vector, surveys_difficulty_vector);
    correlation_stride_time_pspeed = Correlation_prompts(stride_time_byprompt, "s0", "s1", "s2");
    correlation_stride_time_paccuracy = Correlation_prompts(stride_time_byprompt, "a0", "a1", "a2");
    correlation_stride_time_pbalance = Correlation_prompts(stride_time_byprompt, "b0", "b1", "b2");

    correlation_stride_time_var_step_error = corr(stride_time_var_vector, step_error_vector);
    correlation_stride_time_var_step_error_var = corr(stride_time_var_vector, step_error_var_vector);
    correlation_stride_time_var_head_angle = corr(stride_time_var_vector, head_angle_vector);
    correlation_stride_time_var_trunk_angle = corr(stride_time_var_vector, trunk_angle_vector);
    correlation_stride_time_var_dist_gaze = corr(stride_time_var_vector, dist_gaze_vector);
    correlation_stride_time_var_surveys_balance = corr(stride_time_var_vector, surveys_balance_vector);
    correlation_stride_time_var_surveys_foot_placement = corr(stride_time_var_vector, surveys_foot_placement_vector);
    correlation_stride_time_var_surveys_walking_speed = corr(stride_time_var_vector, surveys_walking_speed_vector);
    correlation_stride_time_var_surveys_metabolics = corr(stride_time_var_vector, surveys_metabolics_vector);
    correlation_stride_time_var_surveys_difficulty = corr(stride_time_var_vector, surveys_difficulty_vector);
    correlation_stride_time_var_pspeed = Correlation_prompts(stride_time_var_byprompt, "s0", "s1", "s2");
    correlation_stride_time_var_paccuracy = Correlation_prompts(stride_time_var_byprompt, "a0", "a1", "a2");
    correlation_stride_time_var_pbalance = Correlation_prompts(stride_time_var_byprompt, "b0", "b1", "b2");

    correlation_step_error_step_error_var = corr(step_error_vector, step_error_var_vector);
    correlation_step_error_head_angle = corr(step_error_vector, head_angle_vector);
    correlation_step_error_trunk_angle = corr(step_error_vector, trunk_angle_vector);
    correlation_step_error_dist_gaze = corr(step_error_vector, dist_gaze_vector);
    correlation_step_error_surveys_balance = corr(step_error_vector, surveys_balance_vector);
    correlation_step_error_surveys_foot_placement = corr(step_error_vector, surveys_foot_placement_vector);
    correlation_step_error_surveys_walking_speed = corr(step_error_vector, surveys_walking_speed_vector);
    correlation_step_error_surveys_metabolics = corr(step_error_vector, surveys_metabolics_vector);
    correlation_step_error_surveys_difficulty = corr(step_error_vector, surveys_difficulty_vector);
    correlation_step_error_pspeed = Correlation_prompts(step_error_byprompt, "s0", "s1", "s2");
    correlation_step_error_paccuracy = Correlation_prompts(step_error_byprompt, "a0", "a1", "a2");
    correlation_step_error_pbalance = Correlation_prompts(step_error_byprompt, "b0", "b1", "b2");

    correlation_step_error_var_head_angle = corr(step_error_var_vector, head_angle_vector);
    correlation_step_error_var_trunk_angle = corr(step_error_var_vector, trunk_angle_vector);
    correlation_step_error_var_dist_gaze = corr(step_error_var_vector, dist_gaze_vector);
    correlation_step_error_var_surveys_balance = corr(step_error_var_vector, surveys_balance_vector);
    correlation_step_error_var_surveys_foot_placement = corr(step_error_var_vector, surveys_foot_placement_vector);
    correlation_step_error_var_surveys_walking_speed = corr(step_error_var_vector, surveys_walking_speed_vector);
    correlation_step_error_var_surveys_metabolics = corr(step_error_var_vector, surveys_metabolics_vector);
    correlation_step_error_var_surveys_difficulty = corr(step_error_var_vector, surveys_difficulty_vector);
    correlation_step_error_var_pspeed = Correlation_prompts(step_error_var_byprompt, "s0", "s1", "s2");
    correlation_step_error_var_paccuracy = Correlation_prompts(step_error_var_byprompt, "a0", "a1", "a2");
    correlation_step_error_var_pbalance = Correlation_prompts(step_error_var_byprompt, "b0", "b1", "b2");

    correlation_head_angle_trunk_angle = corr(head_angle_vector, trunk_angle_vector);
    correlation_head_angle_dist_gaze = corr(head_angle_vector, dist_gaze_vector);
    correlation_head_angle_surveys_balance = corr(head_angle_vector, surveys_balance_vector);
    correlation_head_angle_surveys_foot_placement = corr(head_angle_vector, surveys_foot_placement_vector);
    correlation_head_angle_surveys_walking_speed = corr(head_angle_vector, surveys_walking_speed_vector);
    correlation_head_angle_surveys_metabolics = corr(head_angle_vector, surveys_metabolics_vector);
    correlation_head_angle_surveys_difficulty = corr(head_angle_vector, surveys_difficulty_vector);
    correlation_head_angle_pspeed = Correlation_prompts(head_angle_byprompt, "s0", "s1", "s2");
    correlation_head_angle_paccuracy = Correlation_prompts(head_angle_byprompt, "a0", "a1", "a2");
    correlation_head_angle_pbalance = Correlation_prompts(head_angle_byprompt, "b0", "b1", "b2");

    correlation_trunk_angle_dist_gaze = corr(trunk_angle_vector, dist_gaze_vector);
    correlation_trunk_angle_surveys_balance = corr(trunk_angle_vector, surveys_balance_vector);
    correlation_trunk_angle_surveys_foot_placement = corr(trunk_angle_vector, surveys_foot_placement_vector);
    correlation_trunk_angle_surveys_walking_speed = corr(trunk_angle_vector, surveys_walking_speed_vector);
    correlation_trunk_angle_surveys_metabolics = corr(trunk_angle_vector, surveys_metabolics_vector);
    correlation_trunk_angle_surveys_difficulty = corr(trunk_angle_vector, surveys_difficulty_vector);
    correlation_trunk_angle_pspeed = Correlation_prompts(trunk_angle_byprompt, "s0", "s1", "s2");
    correlation_trunk_angle_paccuracy = Correlation_prompts(trunk_angle_byprompt, "a0", "a1", "a2");
    correlation_trunk_angle_pbalance = Correlation_prompts(trunk_angle_byprompt, "b0", "b1", "b2");

    correlation_dist_gaze_surveys_balance = corr(dist_gaze_vector, surveys_balance_vector);
    correlation_dist_gaze_surveys_foot_placement = corr(dist_gaze_vector, surveys_foot_placement_vector);
    correlation_dist_gaze_surveys_walking_speed = corr(dist_gaze_vector, surveys_walking_speed_vector);
    correlation_dist_gaze_surveys_metabolics = corr(dist_gaze_vector, surveys_metabolics_vector);
    correlation_dist_gaze_surveys_difficulty = corr(dist_gaze_vector, surveys_difficulty_vector);
    correlation_dist_gaze_pspeed = Correlation_prompts(dist_gaze_byprompt, "s0", "s1", "s2");
    correlation_dist_gaze_paccuracy = Correlation_prompts(dist_gaze_byprompt, "a0", "a1", "a2");
    correlation_dist_gaze_pbalance = Correlation_prompts(dist_gaze_byprompt, "b0", "b1", "b2");

    correlation_surveys_balance_surveys_foot_placement = corr(surveys_balance_vector, surveys_foot_placement_vector);
    correlation_surveys_balance_surveys_walking_speed = corr(surveys_balance_vector, surveys_walking_speed_vector);
    correlation_surveys_balance_surveys_metabolics = corr(surveys_balance_vector, surveys_metabolics_vector);
    correlation_surveys_balance_surveys_difficulty = corr(surveys_balance_vector, surveys_difficulty_vector);
    correlation_surveys_balance_pspeed = Correlation_prompts(surveys_balance_byprompt, "s0", "s1", "s2");
    correlation_surveys_balance_paccuracy = Correlation_prompts(surveys_balance_byprompt, "a0", "a1", "a2");
    correlation_surveys_balance_pbalance = Correlation_prompts(surveys_balance_byprompt, "b0", "b1", "b2");

    correlation_surveys_foot_placement_surveys_walking_speed = corr(surveys_foot_placement_vector, surveys_walking_speed_vector);
    correlation_surveys_foot_placement_surveys_metabolics = corr(surveys_foot_placement_vector, surveys_metabolics_vector);
    correlation_surveys_foot_placement_surveys_difficulty = corr(surveys_foot_placement_vector, surveys_difficulty_vector);
    correlation_surveys_foot_placement_pspeed = Correlation_prompts(surveys_foot_placement_byprompt, "s0", "s1", "s2");
    correlation_surveys_foot_placement_paccuracy = Correlation_prompts(surveys_foot_placement_byprompt, "a0", "a1", "a2");
    correlation_surveys_foot_placement_pbalance = Correlation_prompts(surveys_foot_placement_byprompt, "b0", "b1", "b2");

    correlation_surveys_walking_speed_surveys_metabolics = corr(surveys_walking_speed_vector, surveys_metabolics_vector);
    correlation_surveys_walking_speed_surveys_difficulty = corr(surveys_walking_speed_vector, surveys_difficulty_vector);
    correlation_surveys_walking_speed_pspeed = Correlation_prompts(surveys_walking_speed_byprompt, "s0", "s1", "s2");
    correlation_surveys_walking_speed_paccuracy = Correlation_prompts(surveys_walking_speed_byprompt, "a0", "a1", "a2");
    correlation_surveys_walking_speed_pbalance = Correlation_prompts(surveys_walking_speed_byprompt, "b0", "b1", "b2");

    correlation_surveys_metabolics_surveys_difficulty = corr(surveys_metabolics_vector, surveys_difficulty_vector);
    correlation_surveys_metabolics_pspeed = Correlation_prompts(surveys_metabolics_byprompt, "s0", "s1", "s2");
    correlation_surveys_metabolics_paccuracy = Correlation_prompts(surveys_metabolics_byprompt, "a0", "a1", "a2");
    correlation_surveys_metabolics_pbalance = Correlation_prompts(surveys_metabolics_byprompt, "b0", "b1", "b2");

    correlation_surveys_difficulty_pspeed = Correlation_prompts(surveys_difficulty_byprompt, "s0", "s1", "s2");
    correlation_surveys_difficulty_paccuracy = Correlation_prompts(surveys_difficulty_byprompt, "a0", "a1", "a2");
    correlation_surveys_difficulty_pbalance = Correlation_prompts(surveys_difficulty_byprompt, "b0", "b1", "b2");

    % Plot
    variables = {'Walking Speed', 'Walking Speed Var', 'Energy Expenditure', 'Step Width', 'Step Width Var', 'Step Length', 'Step Length Var', ...
        'Stride Duration', 'Stride Duration Var', 'Step Error', 'Step Error Var', 'Head Angle', 'Trunk Angle', 'Look-Ahead Distance', 'Balance Priority', ...
        'Foot Placement Priority', 'Walking Speed Priority', 'Energy Expenditure Priority', 'Self-Perceived Difficulty', '\it Visual Disturbance\it', '\it Speed Prompt\it', ...
        '\it Accuracy Prompt\it'};

    variables_x = variables(1:end-3); % Remove 'Speed Prompt', 'Accuracy Prompt', 'Balance Prompt' from the list
    
    % Without dist gaze
    % correlation_matrix = [
    %     1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN; 
    %     correlation_walkingspeed_walkingspeed_std, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
    %     correlation_walkingspeed_metabolics, correlation_walkingspeed_std_metabolics, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
    %     correlation_walkingspeed_step_width, correlation_walkingspeed_std_step_width, correlation_metabolics_step_width, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
    %     correlation_walkingspeed_step_width_var, correlation_walkingspeed_std_step_width_var, correlation_metabolics_step_width_var, correlation_step_width_step_width_var, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
    %     correlation_walkingspeed_step_length, correlation_walkingspeed_std_step_length, correlation_metabolics_step_length, correlation_step_width_step_length, correlation_step_width_var_step_length, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
    %     correlation_walkingspeed_step_length_var, correlation_walkingspeed_std_step_length_var, correlation_metabolics_step_length_var, correlation_step_width_step_length_var, correlation_step_width_var_step_length_var, correlation_step_length_step_length_var, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
    %     correlation_walkingspeed_stride_time, correlation_walkingspeed_std_stride_time, correlation_metabolics_stride_time, correlation_step_width_stride_time, correlation_step_width_var_stride_time, correlation_step_length_stride_time, correlation_step_length_var_stride_time, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
    %     correlation_walkingspeed_stride_time_var, correlation_walkingspeed_std_stride_time_var, correlation_metabolics_stride_time_var, correlation_step_width_stride_time_var, correlation_step_width_var_stride_time_var, correlation_step_length_stride_time_var, correlation_step_length_var_stride_time_var, correlation_stride_time_stride_time_var, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
    %     correlation_walkingspeed_step_error, correlation_walkingspeed_std_step_error, correlation_metabolics_step_error, correlation_step_width_step_error, correlation_step_width_var_step_error, correlation_step_length_step_error, correlation_step_length_var_step_error, correlation_stride_time_step_error, correlation_stride_time_var_step_error, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
    %     correlation_walkingspeed_step_error_var, correlation_walkingspeed_std_step_error_var, correlation_metabolics_step_error_var, correlation_step_width_step_error_var, correlation_step_width_var_step_error_var, correlation_step_length_step_error_var, correlation_step_length_var_step_error_var, correlation_stride_time_step_error_var, correlation_stride_time_var_step_error_var, correlation_step_error_step_error_var, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
    %     correlation_walkingspeed_head_angle, correlation_walkingspeed_std_head_angle, correlation_metabolics_head_angle, correlation_step_width_head_angle, correlation_step_width_var_head_angle, correlation_step_length_head_angle, correlation_step_length_var_head_angle, correlation_stride_time_head_angle, correlation_stride_time_var_head_angle, correlation_step_error_head_angle, correlation_step_error_var_head_angle, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
    %     correlation_walkingspeed_trunk_angle, correlation_walkingspeed_std_trunk_angle, correlation_metabolics_trunk_angle, correlation_step_width_trunk_angle, correlation_step_width_var_trunk_angle, correlation_step_length_trunk_angle, correlation_step_length_var_trunk_angle, correlation_stride_time_trunk_angle, correlation_stride_time_var_trunk_angle, correlation_step_error_trunk_angle, correlation_step_error_var_trunk_angle, correlation_head_angle_trunk_angle, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
    %     correlation_walkingspeed_surveys_balance, correlation_walkingspeed_std_surveys_balance, correlation_metabolics_surveys_balance, correlation_step_width_surveys_balance, correlation_step_width_var_surveys_balance, correlation_step_length_surveys_balance, correlation_step_length_var_surveys_balance, correlation_stride_time_surveys_balance, correlation_stride_time_var_surveys_balance, correlation_step_error_surveys_balance, correlation_step_error_var_surveys_balance, correlation_head_angle_surveys_balance, correlation_trunk_angle_surveys_balance, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
    %     correlation_walkingspeed_surveys_foot_placement, correlation_walkingspeed_std_surveys_foot_placement, correlation_metabolics_surveys_foot_placement, correlation_step_width_surveys_foot_placement, correlation_step_width_var_surveys_foot_placement, correlation_step_length_surveys_foot_placement, correlation_step_length_var_surveys_foot_placement, correlation_stride_time_surveys_foot_placement, correlation_stride_time_var_surveys_foot_placement, correlation_step_error_surveys_foot_placement, correlation_step_error_var_surveys_foot_placement, correlation_head_angle_surveys_foot_placement, correlation_trunk_angle_surveys_foot_placement, correlation_surveys_balance_surveys_foot_placement, 1, NaN, NaN, NaN, NaN, NaN, NaN;
    %     correlation_walkingspeed_surveys_walking_speed, correlation_walkingspeed_std_surveys_walking_speed, correlation_metabolics_surveys_walking_speed, correlation_step_width_surveys_walking_speed, correlation_step_width_var_surveys_walking_speed, correlation_step_length_surveys_walking_speed, correlation_step_length_var_surveys_walking_speed, correlation_stride_time_surveys_walking_speed, correlation_stride_time_var_surveys_walking_speed, correlation_step_error_surveys_walking_speed, correlation_step_error_var_surveys_walking_speed, correlation_head_angle_surveys_walking_speed, correlation_trunk_angle_surveys_walking_speed, correlation_surveys_balance_surveys_walking_speed, correlation_surveys_foot_placement_surveys_walking_speed, 1, NaN, NaN, NaN, NaN, NaN;
    %     correlation_walkingspeed_surveys_metabolics, correlation_walkingspeed_std_surveys_metabolics, correlation_metabolics_surveys_metabolics, correlation_step_width_surveys_metabolics, correlation_step_width_var_surveys_metabolics, correlation_step_length_surveys_metabolics, correlation_step_length_var_surveys_metabolics, correlation_stride_time_surveys_metabolics, correlation_stride_time_var_surveys_metabolics, correlation_step_error_surveys_metabolics, correlation_step_error_var_surveys_metabolics, correlation_head_angle_surveys_metabolics, correlation_trunk_angle_surveys_metabolics, correlation_surveys_balance_surveys_metabolics, correlation_surveys_foot_placement_surveys_metabolics, correlation_surveys_walking_speed_surveys_metabolics, 1, NaN, NaN, NaN, NaN;
    %     correlation_walkingspeed_surveys_difficulty, correlation_walkingspeed_std_surveys_difficulty, correlation_metabolics_surveys_difficulty, correlation_step_width_surveys_difficulty, correlation_step_width_var_surveys_difficulty, correlation_step_length_surveys_difficulty, correlation_step_length_var_surveys_difficulty, correlation_stride_time_surveys_difficulty, correlation_stride_time_var_surveys_difficulty, correlation_step_error_surveys_difficulty, correlation_step_error_var_surveys_difficulty, correlation_head_angle_surveys_difficulty, correlation_trunk_angle_surveys_difficulty, correlation_surveys_balance_surveys_difficulty, correlation_surveys_foot_placement_surveys_difficulty, correlation_surveys_walking_speed_surveys_difficulty, correlation_surveys_metabolics_surveys_difficulty, 1, NaN, NaN, NaN;        
    %     correlation_walkingspeed_pbalance, correlation_walkingspeed_std_pbalance, correlation_metabolics_pbalance, correlation_step_width_pbalance, correlation_step_width_var_pbalance, correlation_step_length_pbalance, correlation_step_length_var_pbalance, correlation_stride_time_pbalance, correlation_stride_time_var_pbalance, correlation_step_error_pbalance, correlation_step_error_var_pbalance, correlation_head_angle_pbalance, correlation_trunk_angle_pbalance, correlation_surveys_balance_pbalance, correlation_surveys_foot_placement_pbalance, correlation_surveys_walking_speed_pbalance, correlation_surveys_metabolics_pbalance, correlation_surveys_difficulty_pbalance, NaN, NaN, NaN;
    %     correlation_walkingspeed_pspeed, correlation_walkingspeed_std_pspeed, correlation_metabolics_pspeed, correlation_step_width_pspeed, correlation_step_width_var_pspeed, correlation_step_length_pspeed, correlation_step_length_var_pspeed, correlation_stride_time_pspeed, correlation_stride_time_var_pspeed, correlation_step_error_pspeed, correlation_step_error_var_pspeed, correlation_head_angle_pspeed, correlation_trunk_angle_pspeed, correlation_surveys_balance_pspeed, correlation_surveys_foot_placement_pspeed, correlation_surveys_walking_speed_pspeed, correlation_surveys_metabolics_pspeed, correlation_surveys_difficulty_pspeed, NaN, NaN, NaN;
    %     correlation_walkingspeed_paccuracy, correlation_walkingspeed_std_paccuracy, correlation_metabolics_paccuracy, correlation_step_width_paccuracy, correlation_step_width_var_paccuracy, correlation_step_length_paccuracy, correlation_step_length_var_paccuracy, correlation_stride_time_paccuracy, correlation_stride_time_var_paccuracy, correlation_step_error_paccuracy, correlation_step_error_var_paccuracy, correlation_head_angle_paccuracy, correlation_trunk_angle_paccuracy, correlation_surveys_balance_paccuracy, correlation_surveys_foot_placement_paccuracy, correlation_surveys_walking_speed_paccuracy, correlation_surveys_metabolics_paccuracy, correlation_surveys_difficulty_paccuracy, NaN, NaN, NaN;
    % ];

    correlation_matrix = [
        1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN; 
        correlation_walkingspeed_walkingspeed_std, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
        correlation_walkingspeed_metabolics, correlation_walkingspeed_std_metabolics, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
        correlation_walkingspeed_step_width, correlation_walkingspeed_std_step_width, correlation_metabolics_step_width, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
        correlation_walkingspeed_step_width_var, correlation_walkingspeed_std_step_width_var, correlation_metabolics_step_width_var, correlation_step_width_step_width_var, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
        correlation_walkingspeed_step_length, correlation_walkingspeed_std_step_length, correlation_metabolics_step_length, correlation_step_width_step_length, correlation_step_width_var_step_length, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
        correlation_walkingspeed_step_length_var, correlation_walkingspeed_std_step_length_var, correlation_metabolics_step_length_var, correlation_step_width_step_length_var, correlation_step_width_var_step_length_var, correlation_step_length_step_length_var, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
        correlation_walkingspeed_stride_time, correlation_walkingspeed_std_stride_time, correlation_metabolics_stride_time, correlation_step_width_stride_time, correlation_step_width_var_stride_time, correlation_step_length_stride_time, correlation_step_length_var_stride_time, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
        correlation_walkingspeed_stride_time_var, correlation_walkingspeed_std_stride_time_var, correlation_metabolics_stride_time_var, correlation_step_width_stride_time_var, correlation_step_width_var_stride_time_var, correlation_step_length_stride_time_var, correlation_step_length_var_stride_time_var, correlation_stride_time_stride_time_var, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
        correlation_walkingspeed_step_error, correlation_walkingspeed_std_step_error, correlation_metabolics_step_error, correlation_step_width_step_error, correlation_step_width_var_step_error, correlation_step_length_step_error, correlation_step_length_var_step_error, correlation_stride_time_step_error, correlation_stride_time_var_step_error, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
        correlation_walkingspeed_step_error_var, correlation_walkingspeed_std_step_error_var, correlation_metabolics_step_error_var, correlation_step_width_step_error_var, correlation_step_width_var_step_error_var, correlation_step_length_step_error_var, correlation_step_length_var_step_error_var, correlation_stride_time_step_error_var, correlation_stride_time_var_step_error_var, correlation_step_error_step_error_var, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
        correlation_walkingspeed_head_angle, correlation_walkingspeed_std_head_angle, correlation_metabolics_head_angle, correlation_step_width_head_angle, correlation_step_width_var_head_angle, correlation_step_length_head_angle, correlation_step_length_var_head_angle, correlation_stride_time_head_angle, correlation_stride_time_var_head_angle, correlation_step_error_head_angle, correlation_step_error_var_head_angle, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
        correlation_walkingspeed_trunk_angle, correlation_walkingspeed_std_trunk_angle, correlation_metabolics_trunk_angle, correlation_step_width_trunk_angle, correlation_step_width_var_trunk_angle, correlation_step_length_trunk_angle, correlation_step_length_var_trunk_angle, correlation_stride_time_trunk_angle, correlation_stride_time_var_trunk_angle, correlation_step_error_trunk_angle, correlation_step_error_var_trunk_angle, correlation_head_angle_trunk_angle, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
        correlation_walkingspeed_dist_gaze, correlation_walkingspeed_std_dist_gaze, correlation_metabolics_dist_gaze, correlation_step_width_dist_gaze, correlation_step_width_var_dist_gaze, correlation_step_length_dist_gaze, correlation_step_length_var_dist_gaze, correlation_stride_time_dist_gaze, correlation_stride_time_var_dist_gaze, correlation_step_error_dist_gaze, correlation_step_error_var_dist_gaze, correlation_head_angle_dist_gaze, correlation_trunk_angle_dist_gaze 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
        correlation_walkingspeed_surveys_balance, correlation_walkingspeed_std_surveys_balance, correlation_metabolics_surveys_balance, correlation_step_width_surveys_balance, correlation_step_width_var_surveys_balance, correlation_step_length_surveys_balance, correlation_step_length_var_surveys_balance, correlation_stride_time_surveys_balance, correlation_stride_time_var_surveys_balance, correlation_step_error_surveys_balance, correlation_step_error_var_surveys_balance, correlation_head_angle_surveys_balance, correlation_trunk_angle_surveys_balance, correlation_dist_gaze_surveys_balance, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
        correlation_walkingspeed_surveys_foot_placement, correlation_walkingspeed_std_surveys_foot_placement, correlation_metabolics_surveys_foot_placement, correlation_step_width_surveys_foot_placement, correlation_step_width_var_surveys_foot_placement, correlation_step_length_surveys_foot_placement, correlation_step_length_var_surveys_foot_placement, correlation_stride_time_surveys_foot_placement, correlation_stride_time_var_surveys_foot_placement, correlation_step_error_surveys_foot_placement, correlation_step_error_var_surveys_foot_placement, correlation_head_angle_surveys_foot_placement, correlation_trunk_angle_surveys_foot_placement, correlation_dist_gaze_surveys_foot_placement, correlation_surveys_balance_surveys_foot_placement, 1, NaN, NaN, NaN, NaN, NaN, NaN;
        correlation_walkingspeed_surveys_walking_speed, correlation_walkingspeed_std_surveys_walking_speed, correlation_metabolics_surveys_walking_speed, correlation_step_width_surveys_walking_speed, correlation_step_width_var_surveys_walking_speed, correlation_step_length_surveys_walking_speed, correlation_step_length_var_surveys_walking_speed, correlation_stride_time_surveys_walking_speed, correlation_stride_time_var_surveys_walking_speed, correlation_step_error_surveys_walking_speed, correlation_step_error_var_surveys_walking_speed, correlation_head_angle_surveys_walking_speed, correlation_trunk_angle_surveys_walking_speed, correlation_dist_gaze_surveys_walking_speed, correlation_surveys_balance_surveys_walking_speed, correlation_surveys_foot_placement_surveys_walking_speed, 1, NaN, NaN, NaN, NaN, NaN;
        correlation_walkingspeed_surveys_metabolics, correlation_walkingspeed_std_surveys_metabolics, correlation_metabolics_surveys_metabolics, correlation_step_width_surveys_metabolics, correlation_step_width_var_surveys_metabolics, correlation_step_length_surveys_metabolics, correlation_step_length_var_surveys_metabolics, correlation_stride_time_surveys_metabolics, correlation_stride_time_var_surveys_metabolics, correlation_step_error_surveys_metabolics, correlation_step_error_var_surveys_metabolics, correlation_head_angle_surveys_metabolics, correlation_trunk_angle_surveys_metabolics, correlation_dist_gaze_surveys_metabolics, correlation_surveys_balance_surveys_metabolics, correlation_surveys_foot_placement_surveys_metabolics, correlation_surveys_walking_speed_surveys_metabolics, 1, NaN, NaN, NaN, NaN;
        correlation_walkingspeed_surveys_difficulty, correlation_walkingspeed_std_surveys_difficulty, correlation_metabolics_surveys_difficulty, correlation_step_width_surveys_difficulty, correlation_step_width_var_surveys_difficulty, correlation_step_length_surveys_difficulty, correlation_step_length_var_surveys_difficulty, correlation_stride_time_surveys_difficulty, correlation_stride_time_var_surveys_difficulty, correlation_step_error_surveys_difficulty, correlation_step_error_var_surveys_difficulty, correlation_head_angle_surveys_difficulty, correlation_trunk_angle_surveys_difficulty, correlation_dist_gaze_surveys_difficulty, correlation_surveys_balance_surveys_difficulty, correlation_surveys_foot_placement_surveys_difficulty, correlation_surveys_walking_speed_surveys_difficulty, correlation_surveys_metabolics_surveys_difficulty, 1, NaN, NaN, NaN;        
        correlation_walkingspeed_pbalance, correlation_walkingspeed_std_pbalance, correlation_metabolics_pbalance, correlation_step_width_pbalance, correlation_step_width_var_pbalance, correlation_step_length_pbalance, correlation_step_length_var_pbalance, correlation_stride_time_pbalance, correlation_stride_time_var_pbalance, correlation_step_error_pbalance, correlation_step_error_var_pbalance, correlation_head_angle_pbalance, correlation_trunk_angle_pbalance, correlation_dist_gaze_pbalance, correlation_surveys_balance_pbalance, correlation_surveys_foot_placement_pbalance, correlation_surveys_walking_speed_pbalance, correlation_surveys_metabolics_pbalance, correlation_surveys_difficulty_pbalance, NaN, NaN, NaN;
        correlation_walkingspeed_pspeed, correlation_walkingspeed_std_pspeed, correlation_metabolics_pspeed, correlation_step_width_pspeed, correlation_step_width_var_pspeed, correlation_step_length_pspeed, correlation_step_length_var_pspeed, correlation_stride_time_pspeed, correlation_stride_time_var_pspeed, correlation_step_error_pspeed, correlation_step_error_var_pspeed, correlation_head_angle_pspeed, correlation_trunk_angle_pspeed, correlation_dist_gaze_pspeed, correlation_surveys_balance_pspeed, correlation_surveys_foot_placement_pspeed, correlation_surveys_walking_speed_pspeed, correlation_surveys_metabolics_pspeed, correlation_surveys_difficulty_pspeed, NaN, NaN, NaN;
        correlation_walkingspeed_paccuracy, correlation_walkingspeed_std_paccuracy, correlation_metabolics_paccuracy, correlation_step_width_paccuracy, correlation_step_width_var_paccuracy, correlation_step_length_paccuracy, correlation_step_length_var_paccuracy, correlation_stride_time_paccuracy, correlation_stride_time_var_paccuracy, correlation_step_error_paccuracy, correlation_step_error_var_paccuracy, correlation_head_angle_paccuracy, correlation_trunk_angle_paccuracy, correlation_dist_gaze_paccuracy, correlation_surveys_balance_paccuracy, correlation_surveys_foot_placement_paccuracy, correlation_surveys_walking_speed_paccuracy, correlation_surveys_metabolics_paccuracy, correlation_surveys_difficulty_paccuracy, NaN, NaN, NaN;
    ];

    % Store in an excel file
    headings = {'Variables', 'All conditions', 'Slow Speed', 'Medium Speed', 'High Speed', 'Null Accuracy', 'Medium Accuracy', 'High Accuracy', 'No Disturbance', 'Medium Disturbance', 'High Disturbance'};

    number_columns = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,19,19]; % number of columns of each row
    variables_string = string(variables);
    variables_x_string = string(variables_x);
    count = 1;

    for i=1:size(variables_string,2)
        for j=1:number_columns(i)
            name1 = variables_string(i);
            name2 = variables_x_string(j);
            list_var(count) = name1 + " vs " + name2;
            count = count + 1;
        end
    end

    count = 1;
    for i=1:size(correlation_matrix,1)
        for j=1:size(correlation_matrix,2)
            if correlation_matrix(i,j) == 1 || isnan(correlation_matrix(i,j))
                continue
            end
            all_correlations(count) = correlation_matrix(i,j);
            count = count + 1;
        end
    end

    filtered_matrix = correlation_matrix;
    filtered_matrix(filtered_matrix == 1 | isnan(filtered_matrix)) = [];
    flattened_row = filtered_matrix(:)';
    

    writecell(headings, filename2, 'Sheet', 1, 'Range', 'A1');
    writematrix(list_var', filename2, 'Sheet', 1, 'Range', 'A2');
    writematrix(all_correlations', filename2, 'Sheet', 1, 'Range', 'B2');

    figure

    % Set the colormap with -1 as blue, 0 as white, and 1 as red
    custom_cmap = [linspace(0, 1, 128)', linspace(0, 1, 128)', ones(128, 1); % Blue to white
                   ones(128, 1), linspace(1, 0, 128)', linspace(1, 0, 128)']; % White to red
    
    imagesc(correlation_matrix, 'AlphaData', ~isnan(correlation_matrix)); % Mask NaN values
    colormap(custom_cmap); % Apply custom colormap
    colorbar; % Display the color scale
    clim([-1, 1]); % Set color axis limits
    axis square; % Make the plot square-shaped
    
    set(gca, 'Color', 'white');  % Set the background color for NaN cells to white
    
    xticks(1:length(variables_x));
    yticks(1:length(variables));
    set(gca, 'TickLength', [0.002 0]);
    xticklabels(variables_x);
    yticklabels(variables);
    
    xlim([0.5 19.5])
    set(gca, 'FontSize', 13.5);
    xtickangle(45); % Rotate the x-axis labels for better readability

    title('Linear Correlation Between Gait Variables');
    
    % Add correlation values inside the plot boxes
    for i = 1:length(variables)
        for j = 1:length(variables)
            if ~isnan(correlation_matrix(i, j)) % Only display text for non-NaN values
                if option_correlations == 1 % all
                    text(j, i, sprintf('%.2f', correlation_matrix(i, j)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'black');
                elseif option_correlations == 2 % close to +1
                    if correlation_matrix(i, j) >= 0.845 && correlation_matrix(i,j) ~= 1 % Modify this to plot only certain values
                        text(j, i, sprintf('%.2f', correlation_matrix(i, j)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'black');
                    end
                elseif option_correlations == 3 % close to -1
                    if correlation_matrix(i, j) <= -0.845
                        text(j, i, sprintf('%.2f', correlation_matrix(i, j)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'black');
                    end
                elseif option_correlations == 4 % close to 0
                    if correlation_matrix(i, j) >= -0.105 && correlation_matrix(i,j) <= 0.105
                        text(j, i, sprintf('%.2f', correlation_matrix(i, j)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'black');
                    end
                end
            end
        end
    end

    %% Plots for a given condition
    for k=1:9 % Plot for one condition: 2: s0, 3: s1, 4: s2, 5: a0, 6: a1, 7: a2, 8: b0, 9: b1, 10: b2
        if k == 1
            positions = positions_s0;
            condition = "s0";
        elseif k == 2
            positions = positions_s1;
            condition = "s1";
        elseif k == 3
            positions = positions_s2;
            condition = "s2";
        elseif k == 4
            positions = positions_a0;
            condition = "a0";
        elseif k == 5
            positions = positions_a1;
            condition = "a1";
        elseif k == 6
            positions = positions_a2;
            condition = "a2";
        elseif k == 7
            positions = positions_b0;
            condition = "b0";
        elseif k == 8
            positions = positions_b1;
            condition = "b1";
        elseif k == 9
            positions = positions_b2;
            condition = "b2";
        end

        % Take only the labels that we want
        % walking_speed_vector = data_ordered_bylabel.walking_speed(positions, :);
        % walking_speed_vector = walking_speed_vector(:);
        % 
        % walking_speed_std_vector = data_ordered_bylabel.walking_speed_std(positions, :);
        % walking_speed_std_vector = walking_speed_std_vector(:);
        % 
        % metabolics_vector = data_ordered_bylabel.metabolics.normalized(positions, :);
        % metabolics_vector = metabolics_vector(:);
        % 
        % step_width_vector = data_ordered_bylabel.step_width(positions, :);
        % step_width_vector = step_width_vector(:);
        % 
        % step_width_var_vector = data_ordered_bylabel.step_width_var(positions, :);
        % step_width_var_vector = step_width_var_vector(:);
        % 
        % step_length_vector = data_ordered_bylabel.step_length(positions, :);
        % step_length_vector = step_length_vector(:);
        % 
        % step_length_var_vector = data_ordered_bylabel.step_length_var(positions, :);
        % step_length_var_vector = step_length_var_vector(:);
        % 
        % stride_time_vector = data_ordered_bylabel.stride_time(positions, :);
        % stride_time_vector = stride_time_vector(:);
        % 
        % stride_time_var_vector = data_ordered_bylabel.stride_time_var(positions, :);
        % stride_time_var_vector = stride_time_var_vector(:);
        % 
        % step_error_vector = data_ordered_bylabel.step_error(positions, :);
        % step_error_vector = step_error_vector(:);
        % 
        % step_error_var_vector = data_ordered_bylabel.step_error_std(positions, :);
        % step_error_var_vector = step_error_var_vector(:);
        % 
        % head_angle_vector = data_ordered_bylabel.head_angle(positions, :);
        % head_angle_vector = head_angle_vector(:);
        % 
        % trunk_angle_vector = data_ordered_bylabel.trunk_angle(positions, :);
        % trunk_angle_vector = trunk_angle_vector(:);
        % 
        % dist_gaze_vector = data_ordered_bylabel.dist_gaze(positions, :);
        % dist_gaze_vector = dist_gaze_vector(:);
        % 
        % surveys_balance_vector = data_ordered_bylabel.surveys_balance(positions, :);
        % surveys_balance_vector = surveys_balance_vector(:);
        % 
        % surveys_foot_placement_vector = data_ordered_bylabel.surveys_foot_placement(positions, :);
        % surveys_foot_placement_vector = surveys_foot_placement_vector(:);
        % 
        % surveys_walking_speed_vector = data_ordered_bylabel.surveys_walking_speed(positions, :);
        % surveys_walking_speed_vector = surveys_walking_speed_vector(:);
        % 
        % surveys_metabolics_vector = data_ordered_bylabel.surveys_metabolics(positions, :);
        % surveys_metabolics_vector = surveys_metabolics_vector(:);
        % 
        % surveys_difficulty_vector = data_ordered_bylabel.surveys_difficulty(positions, :);
        % surveys_difficulty_vector = surveys_difficulty_vector(:);

        % Average to eliminate noise
        walking_speed_vector = mean(data_ordered_bylabel.walking_speed(positions, :), 2);
        walking_speed_vector = walking_speed_vector(:);

        walking_speed_std_vector = mean(data_ordered_bylabel.walking_speed_std(positions, :), 2);
        walking_speed_std_vector = walking_speed_std_vector(:);

        metabolics_vector = mean(data_ordered_bylabel.metabolics.normalized(positions, :), 2);
        metabolics_vector = metabolics_vector(:);

        step_width_vector = mean(data_ordered_bylabel.step_width(positions, :), 2);
        step_width_vector = step_width_vector(:);

        step_width_var_vector = mean(data_ordered_bylabel.step_width_var(positions, :), 2);
        step_width_var_vector = step_width_var_vector(:);

        step_length_vector = mean(data_ordered_bylabel.step_length(positions, :), 2);
        step_length_vector = step_length_vector(:);

        step_length_var_vector = mean(data_ordered_bylabel.step_length_var(positions, :), 2);
        step_length_var_vector = step_length_var_vector(:);

        stride_time_vector = mean(data_ordered_bylabel.stride_time(positions, :), 2);
        stride_time_vector = stride_time_vector(:);

        stride_time_var_vector = mean(data_ordered_bylabel.stride_time_var(positions, :), 2);
        stride_time_var_vector = stride_time_var_vector(:);

        step_error_vector = mean(data_ordered_bylabel.step_error(positions, :), 2);
        step_error_vector = step_error_vector(:);

        step_error_var_vector = mean(data_ordered_bylabel.step_error_std(positions, :), 2);
        step_error_var_vector = step_error_var_vector(:);

        head_angle_vector = mean(data_ordered_bylabel.head_angle(positions, :), 2);
        head_angle_vector = head_angle_vector(:);

        trunk_angle_vector = mean(data_ordered_bylabel.trunk_angle(positions, :), 2);
        trunk_angle_vector = trunk_angle_vector(:);

        dist_gaze_vector = mean(data_ordered_bylabel.dist_gaze(positions, :), 2);
        dist_gaze_vector = dist_gaze_vector(:);

        surveys_balance_vector = mean(data_ordered_bylabel.surveys_balance(positions, :), 2);
        surveys_balance_vector = surveys_balance_vector(:);

        surveys_foot_placement_vector = mean(data_ordered_bylabel.surveys_foot_placement(positions, :), 2);
        surveys_foot_placement_vector = surveys_foot_placement_vector(:);

        surveys_walking_speed_vector = mean(data_ordered_bylabel.surveys_walking_speed(positions, :), 2);
        surveys_walking_speed_vector = surveys_walking_speed_vector(:);

        surveys_metabolics_vector = mean(data_ordered_bylabel.surveys_metabolics(positions, :), 2);
        surveys_metabolics_vector = surveys_metabolics_vector(:);

        surveys_difficulty_vector = mean(data_ordered_bylabel.surveys_difficulty(positions, :), 2);
        surveys_difficulty_vector = surveys_difficulty_vector(:);
        
        % Prompts
        positions_s0_filtered = intersect(positions, positions_s0);
        positions_s1_filtered = intersect(positions, positions_s1);
        positions_s2_filtered = intersect(positions, positions_s2);
        positions_a0_filtered = intersect(positions, positions_a0);
        positions_a1_filtered = intersect(positions, positions_a1);
        positions_a2_filtered = intersect(positions, positions_a2);
        positions_b0_filtered = intersect(positions, positions_b0);
        positions_b1_filtered = intersect(positions, positions_b1);
        positions_b2_filtered = intersect(positions, positions_b2);

        walking_speed_byprompt.s0 = data_ordered_bylabel.walking_speed(positions_s0_filtered, :);
        walking_speed_byprompt.s1 = data_ordered_bylabel.walking_speed(positions_s1_filtered, :);
        walking_speed_byprompt.s2 = data_ordered_bylabel.walking_speed(positions_s2_filtered, :);
        walking_speed_byprompt.a0 = data_ordered_bylabel.walking_speed(positions_a0_filtered, :);
        walking_speed_byprompt.a1 = data_ordered_bylabel.walking_speed(positions_a1_filtered, :);
        walking_speed_byprompt.a2 = data_ordered_bylabel.walking_speed(positions_a2_filtered, :);
        walking_speed_byprompt.b0 = data_ordered_bylabel.walking_speed(positions_b0_filtered, :);
        walking_speed_byprompt.b1 = data_ordered_bylabel.walking_speed(positions_b1_filtered, :);
        walking_speed_byprompt.b2 = data_ordered_bylabel.walking_speed(positions_b2_filtered, :);

        walking_speed_std_byprompt.s0 = data_ordered_bylabel.walking_speed_std(positions_s0_filtered, :);
        walking_speed_std_byprompt.s1 = data_ordered_bylabel.walking_speed_std(positions_s1_filtered, :);
        walking_speed_std_byprompt.s2 = data_ordered_bylabel.walking_speed_std(positions_s2_filtered, :);
        walking_speed_std_byprompt.a0 = data_ordered_bylabel.walking_speed_std(positions_a0_filtered, :);
        walking_speed_std_byprompt.a1 = data_ordered_bylabel.walking_speed_std(positions_a1_filtered, :);
        walking_speed_std_byprompt.a2 = data_ordered_bylabel.walking_speed_std(positions_a2_filtered, :);
        walking_speed_std_byprompt.b0 = data_ordered_bylabel.walking_speed_std(positions_b0_filtered, :);
        walking_speed_std_byprompt.b1 = data_ordered_bylabel.walking_speed_std(positions_b1_filtered, :);
        walking_speed_std_byprompt.b2 = data_ordered_bylabel.walking_speed_std(positions_b2_filtered, :);

        metabolics_byprompt.s0 = data_ordered_bylabel.metabolics.normalized(positions_s0_filtered, :);
        metabolics_byprompt.s1 = data_ordered_bylabel.metabolics.normalized(positions_s1_filtered, :);
        metabolics_byprompt.s2 = data_ordered_bylabel.metabolics.normalized(positions_s2_filtered, :);
        metabolics_byprompt.a0 = data_ordered_bylabel.metabolics.normalized(positions_a0_filtered, :);
        metabolics_byprompt.a1 = data_ordered_bylabel.metabolics.normalized(positions_a1_filtered, :);
        metabolics_byprompt.a2 = data_ordered_bylabel.metabolics.normalized(positions_a2_filtered, :);
        metabolics_byprompt.b0 = data_ordered_bylabel.metabolics.normalized(positions_b0_filtered, :);
        metabolics_byprompt.b1 = data_ordered_bylabel.metabolics.normalized(positions_b1_filtered, :);
        metabolics_byprompt.b2 = data_ordered_bylabel.metabolics.normalized(positions_b2_filtered, :);

        step_width_byprompt.s0 = data_ordered_bylabel.step_width(positions_s0_filtered, :);
        step_width_byprompt.s1 = data_ordered_bylabel.step_width(positions_s1_filtered, :);
        step_width_byprompt.s2 = data_ordered_bylabel.step_width(positions_s2_filtered, :);
        step_width_byprompt.a0 = data_ordered_bylabel.step_width(positions_a0_filtered, :);
        step_width_byprompt.a1 = data_ordered_bylabel.step_width(positions_a1_filtered, :);
        step_width_byprompt.a2 = data_ordered_bylabel.step_width(positions_a2_filtered, :);
        step_width_byprompt.b0 = data_ordered_bylabel.step_width(positions_b0_filtered, :);
        step_width_byprompt.b1 = data_ordered_bylabel.step_width(positions_b1_filtered, :);
        step_width_byprompt.b2 = data_ordered_bylabel.step_width(positions_b2_filtered, :);

        step_width_var_byprompt.s0 = data_ordered_bylabel.step_width_var(positions_s0_filtered, :);
        step_width_var_byprompt.s1 = data_ordered_bylabel.step_width_var(positions_s1_filtered, :);
        step_width_var_byprompt.s2 = data_ordered_bylabel.step_width_var(positions_s2_filtered, :);
        step_width_var_byprompt.a0 = data_ordered_bylabel.step_width_var(positions_a0_filtered, :);
        step_width_var_byprompt.a1 = data_ordered_bylabel.step_width_var(positions_a1_filtered, :);
        step_width_var_byprompt.a2 = data_ordered_bylabel.step_width_var(positions_a2_filtered, :);
        step_width_var_byprompt.b0 = data_ordered_bylabel.step_width_var(positions_b0_filtered, :);
        step_width_var_byprompt.b1 = data_ordered_bylabel.step_width_var(positions_b1_filtered, :);
        step_width_var_byprompt.b2 = data_ordered_bylabel.step_width_var(positions_b2_filtered, :);

        step_length_byprompt.s0 = data_ordered_bylabel.step_length(positions_s0_filtered, :);
        step_length_byprompt.s1 = data_ordered_bylabel.step_length(positions_s1_filtered, :);
        step_length_byprompt.s2 = data_ordered_bylabel.step_length(positions_s2_filtered, :);
        step_length_byprompt.a0 = data_ordered_bylabel.step_length(positions_a0_filtered, :);
        step_length_byprompt.a1 = data_ordered_bylabel.step_length(positions_a1_filtered, :);
        step_length_byprompt.a2 = data_ordered_bylabel.step_length(positions_a2_filtered, :);
        step_length_byprompt.b0 = data_ordered_bylabel.step_length(positions_b0_filtered, :);
        step_length_byprompt.b1 = data_ordered_bylabel.step_length(positions_b1_filtered, :);
        step_length_byprompt.b2 = data_ordered_bylabel.step_length(positions_b2_filtered, :);

        step_length_var_byprompt.s0 = data_ordered_bylabel.step_length_var(positions_s0_filtered, :);
        step_length_var_byprompt.s1 = data_ordered_bylabel.step_length_var(positions_s1_filtered, :);
        step_length_var_byprompt.s2 = data_ordered_bylabel.step_length_var(positions_s2_filtered, :);
        step_length_var_byprompt.a0 = data_ordered_bylabel.step_length_var(positions_a0_filtered, :);
        step_length_var_byprompt.a1 = data_ordered_bylabel.step_length_var(positions_a1_filtered, :);
        step_length_var_byprompt.a2 = data_ordered_bylabel.step_length_var(positions_a2_filtered, :);
        step_length_var_byprompt.b0 = data_ordered_bylabel.step_length_var(positions_b0_filtered, :);
        step_length_var_byprompt.b1 = data_ordered_bylabel.step_length_var(positions_b1_filtered, :);
        step_length_var_byprompt.b2 = data_ordered_bylabel.step_length_var(positions_b2_filtered, :);

        step_error_byprompt.s0 = data_ordered_bylabel.step_error(positions_s0_filtered, :);
        step_error_byprompt.s1 = data_ordered_bylabel.step_error(positions_s1_filtered, :);
        step_error_byprompt.s2 = data_ordered_bylabel.step_error(positions_s2_filtered, :);
        step_error_byprompt.a0 = data_ordered_bylabel.step_error(positions_a0_filtered, :);
        step_error_byprompt.a1 = data_ordered_bylabel.step_error(positions_a1_filtered, :);
        step_error_byprompt.a2 = data_ordered_bylabel.step_error(positions_a2_filtered, :);
        step_error_byprompt.b0 = data_ordered_bylabel.step_error(positions_b0_filtered, :);
        step_error_byprompt.b1 = data_ordered_bylabel.step_error(positions_b1_filtered, :);
        step_error_byprompt.b2 = data_ordered_bylabel.step_error(positions_b2_filtered, :);

        step_error_var_byprompt.s0 = data_ordered_bylabel.step_error_std(positions_s0_filtered, :);
        step_error_var_byprompt.s1 = data_ordered_bylabel.step_error_std(positions_s1_filtered, :);
        step_error_var_byprompt.s2 = data_ordered_bylabel.step_error_std(positions_s2_filtered, :);
        step_error_var_byprompt.a0 = data_ordered_bylabel.step_error_std(positions_a0_filtered, :);
        step_error_var_byprompt.a1 = data_ordered_bylabel.step_error_std(positions_a1_filtered, :);
        step_error_var_byprompt.a2 = data_ordered_bylabel.step_error_std(positions_a2_filtered, :);
        step_error_var_byprompt.b0 = data_ordered_bylabel.step_error_std(positions_b0_filtered, :);
        step_error_var_byprompt.b1 = data_ordered_bylabel.step_error_std(positions_b1_filtered, :);
        step_error_var_byprompt.b2 = data_ordered_bylabel.step_error_std(positions_b2_filtered, :);

        stride_time_byprompt.s0 = data_ordered_bylabel.stride_time(positions_s0_filtered, :);
        stride_time_byprompt.s1 = data_ordered_bylabel.stride_time(positions_s1_filtered, :);
        stride_time_byprompt.s2 = data_ordered_bylabel.stride_time(positions_s2_filtered, :);
        stride_time_byprompt.a0 = data_ordered_bylabel.stride_time(positions_a0_filtered, :);
        stride_time_byprompt.a1 = data_ordered_bylabel.stride_time(positions_a1_filtered, :);
        stride_time_byprompt.a2 = data_ordered_bylabel.stride_time(positions_a2_filtered, :);
        stride_time_byprompt.b0 = data_ordered_bylabel.stride_time(positions_b0_filtered, :);
        stride_time_byprompt.b1 = data_ordered_bylabel.stride_time(positions_b1_filtered, :);
        stride_time_byprompt.b2 = data_ordered_bylabel.stride_time(positions_b2_filtered, :);

        stride_time_var_byprompt.s0 = data_ordered_bylabel.stride_time_var(positions_s0_filtered, :);
        stride_time_var_byprompt.s1 = data_ordered_bylabel.stride_time_var(positions_s1_filtered, :);
        stride_time_var_byprompt.s2 = data_ordered_bylabel.stride_time_var(positions_s2_filtered, :);
        stride_time_var_byprompt.a0 = data_ordered_bylabel.stride_time_var(positions_a0_filtered, :);
        stride_time_var_byprompt.a1 = data_ordered_bylabel.stride_time_var(positions_a1_filtered, :);
        stride_time_var_byprompt.a2 = data_ordered_bylabel.stride_time_var(positions_a2_filtered, :);
        stride_time_var_byprompt.b0 = data_ordered_bylabel.stride_time_var(positions_b0_filtered, :);
        stride_time_var_byprompt.b1 = data_ordered_bylabel.stride_time_var(positions_b1_filtered, :);
        stride_time_var_byprompt.b2 = data_ordered_bylabel.stride_time_var(positions_b2_filtered, :);

        head_angle_byprompt.s0 = data_ordered_bylabel.head_angle(positions_s0_filtered, :);
        head_angle_byprompt.s1 = data_ordered_bylabel.head_angle(positions_s1_filtered, :);
        head_angle_byprompt.s2 = data_ordered_bylabel.head_angle(positions_s2_filtered, :);
        head_angle_byprompt.a0 = data_ordered_bylabel.head_angle(positions_a0_filtered, :);
        head_angle_byprompt.a1 = data_ordered_bylabel.head_angle(positions_a1_filtered, :);
        head_angle_byprompt.a2 = data_ordered_bylabel.head_angle(positions_a2_filtered, :);
        head_angle_byprompt.b0 = data_ordered_bylabel.head_angle(positions_b0_filtered, :);
        head_angle_byprompt.b1 = data_ordered_bylabel.head_angle(positions_b1_filtered, :);
        head_angle_byprompt.b2 = data_ordered_bylabel.head_angle(positions_b2_filtered, :);

        trunk_angle_byprompt.s0 = data_ordered_bylabel.trunk_angle(positions_s0_filtered, :);
        trunk_angle_byprompt.s1 = data_ordered_bylabel.trunk_angle(positions_s1_filtered, :);
        trunk_angle_byprompt.s2 = data_ordered_bylabel.trunk_angle(positions_s2_filtered, :);
        trunk_angle_byprompt.a0 = data_ordered_bylabel.trunk_angle(positions_a0_filtered, :);
        trunk_angle_byprompt.a1 = data_ordered_bylabel.trunk_angle(positions_a1_filtered, :);
        trunk_angle_byprompt.a2 = data_ordered_bylabel.trunk_angle(positions_a2_filtered, :);
        trunk_angle_byprompt.b0 = data_ordered_bylabel.trunk_angle(positions_b0_filtered, :);
        trunk_angle_byprompt.b1 = data_ordered_bylabel.trunk_angle(positions_b1_filtered, :);
        trunk_angle_byprompt.b2 = data_ordered_bylabel.trunk_angle(positions_b2_filtered, :);

        dist_gaze_byprompt.s0 = data_ordered_bylabel.dist_gaze(positions_s0_filtered, :);
        dist_gaze_byprompt.s1 = data_ordered_bylabel.dist_gaze(positions_s1_filtered, :);
        dist_gaze_byprompt.s2 = data_ordered_bylabel.dist_gaze(positions_s2_filtered, :);
        dist_gaze_byprompt.a0 = data_ordered_bylabel.dist_gaze(positions_a0_filtered, :);
        dist_gaze_byprompt.a1 = data_ordered_bylabel.dist_gaze(positions_a1_filtered, :);
        dist_gaze_byprompt.a2 = data_ordered_bylabel.dist_gaze(positions_a2_filtered, :);
        dist_gaze_byprompt.b0 = data_ordered_bylabel.dist_gaze(positions_b0_filtered, :);
        dist_gaze_byprompt.b1 = data_ordered_bylabel.dist_gaze(positions_b1_filtered, :);
        dist_gaze_byprompt.b2 = data_ordered_bylabel.dist_gaze(positions_b2_filtered, :);

        surveys_balance_byprompt.s0 = data_ordered_bylabel.surveys_balance(positions_s0_filtered, :);
        surveys_balance_byprompt.s1 = data_ordered_bylabel.surveys_balance(positions_s1_filtered, :);
        surveys_balance_byprompt.s2 = data_ordered_bylabel.surveys_balance(positions_s2_filtered, :);
        surveys_balance_byprompt.a0 = data_ordered_bylabel.surveys_balance(positions_a0_filtered, :);
        surveys_balance_byprompt.a1 = data_ordered_bylabel.surveys_balance(positions_a1_filtered, :);
        surveys_balance_byprompt.a2 = data_ordered_bylabel.surveys_balance(positions_a2_filtered, :);
        surveys_balance_byprompt.b0 = data_ordered_bylabel.surveys_balance(positions_b0_filtered, :);
        surveys_balance_byprompt.b1 = data_ordered_bylabel.surveys_balance(positions_b1_filtered, :);
        surveys_balance_byprompt.b2 = data_ordered_bylabel.surveys_balance(positions_b2_filtered, :);

        surveys_foot_placement_byprompt.s0 = data_ordered_bylabel.surveys_foot_placement(positions_s0_filtered, :);
        surveys_foot_placement_byprompt.s1 = data_ordered_bylabel.surveys_foot_placement(positions_s1_filtered, :);
        surveys_foot_placement_byprompt.s2 = data_ordered_bylabel.surveys_foot_placement(positions_s2_filtered, :);
        surveys_foot_placement_byprompt.a0 = data_ordered_bylabel.surveys_foot_placement(positions_a0_filtered, :);
        surveys_foot_placement_byprompt.a1 = data_ordered_bylabel.surveys_foot_placement(positions_a1_filtered, :);
        surveys_foot_placement_byprompt.a2 = data_ordered_bylabel.surveys_foot_placement(positions_a2_filtered, :);
        surveys_foot_placement_byprompt.b0 = data_ordered_bylabel.surveys_foot_placement(positions_b0_filtered, :);
        surveys_foot_placement_byprompt.b1 = data_ordered_bylabel.surveys_foot_placement(positions_b1_filtered, :);
        surveys_foot_placement_byprompt.b2 = data_ordered_bylabel.surveys_foot_placement(positions_b2_filtered, :);

        surveys_walking_speed_byprompt.s0 = data_ordered_bylabel.surveys_walking_speed(positions_s0_filtered, :);
        surveys_walking_speed_byprompt.s1 = data_ordered_bylabel.surveys_walking_speed(positions_s1_filtered, :);
        surveys_walking_speed_byprompt.s2 = data_ordered_bylabel.surveys_walking_speed(positions_s2_filtered, :);
        surveys_walking_speed_byprompt.a0 = data_ordered_bylabel.surveys_walking_speed(positions_a0_filtered, :);
        surveys_walking_speed_byprompt.a1 = data_ordered_bylabel.surveys_walking_speed(positions_a1_filtered, :);
        surveys_walking_speed_byprompt.a2 = data_ordered_bylabel.surveys_walking_speed(positions_a2_filtered, :);
        surveys_walking_speed_byprompt.b0 = data_ordered_bylabel.surveys_walking_speed(positions_b0_filtered, :);
        surveys_walking_speed_byprompt.b1 = data_ordered_bylabel.surveys_walking_speed(positions_b1_filtered, :);
        surveys_walking_speed_byprompt.b2 = data_ordered_bylabel.surveys_walking_speed(positions_b2_filtered, :);

        surveys_metabolics_byprompt.s0 = data_ordered_bylabel.surveys_metabolics(positions_s0_filtered, :);
        surveys_metabolics_byprompt.s1 = data_ordered_bylabel.surveys_metabolics(positions_s1_filtered, :);
        surveys_metabolics_byprompt.s2 = data_ordered_bylabel.surveys_metabolics(positions_s2_filtered, :);
        surveys_metabolics_byprompt.a0 = data_ordered_bylabel.surveys_metabolics(positions_a0_filtered, :);
        surveys_metabolics_byprompt.a1 = data_ordered_bylabel.surveys_metabolics(positions_a1_filtered, :);
        surveys_metabolics_byprompt.a2 = data_ordered_bylabel.surveys_metabolics(positions_a2_filtered, :);
        surveys_metabolics_byprompt.b0 = data_ordered_bylabel.surveys_metabolics(positions_b0_filtered, :);
        surveys_metabolics_byprompt.b1 = data_ordered_bylabel.surveys_metabolics(positions_b1_filtered, :);
        surveys_metabolics_byprompt.b2 = data_ordered_bylabel.surveys_metabolics(positions_b2_filtered, :);

        surveys_difficulty_byprompt.s0 = data_ordered_bylabel.surveys_difficulty(positions_s0_filtered, :);
        surveys_difficulty_byprompt.s1 = data_ordered_bylabel.surveys_difficulty(positions_s1_filtered, :);
        surveys_difficulty_byprompt.s2 = data_ordered_bylabel.surveys_difficulty(positions_s2_filtered, :);
        surveys_difficulty_byprompt.a0 = data_ordered_bylabel.surveys_difficulty(positions_a0_filtered, :);
        surveys_difficulty_byprompt.a1 = data_ordered_bylabel.surveys_difficulty(positions_a1_filtered, :);
        surveys_difficulty_byprompt.a2 = data_ordered_bylabel.surveys_difficulty(positions_a2_filtered, :);
        surveys_difficulty_byprompt.b0 = data_ordered_bylabel.surveys_difficulty(positions_b0_filtered, :);
        surveys_difficulty_byprompt.b1 = data_ordered_bylabel.surveys_difficulty(positions_b1_filtered, :);
        surveys_difficulty_byprompt.b2 = data_ordered_bylabel.surveys_difficulty(positions_b2_filtered, :);

       % Calculate the correlations: same as in general plot
        correlation_walkingspeed_walkingspeed_std = corr(walking_speed_vector, walking_speed_std_vector);
        correlation_walkingspeed_metabolics = corr(walking_speed_vector, metabolics_vector);
        correlation_walkingspeed_step_width = corr(walking_speed_vector, step_width_vector);
        correlation_walkingspeed_step_width_var = corr(walking_speed_vector, step_width_var_vector);
        correlation_walkingspeed_step_length = corr(walking_speed_vector, step_length_vector);
        correlation_walkingspeed_step_length_var = corr(walking_speed_vector, step_length_var_vector);
        correlation_walkingspeed_stride_time = corr(walking_speed_vector, stride_time_vector);
        correlation_walkingspeed_stride_time_var = corr(walking_speed_vector, stride_time_var_vector);
        correlation_walkingspeed_step_error = corr(walking_speed_vector, step_error_vector);
        correlation_walkingspeed_step_error_var = corr(walking_speed_vector, step_error_var_vector);
        correlation_walkingspeed_head_angle = corr(walking_speed_vector, head_angle_vector);
        correlation_walkingspeed_trunk_angle = corr(walking_speed_vector, trunk_angle_vector);
        correlation_walkingspeed_dist_gaze = corr(walking_speed_vector, dist_gaze_vector);
        correlation_walkingspeed_surveys_balance = corr(walking_speed_vector, surveys_balance_vector);
        correlation_walkingspeed_surveys_foot_placement = corr(walking_speed_vector, surveys_foot_placement_vector);
        correlation_walkingspeed_surveys_walking_speed = corr(walking_speed_vector, surveys_walking_speed_vector);
        correlation_walkingspeed_surveys_metabolics = corr(walking_speed_vector, surveys_metabolics_vector);
        correlation_walkingspeed_surveys_difficulty = corr(walking_speed_vector, surveys_difficulty_vector);
        correlation_walkingspeed_pspeed = Correlation_prompts(walking_speed_byprompt, "s0", "s1", "s2");
        correlation_walkingspeed_paccuracy = Correlation_prompts(walking_speed_byprompt, "a0", "a1", "a2");
        correlation_walkingspeed_pbalance = Correlation_prompts(walking_speed_byprompt, "b0", "b1", "b2");
    
        correlation_walkingspeed_std_metabolics = corr(walking_speed_std_vector, metabolics_vector);
        correlation_walkingspeed_std_step_width = corr(walking_speed_std_vector, step_width_vector);
        correlation_walkingspeed_std_step_width_var = corr(walking_speed_std_vector, step_width_var_vector);
        correlation_walkingspeed_std_step_length = corr(walking_speed_std_vector, step_length_vector);
        correlation_walkingspeed_std_step_length_var = corr(walking_speed_std_vector, step_length_var_vector);
        correlation_walkingspeed_std_stride_time = corr(walking_speed_std_vector, stride_time_vector);
        correlation_walkingspeed_std_stride_time_var = corr(walking_speed_std_vector, stride_time_var_vector);
        correlation_walkingspeed_std_step_error = corr(walking_speed_std_vector, step_error_vector);
        correlation_walkingspeed_std_step_error_var = corr(walking_speed_std_vector, step_error_var_vector);
        correlation_walkingspeed_std_head_angle = corr(walking_speed_std_vector, head_angle_vector);
        correlation_walkingspeed_std_trunk_angle = corr(walking_speed_std_vector, trunk_angle_vector);
        correlation_walkingspeed_std_dist_gaze = corr(walking_speed_std_vector, dist_gaze_vector);
        correlation_walkingspeed_std_surveys_balance = corr(walking_speed_std_vector, surveys_balance_vector);
        correlation_walkingspeed_std_surveys_foot_placement = corr(walking_speed_std_vector, surveys_foot_placement_vector);
        correlation_walkingspeed_std_surveys_walking_speed = corr(walking_speed_std_vector, surveys_walking_speed_vector);
        correlation_walkingspeed_std_surveys_metabolics = corr(walking_speed_std_vector, surveys_metabolics_vector);
        correlation_walkingspeed_std_surveys_difficulty = corr(walking_speed_std_vector, surveys_difficulty_vector);
        correlation_walkingspeed_std_pspeed = Correlation_prompts(walking_speed_std_byprompt, "s0", "s1", "s2");
        correlation_walkingspeed_std_paccuracy = Correlation_prompts(walking_speed_std_byprompt, "a0", "a1", "a2");
        correlation_walkingspeed_std_pbalance = Correlation_prompts(walking_speed_std_byprompt, "b0", "b1", "b2");
    
        correlation_metabolics_step_width = corr(metabolics_vector,step_width_vector);
        correlation_metabolics_step_width_var = corr(metabolics_vector, step_width_var_vector);
        correlation_metabolics_step_length = corr(metabolics_vector, step_length_vector);
        correlation_metabolics_step_length_var = corr(metabolics_vector, step_length_var_vector);
        correlation_metabolics_stride_time = corr(metabolics_vector, stride_time_vector);
        correlation_metabolics_stride_time_var = corr(metabolics_vector, stride_time_var_vector);
        correlation_metabolics_step_error = corr(metabolics_vector, step_error_vector);
        correlation_metabolics_step_error_var = corr(metabolics_vector, step_error_var_vector);
        correlation_metabolics_head_angle = corr(metabolics_vector, head_angle_vector);
        correlation_metabolics_trunk_angle = corr(metabolics_vector, trunk_angle_vector);
        correlation_metabolics_dist_gaze = corr(metabolics_vector, dist_gaze_vector);
        correlation_metabolics_surveys_balance = corr(metabolics_vector, surveys_balance_vector);
        correlation_metabolics_surveys_foot_placement = corr(metabolics_vector, surveys_foot_placement_vector);
        correlation_metabolics_surveys_walking_speed = corr(metabolics_vector, surveys_walking_speed_vector);
        correlation_metabolics_surveys_metabolics = corr(metabolics_vector, surveys_metabolics_vector);
        correlation_metabolics_surveys_difficulty = corr(metabolics_vector, surveys_difficulty_vector);
        correlation_metabolics_pspeed = Correlation_prompts(metabolics_byprompt, "s0", "s1", "s2");
        correlation_metabolics_paccuracy = Correlation_prompts(metabolics_byprompt, "a0", "a1", "a2");
        correlation_metabolics_pbalance = Correlation_prompts(metabolics_byprompt, "b0", "b1", "b2");
    
        correlation_step_width_step_width_var = corr(step_width_vector, step_width_var_vector);
        correlation_step_width_step_length = corr(step_width_vector, step_length_vector);
        correlation_step_width_step_length_var = corr(step_width_vector, step_length_var_vector);
        correlation_step_width_stride_time = corr(step_width_vector, stride_time_vector);
        correlation_step_width_stride_time_var = corr(step_width_vector, stride_time_var_vector);
        correlation_step_width_step_error = corr(step_width_vector, step_error_vector);
        correlation_step_width_step_error_var = corr(step_width_vector, step_error_var_vector);
        correlation_step_width_head_angle = corr(step_width_vector, head_angle_vector);
        correlation_step_width_trunk_angle = corr(step_width_vector, trunk_angle_vector);
        correlation_step_width_dist_gaze = corr(step_width_vector, dist_gaze_vector);
        correlation_step_width_surveys_balance = corr(step_width_vector, surveys_balance_vector);
        correlation_step_width_surveys_foot_placement = corr(step_width_vector, surveys_foot_placement_vector);
        correlation_step_width_surveys_walking_speed = corr(step_width_vector, surveys_walking_speed_vector);
        correlation_step_width_surveys_metabolics = corr(step_width_vector, surveys_metabolics_vector);
        correlation_step_width_surveys_difficulty = corr(step_width_vector, surveys_difficulty_vector);
        correlation_step_width_pspeed = Correlation_prompts(step_width_byprompt, "s0", "s1", "s2");
        correlation_step_width_paccuracy = Correlation_prompts(step_width_byprompt, "a0", "a1", "a2");
        correlation_step_width_pbalance = Correlation_prompts(step_width_byprompt, "b0", "b1", "b2");
    
        correlation_step_width_var_step_length = corr(step_width_var_vector, step_length_vector);
        correlation_step_width_var_step_length_var = corr(step_width_var_vector, step_length_var_vector);
        correlation_step_width_var_stride_time = corr(step_width_var_vector, stride_time_vector);
        correlation_step_width_var_stride_time_var = corr(step_width_var_vector, stride_time_var_vector);
        correlation_step_width_var_step_error = corr(step_width_var_vector, step_error_vector);
        correlation_step_width_var_step_error_var = corr(step_width_var_vector, step_error_var_vector);
        correlation_step_width_var_head_angle = corr(step_width_var_vector, head_angle_vector);
        correlation_step_width_var_trunk_angle = corr(step_width_var_vector, trunk_angle_vector);
        correlation_step_width_var_dist_gaze = corr(step_width_var_vector, dist_gaze_vector);
        correlation_step_width_var_surveys_balance = corr(step_width_var_vector, surveys_balance_vector);
        correlation_step_width_var_surveys_foot_placement = corr(step_width_var_vector, surveys_foot_placement_vector);
        correlation_step_width_var_surveys_walking_speed = corr(step_width_var_vector, surveys_walking_speed_vector);
        correlation_step_width_var_surveys_metabolics = corr(step_width_var_vector, surveys_metabolics_vector);
        correlation_step_width_var_surveys_difficulty = corr(step_width_var_vector, surveys_difficulty_vector);
        correlation_step_width_var_pspeed = Correlation_prompts(step_width_var_byprompt, "s0", "s1", "s2");
        correlation_step_width_var_paccuracy = Correlation_prompts(step_width_var_byprompt, "a0", "a1", "a2");
        correlation_step_width_var_pbalance = Correlation_prompts(step_width_var_byprompt, "b0", "b1", "b2");
    
        correlation_step_length_step_length_var = corr(step_length_vector, step_length_var_vector);
        correlation_step_length_stride_time = corr(step_length_vector, stride_time_vector);
        correlation_step_length_stride_time_var = corr(step_length_vector, stride_time_var_vector);
        correlation_step_length_step_error = corr(step_length_vector, step_error_vector);
        correlation_step_length_step_error_var = corr(step_length_vector, step_error_var_vector);
        correlation_step_length_head_angle = corr(step_length_vector, head_angle_vector);
        correlation_step_length_trunk_angle = corr(step_length_vector, trunk_angle_vector);
        correlation_step_length_dist_gaze = corr(step_length_vector, dist_gaze_vector);
        correlation_step_length_surveys_balance = corr(step_length_vector, surveys_balance_vector);
        correlation_step_length_surveys_foot_placement = corr(step_length_vector, surveys_foot_placement_vector);
        correlation_step_length_surveys_walking_speed = corr(step_length_vector, surveys_walking_speed_vector);
        correlation_step_length_surveys_metabolics = corr(step_length_vector, surveys_metabolics_vector);
        correlation_step_length_surveys_difficulty = corr(step_length_vector, surveys_difficulty_vector);
        correlation_step_length_pspeed = Correlation_prompts(step_length_byprompt, "s0", "s1", "s2");
        correlation_step_length_paccuracy = Correlation_prompts(step_length_byprompt, "a0", "a1", "a2");
        correlation_step_length_pbalance = Correlation_prompts(step_length_byprompt, "b0", "b1", "b2");
    
        correlation_step_length_var_stride_time = corr(step_length_var_vector, stride_time_vector);
        correlation_step_length_var_stride_time_var = corr(step_length_var_vector, stride_time_var_vector);
        correlation_step_length_var_step_error = corr(step_length_var_vector, step_error_vector);
        correlation_step_length_var_step_error_var = corr(step_length_var_vector, step_error_var_vector);
        correlation_step_length_var_head_angle = corr(step_length_var_vector, head_angle_vector);
        correlation_step_length_var_trunk_angle = corr(step_length_var_vector, trunk_angle_vector);
        correlation_step_length_var_dist_gaze = corr(step_length_var_vector, dist_gaze_vector);
        correlation_step_length_var_surveys_balance = corr(step_length_var_vector, surveys_balance_vector);
        correlation_step_length_var_surveys_foot_placement = corr(step_length_var_vector, surveys_foot_placement_vector);
        correlation_step_length_var_surveys_walking_speed = corr(step_length_var_vector, surveys_walking_speed_vector);
        correlation_step_length_var_surveys_metabolics = corr(step_length_var_vector, surveys_metabolics_vector);
        correlation_step_length_var_surveys_difficulty = corr(step_length_var_vector, surveys_difficulty_vector);
        correlation_step_length_var_pspeed = Correlation_prompts(step_length_var_byprompt, "s0", "s1", "s2");
        correlation_step_length_var_paccuracy = Correlation_prompts(step_length_var_byprompt, "a0", "a1", "a2");
        correlation_step_length_var_pbalance = Correlation_prompts(step_length_var_byprompt, "b0", "b1", "b2");
    
        correlation_stride_time_stride_time_var = corr(stride_time_vector, stride_time_var_vector);
        correlation_stride_time_step_error = corr(stride_time_vector, step_error_vector);
        correlation_stride_time_step_error_var = corr(stride_time_vector, step_error_var_vector);
        correlation_stride_time_head_angle = corr(stride_time_vector, head_angle_vector);
        correlation_stride_time_trunk_angle = corr(stride_time_vector, trunk_angle_vector);
        correlation_stride_time_dist_gaze = corr(stride_time_vector, dist_gaze_vector);
        correlation_stride_time_surveys_balance = corr(stride_time_vector, surveys_balance_vector);
        correlation_stride_time_surveys_foot_placement = corr(stride_time_vector, surveys_foot_placement_vector);
        correlation_stride_time_surveys_walking_speed = corr(stride_time_vector, surveys_walking_speed_vector);
        correlation_stride_time_surveys_metabolics = corr(stride_time_vector, surveys_metabolics_vector);
        correlation_stride_time_surveys_difficulty = corr(stride_time_vector, surveys_difficulty_vector);
        correlation_stride_time_pspeed = Correlation_prompts(stride_time_byprompt, "s0", "s1", "s2");
        correlation_stride_time_paccuracy = Correlation_prompts(stride_time_byprompt, "a0", "a1", "a2");
        correlation_stride_time_pbalance = Correlation_prompts(stride_time_byprompt, "b0", "b1", "b2");
    
        correlation_stride_time_var_step_error = corr(stride_time_var_vector, step_error_vector);
        correlation_stride_time_var_step_error_var = corr(stride_time_var_vector, step_error_var_vector);
        correlation_stride_time_var_head_angle = corr(stride_time_var_vector, head_angle_vector);
        correlation_stride_time_var_trunk_angle = corr(stride_time_var_vector, trunk_angle_vector);
        correlation_stride_time_var_dist_gaze = corr(stride_time_var_vector, dist_gaze_vector);
        correlation_stride_time_var_surveys_balance = corr(stride_time_var_vector, surveys_balance_vector);
        correlation_stride_time_var_surveys_foot_placement = corr(stride_time_var_vector, surveys_foot_placement_vector);
        correlation_stride_time_var_surveys_walking_speed = corr(stride_time_var_vector, surveys_walking_speed_vector);
        correlation_stride_time_var_surveys_metabolics = corr(stride_time_var_vector, surveys_metabolics_vector);
        correlation_stride_time_var_surveys_difficulty = corr(stride_time_var_vector, surveys_difficulty_vector);
        correlation_stride_time_var_pspeed = Correlation_prompts(stride_time_var_byprompt, "s0", "s1", "s2");
        correlation_stride_time_var_paccuracy = Correlation_prompts(stride_time_var_byprompt, "a0", "a1", "a2");
        correlation_stride_time_var_pbalance = Correlation_prompts(stride_time_var_byprompt, "b0", "b1", "b2");
    
        correlation_step_error_step_error_var = corr(step_error_vector, step_error_var_vector);
        correlation_step_error_head_angle = corr(step_error_vector, head_angle_vector);
        correlation_step_error_trunk_angle = corr(step_error_vector, trunk_angle_vector);
        correlation_step_error_dist_gaze = corr(step_error_vector, dist_gaze_vector);
        correlation_step_error_surveys_balance = corr(step_error_vector, surveys_balance_vector);
        correlation_step_error_surveys_foot_placement = corr(step_error_vector, surveys_foot_placement_vector);
        correlation_step_error_surveys_walking_speed = corr(step_error_vector, surveys_walking_speed_vector);
        correlation_step_error_surveys_metabolics = corr(step_error_vector, surveys_metabolics_vector);
        correlation_step_error_surveys_difficulty = corr(step_error_vector, surveys_difficulty_vector);
        correlation_step_error_pspeed = Correlation_prompts(step_error_byprompt, "s0", "s1", "s2");
        correlation_step_error_paccuracy = Correlation_prompts(step_error_byprompt, "a0", "a1", "a2");
        correlation_step_error_pbalance = Correlation_prompts(step_error_byprompt, "b0", "b1", "b2");
    
        correlation_step_error_var_head_angle = corr(step_error_var_vector, head_angle_vector);
        correlation_step_error_var_trunk_angle = corr(step_error_var_vector, trunk_angle_vector);
        correlation_step_error_var_dist_gaze = corr(step_error_var_vector, dist_gaze_vector);
        correlation_step_error_var_surveys_balance = corr(step_error_var_vector, surveys_balance_vector);
        correlation_step_error_var_surveys_foot_placement = corr(step_error_var_vector, surveys_foot_placement_vector);
        correlation_step_error_var_surveys_walking_speed = corr(step_error_var_vector, surveys_walking_speed_vector);
        correlation_step_error_var_surveys_metabolics = corr(step_error_var_vector, surveys_metabolics_vector);
        correlation_step_error_var_surveys_difficulty = corr(step_error_var_vector, surveys_difficulty_vector);
        correlation_step_error_var_pspeed = Correlation_prompts(step_error_var_byprompt, "s0", "s1", "s2");
        correlation_step_error_var_paccuracy = Correlation_prompts(step_error_var_byprompt, "a0", "a1", "a2");
        correlation_step_error_var_pbalance = Correlation_prompts(step_error_var_byprompt, "b0", "b1", "b2");
    
        correlation_head_angle_trunk_angle = corr(head_angle_vector, trunk_angle_vector);
        correlation_head_angle_dist_gaze = corr(head_angle_vector, dist_gaze_vector);
        correlation_head_angle_surveys_balance = corr(head_angle_vector, surveys_balance_vector);
        correlation_head_angle_surveys_foot_placement = corr(head_angle_vector, surveys_foot_placement_vector);
        correlation_head_angle_surveys_walking_speed = corr(head_angle_vector, surveys_walking_speed_vector);
        correlation_head_angle_surveys_metabolics = corr(head_angle_vector, surveys_metabolics_vector);
        correlation_head_angle_surveys_difficulty = corr(head_angle_vector, surveys_difficulty_vector);
        correlation_head_angle_pspeed = Correlation_prompts(head_angle_byprompt, "s0", "s1", "s2");
        correlation_head_angle_paccuracy = Correlation_prompts(head_angle_byprompt, "a0", "a1", "a2");
        correlation_head_angle_pbalance = Correlation_prompts(head_angle_byprompt, "b0", "b1", "b2");
    
        correlation_trunk_angle_dist_gaze = corr(trunk_angle_vector, dist_gaze_vector);
        correlation_trunk_angle_surveys_balance = corr(trunk_angle_vector, surveys_balance_vector);
        correlation_trunk_angle_surveys_foot_placement = corr(trunk_angle_vector, surveys_foot_placement_vector);
        correlation_trunk_angle_surveys_walking_speed = corr(trunk_angle_vector, surveys_walking_speed_vector);
        correlation_trunk_angle_surveys_metabolics = corr(trunk_angle_vector, surveys_metabolics_vector);
        correlation_trunk_angle_surveys_difficulty = corr(trunk_angle_vector, surveys_difficulty_vector);
        correlation_trunk_angle_pspeed = Correlation_prompts(trunk_angle_byprompt, "s0", "s1", "s2");
        correlation_trunk_angle_paccuracy = Correlation_prompts(trunk_angle_byprompt, "a0", "a1", "a2");
        correlation_trunk_angle_pbalance = Correlation_prompts(trunk_angle_byprompt, "b0", "b1", "b2");

        correlation_dist_gaze_surveys_balance = corr(dist_gaze_vector, surveys_balance_vector);
        correlation_dist_gaze_surveys_foot_placement = corr(dist_gaze_vector, surveys_foot_placement_vector);
        correlation_dist_gaze_surveys_walking_speed = corr(dist_gaze_vector, surveys_walking_speed_vector);
        correlation_dist_gaze_surveys_metabolics = corr(dist_gaze_vector, surveys_metabolics_vector);
        correlation_dist_gaze_surveys_difficulty = corr(dist_gaze_vector, surveys_difficulty_vector);
        correlation_dist_gaze_pspeed = Correlation_prompts(dist_gaze_byprompt, "s0", "s1", "s2");
        correlation_dist_gaze_paccuracy = Correlation_prompts(dist_gaze_byprompt, "a0", "a1", "a2");
        correlation_dist_gaze_pbalance = Correlation_prompts(dist_gaze_byprompt, "b0", "b1", "b2");
    
        correlation_surveys_balance_surveys_foot_placement = corr(surveys_balance_vector, surveys_foot_placement_vector);
        correlation_surveys_balance_surveys_walking_speed = corr(surveys_balance_vector, surveys_walking_speed_vector);
        correlation_surveys_balance_surveys_metabolics = corr(surveys_balance_vector, surveys_metabolics_vector);
        correlation_surveys_balance_surveys_difficulty = corr(surveys_balance_vector, surveys_difficulty_vector);
        correlation_surveys_balance_pspeed = Correlation_prompts(surveys_balance_byprompt, "s0", "s1", "s2"); %%%%
        correlation_surveys_balance_paccuracy = Correlation_prompts(surveys_balance_byprompt, "a0", "a1", "a2");
        correlation_surveys_balance_pbalance = Correlation_prompts(surveys_balance_byprompt, "b0", "b1", "b2");
    
        correlation_surveys_foot_placement_surveys_walking_speed = corr(surveys_foot_placement_vector, surveys_walking_speed_vector);
        correlation_surveys_foot_placement_surveys_metabolics = corr(surveys_foot_placement_vector, surveys_metabolics_vector);
        correlation_surveys_foot_placement_surveys_difficulty = corr(surveys_foot_placement_vector, surveys_difficulty_vector);
        correlation_surveys_foot_placement_pspeed = Correlation_prompts(surveys_foot_placement_byprompt, "s0", "s1", "s2");
        correlation_surveys_foot_placement_paccuracy = Correlation_prompts(surveys_foot_placement_byprompt, "a0", "a1", "a2");
        correlation_surveys_foot_placement_pbalance = Correlation_prompts(surveys_foot_placement_byprompt, "b0", "b1", "b2");
    
        correlation_surveys_walking_speed_surveys_metabolics = corr(surveys_walking_speed_vector, surveys_metabolics_vector);
        correlation_surveys_walking_speed_surveys_difficulty = corr(surveys_walking_speed_vector, surveys_difficulty_vector);
        correlation_surveys_walking_speed_pspeed = Correlation_prompts(surveys_walking_speed_byprompt, "s0", "s1", "s2");
        correlation_surveys_walking_speed_paccuracy = Correlation_prompts(surveys_walking_speed_byprompt, "a0", "a1", "a2");
        correlation_surveys_walking_speed_pbalance = Correlation_prompts(surveys_walking_speed_byprompt, "b0", "b1", "b2");
    
        correlation_surveys_metabolics_surveys_difficulty = corr(surveys_metabolics_vector, surveys_difficulty_vector);
        correlation_surveys_metabolics_pspeed = Correlation_prompts(surveys_metabolics_byprompt, "s0", "s1", "s2");
        correlation_surveys_metabolics_paccuracy = Correlation_prompts(surveys_metabolics_byprompt, "a0", "a1", "a2");
        correlation_surveys_metabolics_pbalance = Correlation_prompts(surveys_metabolics_byprompt, "b0", "b1", "b2");
    
        correlation_surveys_difficulty_pspeed = Correlation_prompts(surveys_difficulty_byprompt, "s0", "s1", "s2");
        correlation_surveys_difficulty_paccuracy = Correlation_prompts(surveys_difficulty_byprompt, "a0", "a1", "a2");
        correlation_surveys_difficulty_pbalance = Correlation_prompts(surveys_difficulty_byprompt, "b0", "b1", "b2");

        % Plot

        columns_excel = {'C2','D2','E2','F2','G2','H2','I2','J2','K2','L2'}; % where to start plotting the values for each k

        if k >= 1 && k <= 3 %Speed
            
            variables = {'Walking Speed', 'Walking Speed Var', 'Energy Expenditure', 'Step Width', 'Step Width Var', 'Step Length', ...
                        'Step Length Var', 'Stride Duration', 'Stride Duration Var', 'Step Error', 'Step Error Var', 'Head Angle', 'Trunk Angle', 'Look-Ahead Distance', 'Balance Priority', 'Foot Placement Priority', ...
                        'Walking Speed Priority', 'Energy Expenditure Priority', 'Self-Perceived Difficulty', '\it Visual Disturbance\it', '\it Accuracy Prompt\it'};

            variables_x = variables(1:end-2); % Remove 'Speed Prompt', 'Accuracy Prompt', 'Balance Prompt' from the list

            correlation_matrix = [
                1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN; 
                correlation_walkingspeed_walkingspeed_std, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_metabolics, correlation_walkingspeed_std_metabolics, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_step_width, correlation_walkingspeed_std_step_width, correlation_metabolics_step_width, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_step_width_var, correlation_walkingspeed_std_step_width_var, correlation_metabolics_step_width_var, correlation_step_width_step_width_var, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_step_length, correlation_walkingspeed_std_step_length, correlation_metabolics_step_length, correlation_step_width_step_length, correlation_step_width_var_step_length, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_step_length_var, correlation_walkingspeed_std_step_length_var, correlation_metabolics_step_length_var, correlation_step_width_step_length_var, correlation_step_width_var_step_length_var, correlation_step_length_step_length_var, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_stride_time, correlation_walkingspeed_std_stride_time, correlation_metabolics_stride_time, correlation_step_width_stride_time, correlation_step_width_var_stride_time, correlation_step_length_stride_time, correlation_step_length_var_stride_time, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_stride_time_var, correlation_walkingspeed_std_stride_time_var, correlation_metabolics_stride_time_var, correlation_step_width_stride_time_var, correlation_step_width_var_stride_time_var, correlation_step_length_stride_time_var, correlation_step_length_var_stride_time_var, correlation_stride_time_stride_time_var, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_step_error, correlation_walkingspeed_std_step_error, correlation_metabolics_step_error, correlation_step_width_step_error, correlation_step_width_var_step_error, correlation_step_length_step_error, correlation_step_length_var_step_error, correlation_stride_time_step_error, correlation_stride_time_var_step_error, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_step_error_var, correlation_walkingspeed_std_step_error_var, correlation_metabolics_step_error_var, correlation_step_width_step_error_var, correlation_step_width_var_step_error_var, correlation_step_length_step_error_var, correlation_step_length_var_step_error_var, correlation_stride_time_step_error_var, correlation_stride_time_var_step_error_var, correlation_step_error_step_error_var, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_head_angle, correlation_walkingspeed_std_head_angle, correlation_metabolics_head_angle, correlation_step_width_head_angle, correlation_step_width_var_head_angle, correlation_step_length_head_angle, correlation_step_length_var_head_angle, correlation_stride_time_head_angle, correlation_stride_time_var_head_angle, correlation_step_error_head_angle, correlation_step_error_var_head_angle, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_trunk_angle, correlation_walkingspeed_std_trunk_angle, correlation_metabolics_trunk_angle, correlation_step_width_trunk_angle, correlation_step_width_var_trunk_angle, correlation_step_length_trunk_angle, correlation_step_length_var_trunk_angle, correlation_stride_time_trunk_angle, correlation_stride_time_var_trunk_angle, correlation_step_error_trunk_angle, correlation_step_error_var_trunk_angle, correlation_head_angle_trunk_angle, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_dist_gaze, correlation_walkingspeed_std_dist_gaze, correlation_metabolics_dist_gaze, correlation_step_width_dist_gaze, correlation_step_width_var_dist_gaze, correlation_step_length_dist_gaze, correlation_step_length_var_dist_gaze, correlation_stride_time_dist_gaze, correlation_stride_time_var_dist_gaze, correlation_step_error_dist_gaze, correlation_step_error_var_dist_gaze, correlation_head_angle_dist_gaze, correlation_trunk_angle_dist_gaze 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_surveys_balance, correlation_walkingspeed_std_surveys_balance, correlation_metabolics_surveys_balance, correlation_step_width_surveys_balance, correlation_step_width_var_surveys_balance, correlation_step_length_surveys_balance, correlation_step_length_var_surveys_balance, correlation_stride_time_surveys_balance, correlation_stride_time_var_surveys_balance, correlation_step_error_surveys_balance, correlation_step_error_var_surveys_balance, correlation_head_angle_surveys_balance, correlation_trunk_angle_surveys_balance, correlation_dist_gaze_surveys_balance, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_surveys_foot_placement, correlation_walkingspeed_std_surveys_foot_placement, correlation_metabolics_surveys_foot_placement, correlation_step_width_surveys_foot_placement, correlation_step_width_var_surveys_foot_placement, correlation_step_length_surveys_foot_placement, correlation_step_length_var_surveys_foot_placement, correlation_stride_time_surveys_foot_placement, correlation_stride_time_var_surveys_foot_placement, correlation_step_error_surveys_foot_placement, correlation_step_error_var_surveys_foot_placement, correlation_head_angle_surveys_foot_placement, correlation_trunk_angle_surveys_foot_placement, correlation_dist_gaze_surveys_foot_placement, correlation_surveys_balance_surveys_foot_placement, 1, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_surveys_walking_speed, correlation_walkingspeed_std_surveys_walking_speed, correlation_metabolics_surveys_walking_speed, correlation_step_width_surveys_walking_speed, correlation_step_width_var_surveys_walking_speed, correlation_step_length_surveys_walking_speed, correlation_step_length_var_surveys_walking_speed, correlation_stride_time_surveys_walking_speed, correlation_stride_time_var_surveys_walking_speed, correlation_step_error_surveys_walking_speed, correlation_step_error_var_surveys_walking_speed, correlation_head_angle_surveys_walking_speed, correlation_trunk_angle_surveys_walking_speed, correlation_dist_gaze_surveys_walking_speed, correlation_surveys_balance_surveys_walking_speed, correlation_surveys_foot_placement_surveys_walking_speed, 1, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_surveys_metabolics, correlation_walkingspeed_std_surveys_metabolics, correlation_metabolics_surveys_metabolics, correlation_step_width_surveys_metabolics, correlation_step_width_var_surveys_metabolics, correlation_step_length_surveys_metabolics, correlation_step_length_var_surveys_metabolics, correlation_stride_time_surveys_metabolics, correlation_stride_time_var_surveys_metabolics, correlation_step_error_surveys_metabolics, correlation_step_error_var_surveys_metabolics, correlation_head_angle_surveys_metabolics, correlation_trunk_angle_surveys_metabolics, correlation_dist_gaze_surveys_metabolics, correlation_surveys_balance_surveys_metabolics, correlation_surveys_foot_placement_surveys_metabolics, correlation_surveys_walking_speed_surveys_metabolics, 1, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_surveys_difficulty, correlation_walkingspeed_std_surveys_difficulty, correlation_metabolics_surveys_difficulty, correlation_step_width_surveys_difficulty, correlation_step_width_var_surveys_difficulty, correlation_step_length_surveys_difficulty, correlation_step_length_var_surveys_difficulty, correlation_stride_time_surveys_difficulty, correlation_stride_time_var_surveys_difficulty, correlation_step_error_surveys_difficulty, correlation_step_error_var_surveys_difficulty, correlation_head_angle_surveys_difficulty, correlation_trunk_angle_surveys_difficulty, correlation_dist_gaze_surveys_difficulty, correlation_surveys_balance_surveys_difficulty, correlation_surveys_foot_placement_surveys_difficulty, correlation_surveys_walking_speed_surveys_difficulty, correlation_surveys_metabolics_surveys_difficulty, 1, NaN, NaN, NaN;        
                correlation_walkingspeed_pbalance, correlation_walkingspeed_std_pbalance, correlation_metabolics_pbalance, correlation_step_width_pbalance, correlation_step_width_var_pbalance, correlation_step_length_pbalance, correlation_step_length_var_pbalance, correlation_stride_time_pbalance, correlation_stride_time_var_pbalance, correlation_step_error_pbalance, correlation_step_error_var_pbalance, correlation_head_angle_pbalance, correlation_trunk_angle_pbalance, correlation_dist_gaze_pbalance, correlation_surveys_balance_pbalance, correlation_surveys_foot_placement_pbalance, correlation_surveys_walking_speed_pbalance, correlation_surveys_metabolics_pbalance, correlation_surveys_difficulty_pbalance, NaN, NaN, NaN;
                %correlation_walkingspeed_pspeed, correlation_walkingspeed_std_pspeed, correlation_metabolics_pspeed, correlation_step_width_pspeed, correlation_step_width_var_pspeed, correlation_step_length_pspeed, correlation_step_length_var_pspeed, correlation_stride_time_pspeed, correlation_stride_time_var_pspeed, correlation_step_error_pspeed, correlation_step_error_var_pspeed, correlation_head_angle_pspeed, correlation_trunk_angle_pspeed, correlation_dist_gaze_pspeed, correlation_surveys_balance_pspeed, correlation_surveys_foot_placement_pspeed, correlation_surveys_walking_speed_pspeed, correlation_surveys_metabolics_pspeed, correlation_surveys_difficulty_pspeed, NaN, NaN, NaN;
                correlation_walkingspeed_paccuracy, correlation_walkingspeed_std_paccuracy, correlation_metabolics_paccuracy, correlation_step_width_paccuracy, correlation_step_width_var_paccuracy, correlation_step_length_paccuracy, correlation_step_length_var_paccuracy, correlation_stride_time_paccuracy, correlation_stride_time_var_paccuracy, correlation_step_error_paccuracy, correlation_step_error_var_paccuracy, correlation_head_angle_paccuracy, correlation_trunk_angle_paccuracy, correlation_dist_gaze_paccuracy, correlation_surveys_balance_paccuracy, correlation_surveys_foot_placement_paccuracy, correlation_surveys_walking_speed_paccuracy, correlation_surveys_metabolics_paccuracy, correlation_surveys_difficulty_paccuracy, NaN, NaN, NaN;
            ];
            count = 1;
            skipped = 0; % to know if the row was already skipped
            for i=1:size(correlation_matrix,1)
                for j=1:size(correlation_matrix,2)
                    if i == 21 && skipped == 0 % the speed prompt row in the original matrix 
                        all_correlations(count:(count+18)) = NaN; % number of variables compared to the speed prompt in the original matrix - 1 
                        count = count + 19; % number above + 1    
                        skipped = 1;
                        % all_correlations(count) = NaN;
                        % count = count + 1;
                    end
                    if correlation_matrix(i,j) == 1 || isnan(correlation_matrix(i,j))
                        continue
                    end
                    all_correlations(count) = correlation_matrix(i,j);
                    count = count + 1;
                end
            end

            filtered_matrix = correlation_matrix;
            filtered_matrix(filtered_matrix == 1 | isnan(filtered_matrix)) = [];
            flattened_row = filtered_matrix(:)';

            writematrix(all_correlations', filename2, 'Sheet', 1, 'Range', string(columns_excel{k}));

        elseif k >= 4 && k <= 6 % Accuracy

            variables = {'Walking Speed', 'Walking Speed Var', 'Energy Expenditure', 'Step Width', 'Step Width Var', 'Step Length', ...
                        'Step Length Var', 'Stride Duration', 'Stride Duration Var', 'Step Error', 'Step Error Var', 'Head Angle', 'Trunk Angle', 'Look-Ahead Distance', 'Balance Priority', 'Foot Placement Priority', ...
                        'Walking Speed Priority', 'Energy Expenditure Priority', 'Self-Perceived Difficulty', '\it Visual Disturbance\it', '\it Speed Prompt\it'};

            variables_x = variables(1:end-2); % Remove 'Speed Prompt', 'Accuracy Prompt', 'Balance Prompt' from the list

            correlation_matrix = [
                1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN; 
                correlation_walkingspeed_walkingspeed_std, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_metabolics, correlation_walkingspeed_std_metabolics, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_step_width, correlation_walkingspeed_std_step_width, correlation_metabolics_step_width, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_step_width_var, correlation_walkingspeed_std_step_width_var, correlation_metabolics_step_width_var, correlation_step_width_step_width_var, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_step_length, correlation_walkingspeed_std_step_length, correlation_metabolics_step_length, correlation_step_width_step_length, correlation_step_width_var_step_length, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_step_length_var, correlation_walkingspeed_std_step_length_var, correlation_metabolics_step_length_var, correlation_step_width_step_length_var, correlation_step_width_var_step_length_var, correlation_step_length_step_length_var, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_stride_time, correlation_walkingspeed_std_stride_time, correlation_metabolics_stride_time, correlation_step_width_stride_time, correlation_step_width_var_stride_time, correlation_step_length_stride_time, correlation_step_length_var_stride_time, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_stride_time_var, correlation_walkingspeed_std_stride_time_var, correlation_metabolics_stride_time_var, correlation_step_width_stride_time_var, correlation_step_width_var_stride_time_var, correlation_step_length_stride_time_var, correlation_step_length_var_stride_time_var, correlation_stride_time_stride_time_var, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_step_error, correlation_walkingspeed_std_step_error, correlation_metabolics_step_error, correlation_step_width_step_error, correlation_step_width_var_step_error, correlation_step_length_step_error, correlation_step_length_var_step_error, correlation_stride_time_step_error, correlation_stride_time_var_step_error, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_step_error_var, correlation_walkingspeed_std_step_error_var, correlation_metabolics_step_error_var, correlation_step_width_step_error_var, correlation_step_width_var_step_error_var, correlation_step_length_step_error_var, correlation_step_length_var_step_error_var, correlation_stride_time_step_error_var, correlation_stride_time_var_step_error_var, correlation_step_error_step_error_var, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_head_angle, correlation_walkingspeed_std_head_angle, correlation_metabolics_head_angle, correlation_step_width_head_angle, correlation_step_width_var_head_angle, correlation_step_length_head_angle, correlation_step_length_var_head_angle, correlation_stride_time_head_angle, correlation_stride_time_var_head_angle, correlation_step_error_head_angle, correlation_step_error_var_head_angle, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_trunk_angle, correlation_walkingspeed_std_trunk_angle, correlation_metabolics_trunk_angle, correlation_step_width_trunk_angle, correlation_step_width_var_trunk_angle, correlation_step_length_trunk_angle, correlation_step_length_var_trunk_angle, correlation_stride_time_trunk_angle, correlation_stride_time_var_trunk_angle, correlation_step_error_trunk_angle, correlation_step_error_var_trunk_angle, correlation_head_angle_trunk_angle, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_dist_gaze, correlation_walkingspeed_std_dist_gaze, correlation_metabolics_dist_gaze, correlation_step_width_dist_gaze, correlation_step_width_var_dist_gaze, correlation_step_length_dist_gaze, correlation_step_length_var_dist_gaze, correlation_stride_time_dist_gaze, correlation_stride_time_var_dist_gaze, correlation_step_error_dist_gaze, correlation_step_error_var_dist_gaze, correlation_head_angle_dist_gaze, correlation_trunk_angle_dist_gaze 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_surveys_balance, correlation_walkingspeed_std_surveys_balance, correlation_metabolics_surveys_balance, correlation_step_width_surveys_balance, correlation_step_width_var_surveys_balance, correlation_step_length_surveys_balance, correlation_step_length_var_surveys_balance, correlation_stride_time_surveys_balance, correlation_stride_time_var_surveys_balance, correlation_step_error_surveys_balance, correlation_step_error_var_surveys_balance, correlation_head_angle_surveys_balance, correlation_trunk_angle_surveys_balance, correlation_dist_gaze_surveys_balance, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_surveys_foot_placement, correlation_walkingspeed_std_surveys_foot_placement, correlation_metabolics_surveys_foot_placement, correlation_step_width_surveys_foot_placement, correlation_step_width_var_surveys_foot_placement, correlation_step_length_surveys_foot_placement, correlation_step_length_var_surveys_foot_placement, correlation_stride_time_surveys_foot_placement, correlation_stride_time_var_surveys_foot_placement, correlation_step_error_surveys_foot_placement, correlation_step_error_var_surveys_foot_placement, correlation_head_angle_surveys_foot_placement, correlation_trunk_angle_surveys_foot_placement, correlation_dist_gaze_surveys_foot_placement, correlation_surveys_balance_surveys_foot_placement, 1, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_surveys_walking_speed, correlation_walkingspeed_std_surveys_walking_speed, correlation_metabolics_surveys_walking_speed, correlation_step_width_surveys_walking_speed, correlation_step_width_var_surveys_walking_speed, correlation_step_length_surveys_walking_speed, correlation_step_length_var_surveys_walking_speed, correlation_stride_time_surveys_walking_speed, correlation_stride_time_var_surveys_walking_speed, correlation_step_error_surveys_walking_speed, correlation_step_error_var_surveys_walking_speed, correlation_head_angle_surveys_walking_speed, correlation_trunk_angle_surveys_walking_speed, correlation_dist_gaze_surveys_walking_speed, correlation_surveys_balance_surveys_walking_speed, correlation_surveys_foot_placement_surveys_walking_speed, 1, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_surveys_metabolics, correlation_walkingspeed_std_surveys_metabolics, correlation_metabolics_surveys_metabolics, correlation_step_width_surveys_metabolics, correlation_step_width_var_surveys_metabolics, correlation_step_length_surveys_metabolics, correlation_step_length_var_surveys_metabolics, correlation_stride_time_surveys_metabolics, correlation_stride_time_var_surveys_metabolics, correlation_step_error_surveys_metabolics, correlation_step_error_var_surveys_metabolics, correlation_head_angle_surveys_metabolics, correlation_trunk_angle_surveys_metabolics, correlation_dist_gaze_surveys_metabolics, correlation_surveys_balance_surveys_metabolics, correlation_surveys_foot_placement_surveys_metabolics, correlation_surveys_walking_speed_surveys_metabolics, 1, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_surveys_difficulty, correlation_walkingspeed_std_surveys_difficulty, correlation_metabolics_surveys_difficulty, correlation_step_width_surveys_difficulty, correlation_step_width_var_surveys_difficulty, correlation_step_length_surveys_difficulty, correlation_step_length_var_surveys_difficulty, correlation_stride_time_surveys_difficulty, correlation_stride_time_var_surveys_difficulty, correlation_step_error_surveys_difficulty, correlation_step_error_var_surveys_difficulty, correlation_head_angle_surveys_difficulty, correlation_trunk_angle_surveys_difficulty, correlation_dist_gaze_surveys_difficulty, correlation_surveys_balance_surveys_difficulty, correlation_surveys_foot_placement_surveys_difficulty, correlation_surveys_walking_speed_surveys_difficulty, correlation_surveys_metabolics_surveys_difficulty, 1, NaN, NaN, NaN;        
                correlation_walkingspeed_pbalance, correlation_walkingspeed_std_pbalance, correlation_metabolics_pbalance, correlation_step_width_pbalance, correlation_step_width_var_pbalance, correlation_step_length_pbalance, correlation_step_length_var_pbalance, correlation_stride_time_pbalance, correlation_stride_time_var_pbalance, correlation_step_error_pbalance, correlation_step_error_var_pbalance, correlation_head_angle_pbalance, correlation_trunk_angle_pbalance, correlation_dist_gaze_pbalance, correlation_surveys_balance_pbalance, correlation_surveys_foot_placement_pbalance, correlation_surveys_walking_speed_pbalance, correlation_surveys_metabolics_pbalance, correlation_surveys_difficulty_pbalance, NaN, NaN, NaN;
                correlation_walkingspeed_pspeed, correlation_walkingspeed_std_pspeed, correlation_metabolics_pspeed, correlation_step_width_pspeed, correlation_step_width_var_pspeed, correlation_step_length_pspeed, correlation_step_length_var_pspeed, correlation_stride_time_pspeed, correlation_stride_time_var_pspeed, correlation_step_error_pspeed, correlation_step_error_var_pspeed, correlation_head_angle_pspeed, correlation_trunk_angle_pspeed, correlation_dist_gaze_pspeed, correlation_surveys_balance_pspeed, correlation_surveys_foot_placement_pspeed, correlation_surveys_walking_speed_pspeed, correlation_surveys_metabolics_pspeed, correlation_surveys_difficulty_pspeed, NaN, NaN, NaN;
                %correlation_walkingspeed_paccuracy, correlation_walkingspeed_std_paccuracy, correlation_metabolics_paccuracy, correlation_step_width_paccuracy, correlation_step_width_var_paccuracy, correlation_step_length_paccuracy, correlation_step_length_var_paccuracy, correlation_stride_time_paccuracy, correlation_stride_time_var_paccuracy, correlation_step_error_paccuracy, correlation_step_error_var_paccuracy, correlation_head_angle_paccuracy, correlation_trunk_angle_paccuracy, correlation_dist_gaze_paccuracy, correlation_surveys_balance_paccuracy, correlation_surveys_foot_placement_paccuracy, correlation_surveys_walking_speed_paccuracy, correlation_surveys_metabolics_paccuracy, correlation_surveys_difficulty_paccuracy, NaN, NaN, NaN;
            ];

            % Write in the Excel file
            count = 1;
            for i=1:size(correlation_matrix,1)
                for j=1:size(correlation_matrix,2)
                    if i == size(correlation_matrix,1) && j == size(correlation_matrix,2) 
                        all_correlations(count:(count+18)) = NaN; %total number of rows - 3  
                    end
                    if correlation_matrix(i,j) == 1 || isnan(correlation_matrix(i,j))
                        continue
                    end
                    all_correlations(count) = correlation_matrix(i,j);
                    count = count + 1;
                end
            end

            filtered_matrix = correlation_matrix;
            filtered_matrix(filtered_matrix == 1 | isnan(filtered_matrix)) = [];
            flattened_row = filtered_matrix(:)';

            writematrix(all_correlations', filename2, 'Sheet', 1, 'Range', string(columns_excel{k}));

        elseif k >= 7 && k <= 9 % Balance

            variables = {'Walking Speed', 'Walking Speed Var', 'Energy Expenditure', 'Step Width', 'Step Width Var', 'Step Length', ...
                        'Step Length Var', 'Stride Duration', 'Stride Duration Var', 'Step Error', 'Step Error Var', 'Head Angle', 'Trunk Angle', 'Look-Ahead Distance', 'Balance Priority', 'Foot Placement Priority', ...
                        'Walking Speed Priority', 'Energy Expenditure Priority', 'Self-Perceived Difficulty', '\it Speed Prompt\it', '\it Accuracy Prompt\it'};

            variables_x = variables(1:end-2); % Remove 'Speed Prompt', 'Accuracy Prompt', 'Balance Prompt' from the list

            correlation_matrix = [
                1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN; 
                correlation_walkingspeed_walkingspeed_std, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_metabolics, correlation_walkingspeed_std_metabolics, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_step_width, correlation_walkingspeed_std_step_width, correlation_metabolics_step_width, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_step_width_var, correlation_walkingspeed_std_step_width_var, correlation_metabolics_step_width_var, correlation_step_width_step_width_var, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_step_length, correlation_walkingspeed_std_step_length, correlation_metabolics_step_length, correlation_step_width_step_length, correlation_step_width_var_step_length, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_step_length_var, correlation_walkingspeed_std_step_length_var, correlation_metabolics_step_length_var, correlation_step_width_step_length_var, correlation_step_width_var_step_length_var, correlation_step_length_step_length_var, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_stride_time, correlation_walkingspeed_std_stride_time, correlation_metabolics_stride_time, correlation_step_width_stride_time, correlation_step_width_var_stride_time, correlation_step_length_stride_time, correlation_step_length_var_stride_time, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_stride_time_var, correlation_walkingspeed_std_stride_time_var, correlation_metabolics_stride_time_var, correlation_step_width_stride_time_var, correlation_step_width_var_stride_time_var, correlation_step_length_stride_time_var, correlation_step_length_var_stride_time_var, correlation_stride_time_stride_time_var, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_step_error, correlation_walkingspeed_std_step_error, correlation_metabolics_step_error, correlation_step_width_step_error, correlation_step_width_var_step_error, correlation_step_length_step_error, correlation_step_length_var_step_error, correlation_stride_time_step_error, correlation_stride_time_var_step_error, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_step_error_var, correlation_walkingspeed_std_step_error_var, correlation_metabolics_step_error_var, correlation_step_width_step_error_var, correlation_step_width_var_step_error_var, correlation_step_length_step_error_var, correlation_step_length_var_step_error_var, correlation_stride_time_step_error_var, correlation_stride_time_var_step_error_var, correlation_step_error_step_error_var, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_head_angle, correlation_walkingspeed_std_head_angle, correlation_metabolics_head_angle, correlation_step_width_head_angle, correlation_step_width_var_head_angle, correlation_step_length_head_angle, correlation_step_length_var_head_angle, correlation_stride_time_head_angle, correlation_stride_time_var_head_angle, correlation_step_error_head_angle, correlation_step_error_var_head_angle, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_trunk_angle, correlation_walkingspeed_std_trunk_angle, correlation_metabolics_trunk_angle, correlation_step_width_trunk_angle, correlation_step_width_var_trunk_angle, correlation_step_length_trunk_angle, correlation_step_length_var_trunk_angle, correlation_stride_time_trunk_angle, correlation_stride_time_var_trunk_angle, correlation_step_error_trunk_angle, correlation_step_error_var_trunk_angle, correlation_head_angle_trunk_angle, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_dist_gaze, correlation_walkingspeed_std_dist_gaze, correlation_metabolics_dist_gaze, correlation_step_width_dist_gaze, correlation_step_width_var_dist_gaze, correlation_step_length_dist_gaze, correlation_step_length_var_dist_gaze, correlation_stride_time_dist_gaze, correlation_stride_time_var_dist_gaze, correlation_step_error_dist_gaze, correlation_step_error_var_dist_gaze, correlation_head_angle_dist_gaze, correlation_trunk_angle_dist_gaze 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_surveys_balance, correlation_walkingspeed_std_surveys_balance, correlation_metabolics_surveys_balance, correlation_step_width_surveys_balance, correlation_step_width_var_surveys_balance, correlation_step_length_surveys_balance, correlation_step_length_var_surveys_balance, correlation_stride_time_surveys_balance, correlation_stride_time_var_surveys_balance, correlation_step_error_surveys_balance, correlation_step_error_var_surveys_balance, correlation_head_angle_surveys_balance, correlation_trunk_angle_surveys_balance, correlation_dist_gaze_surveys_balance, 1, NaN, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_surveys_foot_placement, correlation_walkingspeed_std_surveys_foot_placement, correlation_metabolics_surveys_foot_placement, correlation_step_width_surveys_foot_placement, correlation_step_width_var_surveys_foot_placement, correlation_step_length_surveys_foot_placement, correlation_step_length_var_surveys_foot_placement, correlation_stride_time_surveys_foot_placement, correlation_stride_time_var_surveys_foot_placement, correlation_step_error_surveys_foot_placement, correlation_step_error_var_surveys_foot_placement, correlation_head_angle_surveys_foot_placement, correlation_trunk_angle_surveys_foot_placement, correlation_dist_gaze_surveys_foot_placement, correlation_surveys_balance_surveys_foot_placement, 1, NaN, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_surveys_walking_speed, correlation_walkingspeed_std_surveys_walking_speed, correlation_metabolics_surveys_walking_speed, correlation_step_width_surveys_walking_speed, correlation_step_width_var_surveys_walking_speed, correlation_step_length_surveys_walking_speed, correlation_step_length_var_surveys_walking_speed, correlation_stride_time_surveys_walking_speed, correlation_stride_time_var_surveys_walking_speed, correlation_step_error_surveys_walking_speed, correlation_step_error_var_surveys_walking_speed, correlation_head_angle_surveys_walking_speed, correlation_trunk_angle_surveys_walking_speed, correlation_dist_gaze_surveys_walking_speed, correlation_surveys_balance_surveys_walking_speed, correlation_surveys_foot_placement_surveys_walking_speed, 1, NaN, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_surveys_metabolics, correlation_walkingspeed_std_surveys_metabolics, correlation_metabolics_surveys_metabolics, correlation_step_width_surveys_metabolics, correlation_step_width_var_surveys_metabolics, correlation_step_length_surveys_metabolics, correlation_step_length_var_surveys_metabolics, correlation_stride_time_surveys_metabolics, correlation_stride_time_var_surveys_metabolics, correlation_step_error_surveys_metabolics, correlation_step_error_var_surveys_metabolics, correlation_head_angle_surveys_metabolics, correlation_trunk_angle_surveys_metabolics, correlation_dist_gaze_surveys_metabolics, correlation_surveys_balance_surveys_metabolics, correlation_surveys_foot_placement_surveys_metabolics, correlation_surveys_walking_speed_surveys_metabolics, 1, NaN, NaN, NaN, NaN;
                correlation_walkingspeed_surveys_difficulty, correlation_walkingspeed_std_surveys_difficulty, correlation_metabolics_surveys_difficulty, correlation_step_width_surveys_difficulty, correlation_step_width_var_surveys_difficulty, correlation_step_length_surveys_difficulty, correlation_step_length_var_surveys_difficulty, correlation_stride_time_surveys_difficulty, correlation_stride_time_var_surveys_difficulty, correlation_step_error_surveys_difficulty, correlation_step_error_var_surveys_difficulty, correlation_head_angle_surveys_difficulty, correlation_trunk_angle_surveys_difficulty, correlation_dist_gaze_surveys_difficulty, correlation_surveys_balance_surveys_difficulty, correlation_surveys_foot_placement_surveys_difficulty, correlation_surveys_walking_speed_surveys_difficulty, correlation_surveys_metabolics_surveys_difficulty, 1, NaN, NaN, NaN;        
                %correlation_walkingspeed_pbalance, correlation_walkingspeed_std_pbalance, correlation_metabolics_pbalance, correlation_step_width_pbalance, correlation_step_width_var_pbalance, correlation_step_length_pbalance, correlation_step_length_var_pbalance, correlation_stride_time_pbalance, correlation_stride_time_var_pbalance, correlation_step_error_pbalance, correlation_step_error_var_pbalance, correlation_head_angle_pbalance, correlation_trunk_angle_pbalance, correlation_dist_gaze_pbalance, correlation_surveys_balance_pbalance, correlation_surveys_foot_placement_pbalance, correlation_surveys_walking_speed_pbalance, correlation_surveys_metabolics_pbalance, correlation_surveys_difficulty_pbalance, NaN, NaN, NaN;
                correlation_walkingspeed_pspeed, correlation_walkingspeed_std_pspeed, correlation_metabolics_pspeed, correlation_step_width_pspeed, correlation_step_width_var_pspeed, correlation_step_length_pspeed, correlation_step_length_var_pspeed, correlation_stride_time_pspeed, correlation_stride_time_var_pspeed, correlation_step_error_pspeed, correlation_step_error_var_pspeed, correlation_head_angle_pspeed, correlation_trunk_angle_pspeed, correlation_dist_gaze_pspeed, correlation_surveys_balance_pspeed, correlation_surveys_foot_placement_pspeed, correlation_surveys_walking_speed_pspeed, correlation_surveys_metabolics_pspeed, correlation_surveys_difficulty_pspeed, NaN, NaN, NaN;
                correlation_walkingspeed_paccuracy, correlation_walkingspeed_std_paccuracy, correlation_metabolics_paccuracy, correlation_step_width_paccuracy, correlation_step_width_var_paccuracy, correlation_step_length_paccuracy, correlation_step_length_var_paccuracy, correlation_stride_time_paccuracy, correlation_stride_time_var_paccuracy, correlation_step_error_paccuracy, correlation_step_error_var_paccuracy, correlation_head_angle_paccuracy, correlation_trunk_angle_paccuracy, correlation_dist_gaze_paccuracy, correlation_surveys_balance_paccuracy, correlation_surveys_foot_placement_paccuracy, correlation_surveys_walking_speed_paccuracy, correlation_surveys_metabolics_paccuracy, correlation_surveys_difficulty_paccuracy, NaN, NaN, NaN;
            ];

            count = 1;
            skipped = 0; % to know if the row was already skipped
            for i=1:size(correlation_matrix,1)
                for j=1:size(correlation_matrix,2)
                    if i == 20 && skipped == 0 % the balance prompt row in the original matrix 
                        all_correlations(count:(count+18)) = NaN; % number of variables compared to the balance prompt in the original matrix  
                        count = count + 19; % number above + 1    
                        skipped = 1;
                        % all_correlations(count) = NaN;
                        % count = count + 1;
                    end
                    if correlation_matrix(i,j) == 1 || isnan(correlation_matrix(i,j))
                        continue
                    end
                    all_correlations(count) = correlation_matrix(i,j);
                    count = count + 1;
                end
            end

            filtered_matrix = correlation_matrix;
            filtered_matrix(filtered_matrix == 1 | isnan(filtered_matrix)) = [];
            flattened_row = filtered_matrix(:)';

            writematrix(all_correlations', filename2, 'Sheet', 1, 'Range', string(columns_excel{k}));
        end

        if show_all_plots == 1
            figure
    
            % Set the colormap with -1 as blue, 0 as white, and 1 as red
            custom_cmap = [linspace(0, 1, 128)', linspace(0, 1, 128)', ones(128, 1); % Blue to white
                           ones(128, 1), linspace(1, 0, 128)', linspace(1, 0, 128)']; % White to red
    
            imagesc(correlation_matrix, 'AlphaData', ~isnan(correlation_matrix)); % Mask NaN values
            colormap(custom_cmap); % Apply custom colormap
            colorbar; % Display the color scale
            clim([-1, 1]); % Set color axis limits
            axis square; % Make the plot square-shaped
    
            set(gca, 'Color', 'white');  % Set the background color for NaN cells to white
    
            xticks(1:length(variables_x));
            yticks(1:length(variables));
            set(gca, 'TickLength', [0.002 0]);
            xticklabels(variables_x);
            yticklabels(variables);
    
            xlim([0.5 18.5])
            set(gca, 'FontSize', 13.5);
            xtickangle(45); % Rotate the x-axis labels for better readability
    
            if k == 1
                title('Linear Correlation Between Gait Variables - Slow Speed');
            elseif k == 2
                title('Linear Correlation Between Gait Variables - Medium Speed');
            elseif k == 3
                title('Linear Correlation Between Gait Variables - Fast Speed');
            elseif k == 4
                title('Linear Correlation Between Gait Variables - Null Accuracy');
            elseif k == 5
                title('Linear Correlation Between Gait Variables - Medium Accuracy');
            elseif k == 6
                title('Linear Correlation Between Gait Variables - High Accuracy');
            elseif k == 7
                title('Linear Correlation Between Gait Variables - No Visual Disturbance');
            elseif k == 8
                title('Linear Correlation Between Gait Variables - Medium Visual Disturbance');
            elseif k == 9
                title('Linear Correlation Between Gait Variables - High Visual Disturbance');
            end
    
            % Add correlation values inside the plot boxes
            % for i = 1:length(variables)
            %     for j = 1:length(variables)
            %         if ~isnan(correlation_matrix(i, j)) % Only display text for non-NaN values
            %             text(j, i, sprintf('%.2f', correlation_matrix(i, j)), ...
            %                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            %                 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'black');
            %         end
            %     end
            % end
            % Add correlation values inside the plot boxes
            for i = 1:length(variables)
                for j = 1:length(variables)
                    if ~isnan(correlation_matrix(i, j)) % Only display text for non-NaN values
                        if option_correlations == 1 % all
                            text(j, i, sprintf('%.2f', correlation_matrix(i, j)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 11, 'FontWeight', 'bold', 'Color', 'black');
                        elseif option_correlations == 2 % close to +1
                            if correlation_matrix(i, j) >= 0.695 && correlation_matrix(i,j) ~= 1 % Modify this to plot only certain values
                                text(j, i, sprintf('%.2f', correlation_matrix(i, j)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 11, 'FontWeight', 'bold', 'Color', 'black');
                            end
                        elseif option_correlations == 3 % close to -1
                            if correlation_matrix(i, j) <= -0.695
                                text(j, i, sprintf('%.2f', correlation_matrix(i, j)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 11, 'FontWeight', 'bold', 'Color', 'black');
                            end
                        elseif option_correlations == 4 % close to 0
                            if correlation_matrix(i, j) >= -0.105 && correlation_matrix(i,j) <= 0.105
                                text(j, i, sprintf('%.2f', correlation_matrix(i, j)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 11, 'FontWeight', 'bold', 'Color', 'black');
                            end
                        end
                    end
                end
            end
        end
        clear correlation_matrix
    end
end



