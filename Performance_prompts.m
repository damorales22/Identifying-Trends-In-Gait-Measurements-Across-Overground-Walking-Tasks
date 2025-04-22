%function [] = Performance_prompts(imported_data,error,colors,subtract_all,subtract_speeds,baseline_errors,bmh)
function [] = Performance_prompts(data_labels,imported_data,baseline_errors,colors,subtract_all,subtract_speeds,participants,shaded_plot,plot_eachparticipant,colors_participants,metabolics,bmh)

    if participants == "multiple"
        bmh_text = "Across all Participants";
    else
        bmh_text = string(bmh);
    end

    bmh_list = string(fieldnames(imported_data));
    file_path = "C:\Users\damor\Documents\DAVID\UNIVERSIDAD\Máster (2023-2025)\2º (2024-2025) - Harvard\Data processing\Data\Survey and Prompts.xlsx";

    point_size = 100; % Point size when plot_eachparticipant==1

    %% Identify headset / no headset in trials 1,2 / 34,35

    headset_beginning = [];
    no_headset_beginning = [];
    headset_end = [];
    no_headset_end = [];

    for i = 1:length(bmh_list)
        bmh = bmh_list{i};
        range = 'N1:N37';
        data_excel = readtable(file_path, 'Sheet', bmh, 'Range', range);

        % Trials 1 and 2
        if data_excel{1,1} == "no headset"
            no_headset_beginning = [no_headset_beginning, 1];
            headset_beginning = [headset_beginning, 2];
        elseif data_excel{2,1} == "no headset"
            no_headset_beginning = [no_headset_beginning, 2];
            headset_beginning = [headset_beginning, 1];
        else
            error("<No headset> text is missing in the Excel at column N at participant" + bmh)
        end

        % Trials 34 and 35
        if data_excel{34,1} == "no headset"
            no_headset_end = [no_headset_end, 34];
            headset_end = [headset_end, 35];
        elseif data_excel{35,1} == "no headset"
            no_headset_end = [no_headset_end, 35];
            headset_end = [headset_end, 34];
        else
            error("<No headset> text is missing in the Excel at column N at participant" + bmh)
        end
    end

    %% Speed input vs actual speed - PLOT 1
    
    slow_speed = [];
    medium_speed = [];
    fast_speed = [];
    for i=1:length(fieldnames(imported_data))
        bmh = bmh_list(i);
        for j=1:size(imported_data.(bmh).prompt.speed,1)
            if j>=3 && j<=29
                if double(imported_data.(bmh).prompt.speed(j,2)) == 0
                    slow_speed = [slow_speed, imported_data.(bmh).walkingspeed(1,j)];
                elseif double(imported_data.(bmh).prompt.speed(j,2)) == 1
                    medium_speed = [medium_speed, imported_data.(bmh).walkingspeed(1,j)];             
                elseif double(imported_data.(bmh).prompt.speed(j,2)) == 2
                    fast_speed = [fast_speed, imported_data.(bmh).walkingspeed(1,j)];
                end
            end
        end
    end
    
    avg_slow_speed = mean(slow_speed);
    avg_medium_speed = mean(medium_speed);
    avg_fast_speed = mean(fast_speed);
    
    std_slow_speed = std(slow_speed);
    std_medium_speed = std(medium_speed);
    std_fast_speed = std(fast_speed);
    
    average_speeds = [avg_slow_speed, avg_medium_speed, avg_fast_speed];
    std_speeds = [std_slow_speed, std_medium_speed, std_fast_speed];
    
    figure;
    b = bar(average_speeds, 'FaceColor', 'flat');
    hold on;
    
    for k = 1:length(average_speeds)
        b.CData(k, :) = colors(k, :);
    end
    
    errorbar(1:3, average_speeds, std_speeds, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    
    set(gca, 'xticklabel', {'Slow', 'Medium', 'Fast'});
    xlabel('Prompt Speed');
    ylabel('Average Walking Speed (m/s)');
    title('Average Walking Speed by Prompt Speed - ' + bmh_text);
    ax = gca;               % Get current axis
    ax.FontSize = 14;       % Set axis number font size (e.g., 12)
    set(gca, 'FontSize', 22);
    hold off;
    
    %% Balance vs straight width and straight width variability - PLOT 2,3
    
    balance_0_step_width = [];
    balance_1_step_width = [];
    balance_2_step_width = [];
    
    balance_0_variability = [];
    balance_1_variability = [];
    balance_2_variability = [];
    
    for i=1:length(fieldnames(imported_data))
        bmh = bmh_list(i);
        for j=1:size(imported_data.(bmh).prompt.balance,1)
            if j>=3 && j<=29
                if double(imported_data.(bmh).prompt.balance(j,2)) == 0
                    balance_0_step_width = [balance_0_step_width, imported_data.(bmh).meanwidthstraights(1,j)];
                    balance_0_variability = [balance_0_variability, imported_data.(bmh).straightsvariability(1,j)];
                elseif double(imported_data.(bmh).prompt.balance(j,2)) == 1
                    balance_1_step_width = [balance_1_step_width, imported_data.(bmh).meanwidthstraights(1,j)];    
                    balance_1_variability = [balance_1_variability, imported_data.(bmh).straightsvariability(1,j)];
                elseif double(imported_data.(bmh).prompt.balance(j,2)) == 2
                    balance_2_step_width = [balance_2_step_width, imported_data.(bmh).meanwidthstraights(1,j)];
                    balance_2_variability = [balance_2_variability, imported_data.(bmh).straightsvariability(1,j)];
                end
            end
        end
    end
    
    % Step width and its std (which is the std across all step widths). 
    avg_balance_0_step_width = mean(balance_0_step_width)/10; %From mm to cm
    avg_balance_1_step_width = mean(balance_1_step_width)/10; %From mm to cm
    avg_balance_2_step_width = mean(balance_2_step_width)/10; %From mm to cm
    std_balance_0_step_width = std(balance_0_step_width)/10; %From mm to cm
    std_balance_1_step_width = std(balance_1_step_width)/10; %From mm to cm
    std_balance_2_step_width = std(balance_2_step_width)/10; %From mm to cm
    
    average_step_width_bybalance = [avg_balance_0_step_width, avg_balance_1_step_width, avg_balance_2_step_width];
    std_step_width_bybalance = [std_balance_0_step_width, std_balance_1_step_width, std_balance_2_step_width];
    
    figure;
    b = bar(average_step_width_bybalance, 'FaceColor', 'flat');
    hold on;
    
    for k = 1:length(average_speeds)
        b.CData(k, :) = colors(k, :);
    end
    
    errorbar(1:3, average_step_width_bybalance, std_step_width_bybalance, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    
    set(gca, 'xticklabel', {'None', 'Medium', 'High'});
    xlabel('Visual Disturbance');
    ylabel('Average Step Width (cm)');
    title('Average Step Width by Balance Condition - ' + bmh_text)
    ax = gca;               % Get current axis
    ax.FontSize = 14;       % Set axis number font size (e.g., 12);
    set(gca, 'FontSize', 22);
    hold off;
    
    % Step width variability, which is the std of each single trial
    avg_balance_0_variability = mean(balance_0_variability)/10; %From mm to cm
    avg_balance_1_variability = mean(balance_1_variability)/10; %From mm to cm
    avg_balance_2_variability = mean(balance_2_variability)/10; %From mm to cm
    std_balance_0_variability = std(balance_0_variability)/10; %From mm to cm
    std_balance_1_variability = std(balance_1_variability)/10; %From mm to cm
    std_balance_2_variability = std(balance_2_variability)/10; %From mm to cm
    
    average_step_width_variability_bybalance = [avg_balance_0_variability, avg_balance_1_variability, avg_balance_2_variability];
    std_step_width_variability_bybalance = [std_balance_0_variability, std_balance_1_variability, std_balance_2_variability];
    
    figure;
    b = bar(average_step_width_variability_bybalance, 'FaceColor', 'flat');
    hold on;
    
    for k = 1:length(average_speeds)
        b.CData(k, :) = colors(k, :);
    end
    
    errorbar(1:3, average_step_width_variability_bybalance, std_step_width_variability_bybalance, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    
    set(gca, 'xticklabel', {'None', 'Medium', 'High'});
    xlabel('Visual Disturbance');
    ylabel('Average Step Width Variability (cm)');
    title('Average Step Width Variability by Balance Condition - ' + bmh_text);
    ax = gca;               % Get current axis
    ax.FontSize = 14;       % Set axis number font size (e.g., 12)
    set(gca, 'FontSize', 22);
    hold off;
           
    %% First vs last trial with headset vs no headset  - PLOT 4,5
    %Compare step width, step width variability, step length
    %Trial 1 and 35 (variable): Walk at your typical speed and ignore the red lines on the ground (no headset)
    %Trial 2 and 34 (variable): Walk at your typical speed and ignore the red lines on the ground (headset)
    
    walking_speed = [];
    if participants == "single"
        walking_speed = imported_data.(bmh).walkingspeed;
    elseif participants == "multiple"
        for i=1:size(fieldnames(imported_data),1)
            names = string(fieldnames(imported_data));
            bmh = names(i);
            walking_speed(i,:) = imported_data.(bmh).walkingspeed;
        end
   end

    % First set of data - Without the headset
    step_width1 = [];
    step_width35 = [];
    step_width_var1 = [];
    step_width_var35 = [];
    step_length1 = [];
    step_length35 = [];
    step_length_var1 = [];
    step_length_var35 = [];
    
    step_width1_curves = [];
    step_width35_curves = [];
    step_width_var1_curves = [];
    step_width_var35_curves = [];
    step_length1_curves = [];
    step_length35_curves = [];
    step_length_var1_curves = [];
    step_length_var35_curves = [];

    metabolics1 = [];
    metabolics35 = [];

    walking_speed1 = [];
    walking_speed35 = [];
    
    % No headset
    for i = 1:size(bmh_list,1)
        bmh = bmh_list(i);
    
        step_width1 = [step_width1, imported_data.(bmh).meanwidthstraights(1,no_headset_beginning(i))];
        step_width35 = [step_width35, imported_data.(bmh).meanwidthstraights(1,no_headset_end(i))];
        step_width_var1 = [step_width_var1, imported_data.(bmh).straightsvariability(1,no_headset_beginning(i))];
        step_width_var35 = [step_width_var35, imported_data.(bmh).straightsvariability(1,no_headset_end(i))];
        step_length1 = [step_length1, imported_data.(bmh).meanlength_straights(1,no_headset_beginning(i))];
        step_length35 = [step_length35, imported_data.(bmh).meanlength_straights(1,no_headset_end(i))];
        step_length_var1 = [step_length_var1, imported_data.(bmh).meanlength_straights_variability(1,no_headset_beginning(i))];
        step_length_var35 = [step_length_var35, imported_data.(bmh).meanlength_straights_variability(1,no_headset_end(i))];

        step_width1_curves = [step_width1_curves, imported_data.(bmh).meanwidthcurves(1,no_headset_beginning(i))];
        step_width35_curves = [step_width35_curves, imported_data.(bmh).meanwidthcurves(1,no_headset_end(i))];
        step_width_var1_curves = [step_width_var1_curves, imported_data.(bmh).curvesvariability(1,no_headset_beginning(i))];
        step_width_var35_curves = [step_width_var35_curves, imported_data.(bmh).curvesvariability(1,no_headset_end(i))];
        step_length1_curves = [step_length1_curves, imported_data.(bmh).meanlength_curves(1,no_headset_beginning(i))];
        step_length35_curves = [step_length35_curves, imported_data.(bmh).meanlength_curves(1,no_headset_end(i))];
        step_length_var1_curves = [step_length_var1_curves, imported_data.(bmh).meanlength_curves_variability(1,no_headset_beginning(i))];
        step_length_var35_curves = [step_length_var35_curves, imported_data.(bmh).meanlength_curves_variability(1,no_headset_end(i))];

        metabolics1 = [metabolics1, table2array(metabolics.(bmh)(1,no_headset_beginning(i)))];
        metabolics35 = [metabolics35, table2array(metabolics.(bmh)(1,(no_headset_end(i)-4)))];

        walking_speed1 = [walking_speed1, walking_speed(i,no_headset_beginning(i))];
        walking_speed35 = [walking_speed35, walking_speed(i,no_headset_end(i))];
    end
    
    % Straights
    step_width1_avg = mean(step_width1)/10; %From mm to cm
    step_width35_avg = mean(step_width35)/10; %From mm to cm
    step_width_var1_avg = mean(step_width_var1)/10; %From mm to cm
    step_width_var35_avg = mean(step_width_var35)/10; %From mm to cm
    step_length1_avg = mean(step_length1)/10; %From mm to cm
    step_length35_avg = mean(step_length35)/10; %From mm to cm
    step_length_var1_avg = mean(step_length_var1)/10; %From mm to cm
    step_length_var35_avg = mean(step_length_var35)/10; %From mm to cm

    step_width1_std = std(step_width1)/10; %From mm to cm
    step_width35_std = std(step_width35)/10; %From mm to cm
    step_width_var1_std = std(step_width_var1)/10; %From mm to cm
    step_width_var35_std = std(step_width_var35)/10; %From mm to cm
    step_length1_std = std(step_length1)/10; %From mm to cm
    step_length35_std = std(step_length35)/10; %From mm to cm
    step_length_var1_std = std(step_length_var1)/10; %From mm to cm
    step_length_var35_std = std(step_length_var35)/10; %From mm to cm

    metabolics1_avg = mean(metabolics1);
    metabolics35_avg = mean(metabolics35);

    % Curves
    step_width1_avg_curves = mean(step_width1_curves)/10; %From mm to cm
    step_width35_avg_curves = mean(step_width35_curves)/10; %From mm to cm
    step_width_var1_avg_curves = mean(step_width_var1_curves)/10; %From mm to cm
    step_width_var35_avg_curves = mean(step_width_var35_curves)/10; %From mm to cm
    step_length1_avg_curves = mean(step_length1_curves)/10; %From mm to cm
    step_length35_avg_curves = mean(step_length35_curves)/10; %From mm to cm
    step_length_var1_avg_curves = mean(step_length_var1_curves)/10; %From mm to cm
    step_length_var35_avg_curves = mean(step_length_var35_curves)/10; %From mm to cm

    step_width1_std_curves = std(step_width1_curves)/10; %From mm to cm
    step_width35_std_curves = std(step_width35_curves)/10; %From mm to cm
    step_width_var1_std_curves = std(step_width_var1_curves)/10; %From mm to cm
    step_width_var35_std_curves = std(step_width_var35_curves)/10; %From mm to cm
    step_length1_std_curves = std(step_length1_curves)/10; %From mm to cm
    step_length35_std_curves = std(step_length35_curves)/10; %From mm to cm
    step_length_var1_std_curves = std(step_length_var1_curves)/10; %From mm to cm
    step_length_var35_std_curves = std(step_length_var35_curves)/10; %From mm to cm

    data_without = [step_width1_avg, step_width35_avg; step_width_var1_avg, step_width_var35_avg; step_width1_avg_curves, step_width35_avg_curves; step_width_var1_avg_curves, step_width_var35_avg_curves;
                    step_length1_avg, step_length35_avg; step_length_var1_avg, step_length_var35_avg; step_length1_avg_curves, step_length35_avg_curves; step_length_var1_avg_curves, step_length_var35_avg_curves];
    data_without_std = [step_width1_std, step_width35_std; step_width_var1_std, step_width_var35_std; step_width1_std_curves, step_width35_std_curves; step_width_var1_std_curves, step_width_var35_std_curves;
                        step_length1_std, step_length35_std; step_length_var1_std, step_length_var35_std; step_length1_std_curves, step_length35_std_curves; step_length_var1_std_curves, step_length_var35_std_curves];
        
    % Second set of data - With the headset
    step_width2 = [];
    step_width34 = [];
    step_width_var2 = [];
    step_width_var34 = [];
    step_length2 = [];
    step_length34 = [];
    step_length_var2 = [];
    step_length_var34 = [];

    step_width2_curves = [];
    step_width34_curves = [];
    step_width_var2_curves = [];
    step_width_var34_curves = [];
    step_length2_curves = [];
    step_length34_curves = [];
    step_length_var2_curves = [];
    step_length_var34_curves = [];

    metabolics2 = [];
    metabolics34 = [];

    walking_speed2 = [];
    walking_speed34 = [];
    
    % Headset
    for i = 1:length(fieldnames(imported_data))
        bmh = bmh_list(i);
    
        step_width2 = [step_width2, imported_data.(bmh).meanwidthstraights(1,headset_beginning(i))];
        step_width34 = [step_width34, imported_data.(bmh).meanwidthstraights(1,headset_end(i))];
        step_width_var2 = [step_width_var2, imported_data.(bmh).straightsvariability(1,headset_beginning(i))];
        step_width_var34 = [step_width_var34, imported_data.(bmh).straightsvariability(1,headset_end(i))];
        step_length2 = [step_length2, imported_data.(bmh).meanlength_straights(1,headset_beginning(i))];
        step_length34 = [step_length34, imported_data.(bmh).meanlength_straights(1,headset_end(i))];
        step_length_var2 = [step_length_var2, imported_data.(bmh).meanlength_straights_variability(1,headset_beginning(i))];
        step_length_var34 = [step_length_var34, imported_data.(bmh).meanlength_straights_variability(1,headset_end(i))];

        step_width2_curves = [step_width2_curves, imported_data.(bmh).meanwidthcurves(1,headset_beginning(i))];
        step_width34_curves = [step_width34_curves, imported_data.(bmh).meanwidthcurves(1,headset_end(i))];
        step_width_var2_curves = [step_width_var2_curves, imported_data.(bmh).curvesvariability(1,headset_beginning(i))];
        step_width_var34_curves = [step_width_var34_curves, imported_data.(bmh).curvesvariability(1,headset_end(i))];
        step_length2_curves = [step_length2_curves, imported_data.(bmh).meanlength_curves(1,headset_beginning(i))];
        step_length34_curves = [step_length34_curves, imported_data.(bmh).meanlength_curves(1,headset_end(i))];
        step_length_var2_curves = [step_length_var2_curves, imported_data.(bmh).meanlength_curves_variability(1,headset_beginning(i))];
        step_length_var34_curves = [step_length_var34_curves, imported_data.(bmh).meanlength_curves_variability(1,headset_end(i))];

        metabolics2 = [metabolics2, table2array(metabolics.(bmh)(1,headset_beginning(i)))];
        metabolics34 = [metabolics34, table2array(metabolics.(bmh)(1,(headset_end(i)-4)))];

        walking_speed2 = [walking_speed2, walking_speed(i,headset_beginning(i))];
        walking_speed34 = [walking_speed34, walking_speed(i,headset_end(i))];
    end
    
    % Straights
    step_width2_avg = mean(step_width2)/10; %From mm to cm
    step_width34_avg = mean(step_width34)/10; %From mm to cm
    step_width_var2_avg = mean(step_width_var2)/10; %From mm to cm
    step_width_var34_avg = mean(step_width_var34)/10; %From mm to cm
    step_length2_avg = mean(step_length2)/10; %From mm to cm
    step_length34_avg = mean(step_length34)/10; %From mm to cm
    step_length_var2_avg = mean(step_length_var2)/10; %From mm to cm
    step_length_var34_avg = mean(step_length_var34)/10; %From mm to cm

    step_width2_std = std(step_width2)/10; %From mm to cm
    step_width34_std = std(step_width34)/10; %From mm to cm
    step_width_var2_std = std(step_width_var2)/10; %From mm to cm
    step_width_var34_std = std(step_width_var34)/10; %From mm to cm
    step_length2_std = std(step_length2)/10; %From mm to cm
    step_length34_std = std(step_length34)/10; %From mm to cm
    step_length_var2_std = std(step_length_var2)/10; %From mm to cm
    step_length_var34_std = std(step_length_var34)/10; %From mm to cm

    metabolics2_avg = mean(metabolics2);
    metabolics34_avg = mean(metabolics34);

    % Curves
    step_width2_avg_curves = mean(step_width2_curves)/10; %From mm to cm
    step_width34_avg_curves = mean(step_width34_curves)/10; %From mm to cm
    step_width_var2_avg_curves = mean(step_width_var2_curves)/10; %From mm to cm
    step_width_var34_avg_curves = mean(step_width_var34_curves)/10; %From mm to cm
    step_length2_avg_curves = mean(step_length2_curves)/10; %From mm to cm
    step_length34_avg_curves = mean(step_length34_curves)/10; %From mm to cm
    step_length_var2_avg_curves = mean(step_length_var2_curves)/10; %From mm to cm
    step_length_var34_avg_curves = mean(step_length_var34_curves)/10; %From mm to cm

    step_width2_std_curves = std(step_width2_curves)/10; %From mm to cm
    step_width34_std_curves = std(step_width34_curves)/10; %From mm to cm
    step_width_var2_std_curves = std(step_width_var2_curves)/10; %From mm to cm
    step_width_var34_std_curves = std(step_width_var34_curves)/10; %From mm to cm
    step_length2_std_curves = std(step_length2_curves)/10; %From mm to cm
    step_length34_std_curves = std(step_length34_curves)/10; %From mm to cm
    step_length_var2_std_curves = std(step_length_var2_curves)/10; %From mm to cm
    step_length_var34_std_curves = std(step_length_var34_curves)/10; %From mm to cm

    data_with = [step_width2_avg, step_width34_avg; step_width_var2_avg, step_width_var34_avg; step_width2_avg_curves, step_width34_avg_curves; step_width_var2_avg_curves, step_width_var34_avg_curves;
                 step_length2_avg, step_length34_avg; step_length_var2_avg, step_length_var34_avg; step_length2_avg_curves, step_length34_avg_curves; step_length_var2_avg_curves, step_length_var34_avg_curves];
    data_with_std = [step_width2_std, step_width34_std; step_width_var2_std, step_width_var34_std; step_width2_std_curves, step_width34_std_curves; step_width_var2_std_curves, step_width_var34_std_curves;
                     step_length2_std, step_length34_std; step_length_var2_std, step_length_var34_std; step_length2_std_curves, step_length34_std_curves; step_length_var2_std_curves, step_length_var34_std_curves];

    % Walking speed
    walking_speed1_avg = mean(walking_speed1);
    walking_speed1_std = std(walking_speed1);
    walking_speed35_avg = mean(walking_speed35);
    walking_speed35_std = std(walking_speed35);
    walking_speed2_avg = mean(walking_speed2);
    walking_speed2_std = std(walking_speed2);
    walking_speed34_avg = mean(walking_speed34);
    walking_speed34_std = std(walking_speed34);

    disp('Walking Speed no Headset (at the Beginning): ' + string(walking_speed1_avg) + ' +- ' + string(walking_speed1_std) + ' m/s');
    disp('Walking Speed no Headset (at the End): ' + string(walking_speed35_avg) + ' +- ' + string(walking_speed35_std) + ' m/s');
    disp('Walking Speed Headset (at the Beginning): ' + string(walking_speed2_avg) + ' +- ' + string(walking_speed2_std) + ' m/s');
    disp('Walking Speed Headset (at the End): ' + string(walking_speed34_avg) + ' +- ' + string(walking_speed34_std) + ' m/s');
    
    disp('Walking Speed no Headset: ' + string((walking_speed1_avg+walking_speed35_avg)/2) + ' +- ' + string((walking_speed1_std+walking_speed35_std)/2) + ' m/s');
    disp('Walking Speed Headset: ' + string((walking_speed2_avg+walking_speed34_avg)/2) + ' +- ' + string((walking_speed2_std+walking_speed34_std)/2) + ' m/s');

    disp('EE before and after trials without headset (respectively): ' + string(metabolics1_avg) + '; ' + string(metabolics35_avg) + ' W');
    disp('EE before and after trials with headset (respectively): ' + string(metabolics2_avg) + '; ' + string(metabolics34_avg) + ' W');
    disp('EE no headset vs headset (respectively): ' + string((metabolics1_avg+metabolics35_avg)/2) + '; ' + string((metabolics2_avg+metabolics34_avg)/2) + ' W');

    disp('EE before and after trials without headset (respectively): ' + string(metabolics1_avg/70) + '; ' + string(metabolics35_avg/70) + ' W/kg');
    disp('EE before and after trials with headset (respectively): ' + string(metabolics2_avg/70) + '; ' + string(metabolics34_avg/70) + ' W/kg');
    disp('EE no headset vs headset (respectively): ' + string((metabolics1_avg+metabolics35_avg)/(2*70)) + '; ' + string((metabolics2_avg+metabolics34_avg)/(2*70)) + ' W/kg');

    figure;
    ax = gca;               % Get current axis
    ax.FontSize = 14;       % Set axis number font size (e.g., 12)
    set(gca, 'FontSize', 22);

    groupLabels_left = {'SL Straights','SL Curves'};
    groupLabels_right = {
        'SLV Straights', ...
        'SLV Curves', ...
        'SW Straights', ...
        'SWV Straights', ...
        'SW Curves', ...
        'SWV Curves'        
    };
    conditionLabels = {'Before Trials', 'After Trials'};

    data_without_left(1,:) = data_without(5,:);
    data_without_left(2,:) = data_without(7,:);
    data_without_std_left(1,:) = data_without_std(5,:);
    data_without_std_left(2,:) = data_without_std(7,:);

    data_without_right(1,:) = data_without(6,:);
    data_without_right(2,:) = data_without(8,:);
    data_without_right(3:6,:) = data_without(1:4,:);
    data_without_std_right(1,:) = data_without_std(6,:);
    data_without_std_right(2,:) = data_without_std(8,:);
    data_without_std_right(3:6,:) = data_without_std(1:4,:);
  
    data_with_left(1,:) = data_with(5,:);
    data_with_left(2,:) = data_with(7,:);
    data_with_std_left(1,:) = data_with_std(5,:);
    data_with_std_left(2,:) = data_with_std(7,:);

    data_with_right(1,:) = data_with(6,:);
    data_with_right(2,:) = data_with(8,:);
    data_with_right(3:6,:) = data_with(1:4,:);
    data_with_std_right(1,:) = data_with_std(6,:);
    data_with_std_right(2,:) = data_with_std(8,:);
    data_with_std_right(3:6,:) = data_with_std(1:4,:);
    
    % Without the headset
    subplot('Position', [0.05, 0.2, 0.28, 0.7]); % Left plot
    hold on;
    
    b1 = bar(data_without_left, 'grouped');
    b1(1).FaceColor = [0.2, 0.6, 0.8];
    b1(2).FaceColor = [0.9, 0.4, 0.4];
    
    x_positions = (1:size(data_without_left, 1))' - 0.15;
    
    for i = 1:size(data_without_left, 2)
        errorbar(x_positions + (i - 1) * 0.3, data_without_left(:, i), data_without_std_left(:, i), ...
                 'k', 'LineStyle', 'none', 'CapSize', 5, 'LineWidth', 1);
    end
    
    set(gca, 'FontSize', 22);
    set(gca, 'XTick', 1:size(data_without, 1));
    set(gca, 'XTickLabel', groupLabels_left);
    %legend(conditionLabels);
    if participants == "single"
        title('Without the Headset - ' + bmh_text);
    else
        title('Step Length Without the Headset')
    end
    
    ylabel('Mean Value (cm)');
    grid on;
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.YMinorGrid = 'off';
    xticks(1:size(data_without_left, 1));
    yticks(0:10:max(ylim));
    ylim([0 82])
    
    subplot('Position', [0.38, 0.2, 0.6, 0.7]); % Right plot
    hold on;
    
    b1 = bar(data_without_right, 'grouped');
    b1(1).FaceColor = [0.2, 0.6, 0.8];
    b1(2).FaceColor = [0.9, 0.4, 0.4];
    
    x_positions = (1:size(data_without_right, 1))' - 0.15;
    
    for i = 1:size(data_without_right, 2)
        errorbar(x_positions + (i - 1) * 0.3, data_without_right(:, i), data_without_std_right(:, i), ...
                 'k', 'LineStyle', 'none', 'CapSize', 5, 'LineWidth', 1);
    end
    
    set(gca, 'FontSize', 22);
    set(gca, 'XTick', 1:size(data_without, 1));
    set(gca, 'XTickLabel', groupLabels_right);
    legend(conditionLabels);
    if participants == "single"
        title('Without the Headset - ' + bmh_text);
    else
        title('Various Parameters Without the Headset')
    end
    ylabel('Mean Value (cm)');
    yticks(0:50:max(ylim));
    grid on;
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.YMinorGrid = 'off';
    xticks(1:size(data_without_right, 1));
    yticks(0:2:23);
    ylim([0,23])
    
    % With the headset
    figure;
    subplot('Position', [0.05, 0.2, 0.28, 0.7]); % Left plot
    hold on;
    
    b1 = bar(data_with_left, 'grouped');
    b1(1).FaceColor = [0.2, 0.6, 0.8];
    b1(2).FaceColor = [0.9, 0.4, 0.4];
    
    x_positions = (1:size(data_with_left, 1))' - 0.15;
    
    for i = 1:size(data_with_left, 2)
        errorbar(x_positions + (i - 1) * 0.3, data_with_left(:, i), data_with_std_left(:, i), ...
                 'k', 'LineStyle', 'none', 'CapSize', 5, 'LineWidth', 1);
    end
    
    set(gca, 'FontSize', 22);
    set(gca, 'XTick', 1:size(data_with, 1));
    set(gca, 'XTickLabel', groupLabels_left);
    %legend(conditionLabels);
    if participants == "single"
        title('With the Headset - ' + bmh_text);
    else
        title('Step Length With the Headset')
    end
    ylabel('Mean Value (cm)');
    grid on;
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.YMinorGrid = 'off';
    xticks(1:size(data_with_left, 1));
    yticks(0:10:82);
    ylim([0 82])
    
    subplot('Position', [0.38, 0.2, 0.6, 0.7]); % Right plot
    hold on;
    
    b1 = bar(data_with_right, 'grouped');
    b1(1).FaceColor = [0.2, 0.6, 0.8];
    b1(2).FaceColor = [0.9, 0.4, 0.4];
    
    x_positions = (1:size(data_with_right, 1))' - 0.15;
    
    for i = 1:size(data_with_right, 2)
        errorbar(x_positions + (i - 1) * 0.3, data_with_right(:, i), data_with_std_right(:, i), ...
                 'k', 'LineStyle', 'none', 'CapSize', 5, 'LineWidth', 1);
    end
    
    set(gca, 'FontSize', 22);
    set(gca, 'XTick', 1:size(data_with, 1));
    set(gca, 'XTickLabel', groupLabels_right);
    legend(conditionLabels);
    if participants == "single"
        title('With the Headset - ' + bmh_text);
    else
        title('Various Parameters With the Headset')
    end
    ylabel('Mean Value (cm)');
    yticks(0:50:max(ylim));
    grid on;
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.YMinorGrid = 'off';
    xticks(1:size(data_with_right, 1));
    yticks(0:2:22);
    ylim([0,23])
    
    %% Balance vs step length and step length variability - PLOT 6,7

    balance_0_step_length = [];
    balance_1_step_length = [];
    balance_2_step_length = [];
    
    balance_0_variability = [];
    balance_1_variability = [];
    balance_2_variability = [];
    
    for i=1:length(fieldnames(imported_data))
        bmh = bmh_list(i);
        for j=1:size(imported_data.(bmh).prompt.balance,1)
            if j>=3 && j<=29
                if double(imported_data.(bmh).prompt.balance(j,2)) == 0
                    balance_0_step_length = [balance_0_step_length, imported_data.(bmh).meanlength_straights(1,j)];
                    balance_0_variability = [balance_0_variability, imported_data.(bmh).meanlength_straights_variability(1,j)];
                elseif double(imported_data.(bmh).prompt.balance(j,2)) == 1
                    balance_1_step_length = [balance_1_step_length, imported_data.(bmh).meanlength_straights(1,j)];    
                    balance_1_variability = [balance_1_variability, imported_data.(bmh).meanlength_straights_variability(1,j)];
                elseif double(imported_data.(bmh).prompt.balance(j,2)) == 2
                    balance_2_step_length = [balance_2_step_length, imported_data.(bmh).meanlength_straights(1,j)];
                    balance_2_variability = [balance_2_variability, imported_data.(bmh).meanlength_straights_variability(1,j)];
                end
            end
        end
    end
    
    % Step length and its std (which is the std across all step lengths). 
    avg_balance_0_step_length = mean(balance_0_step_length)/10; %From mm to cm
    avg_balance_1_step_length = mean(balance_1_step_length)/10; %From mm to cm
    avg_balance_2_step_length = mean(balance_2_step_length)/10; %From mm to cm
    % std_balance_0_step_length = std(balance_0_step_length)/10; %From mm to cm
    % std_balance_1_step_length = std(balance_1_step_length)/10; %From mm to cm
    % std_balance_2_step_length = std(balance_2_step_length)/10; %From mm to cm
    std_balance_0_step_length = std(mean(reshape(balance_0_step_length, 9, []), 2), 'omitnan')/10; %do first the mean of the stds each 9 positions, which correspond to the same condition for all participants
    std_balance_1_step_length = std(mean(reshape(balance_1_step_length, 9, []), 2), 'omitnan')/10;
    std_balance_2_step_length = std(mean(reshape(balance_2_step_length, 9, []), 2), 'omitnan')/10;
    
    average_step_length_bybalance = [avg_balance_0_step_length, avg_balance_1_step_length, avg_balance_2_step_length];
    std_step_length_bybalance = [std_balance_0_step_length, std_balance_1_step_length, std_balance_2_step_length];
    
    figure;
    b = bar(average_step_length_bybalance, 'FaceColor', 'flat');
    hold on;
    
    for k = 1:length(average_speeds)
        b.CData(k, :) = colors(k, :);
    end
    
    errorbar(1:3, average_step_length_bybalance, std_step_length_bybalance, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    
    set(gca, 'xticklabel', {'None', 'Medium', 'High'});
    xlabel('Visual Disturbance');
    ylabel('Average Step Length (cm)');
    title('Average Step Length by Balance Condition - ' + bmh_text)
    ax = gca;               % Get current axis
    ax.FontSize = 14;       % Set axis number font size (e.g., 12);
    set(gca, 'FontSize', 22);
    hold off;
    
    % Step width variability, which is the std of each single trial
    avg_balance_0_variability = mean(balance_0_variability)/10; %From mm to cm
    avg_balance_1_variability = mean(balance_1_variability)/10; %From mm to cm
    avg_balance_2_variability = mean(balance_2_variability)/10; %From mm to cm
    % std_balance_0_variability = std(balance_0_variability)/10; %From mm to cm
    % std_balance_1_variability = std(balance_1_variability)/10; %From mm to cm
    % std_balance_2_variability = std(balance_2_variability)/10; %From mm to cm
    std_balance_0_variability = std(mean(reshape(balance_0_variability, 9, []), 2), 'omitnan')/10; %do first the mean of the stds each 9 positions, which correspond to the same condition for all participants
    std_balance_1_variability = std(mean(reshape(balance_1_variability, 9, []), 2), 'omitnan')/10;
    std_balance_2_variability = std(mean(reshape(balance_2_variability, 9, []), 2), 'omitnan')/10;
    
    average_step_length_variability_bybalance = [avg_balance_0_variability, avg_balance_1_variability, avg_balance_2_variability];
    std_step_length_variability_bybalance = [std_balance_0_variability, std_balance_1_variability, std_balance_2_variability];
    
    figure;
    b = bar(average_step_length_variability_bybalance, 'FaceColor', 'flat');
    hold on;
    
    for k = 1:length(average_speeds)
        b.CData(k, :) = colors(k, :);
    end
    
    errorbar(1:3, average_step_length_variability_bybalance, std_step_length_variability_bybalance, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    
    set(gca, 'xticklabel', {'None', 'Medium', 'High'});
    xlabel('Visual Disturbance');
    ylabel('Average Step Length Variability (cm)');
    title('Average Step Length Variability by Balance Condition - ' + bmh_text);
    ax = gca;               % Get current axis
    ax.FontSize = 14;       % Set axis number font size (e.g., 12)
    set(gca, 'FontSize', 22);
    hold off;

    %% Accuracy vs step length and step length variability - PLOT 8,9

    accuracy_0_step_length = [];
    accuracy_1_step_length = [];
    accuracy_2_step_length = [];
    
    accuracy_0_variability = [];
    accuracy_1_variability = [];
    accuracy_2_variability = [];
    
    for i=1:length(fieldnames(imported_data))
        bmh = bmh_list(i);
        for j=1:size(imported_data.(bmh).prompt.accuracy,1)
            if j>=3 && j<=29
                if double(imported_data.(bmh).prompt.accuracy(j,2)) == 0
                    accuracy_0_step_length = [accuracy_0_step_length, imported_data.(bmh).meanlength_straights(1,j)];
                    accuracy_0_variability = [accuracy_0_variability, imported_data.(bmh).meanlength_straights_variability(1,j)];
                elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 1
                    accuracy_1_step_length = [accuracy_1_step_length, imported_data.(bmh).meanlength_straights(1,j)];    
                    accuracy_1_variability = [accuracy_1_variability, imported_data.(bmh).meanlength_straights_variability(1,j)];
                elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 2
                    accuracy_2_step_length = [accuracy_2_step_length, imported_data.(bmh).meanlength_straights(1,j)];
                    accuracy_2_variability = [accuracy_2_variability, imported_data.(bmh).meanlength_straights_variability(1,j)];
                end
            end
        end
    end
    
    % Step length and its std (which is the std across all step lengths). 
    avg_accuracy_0_step_length = mean(accuracy_0_step_length)/10; %From mm to cm
    avg_accuracy_1_step_length = mean(accuracy_1_step_length)/10; %From mm to cm
    avg_accuracy_2_step_length = mean(accuracy_2_step_length)/10; %From mm to cm
    % std_accuracy_0_step_length = std(accuracy_0_step_length)/10; %From mm to cm
    % std_accuracy_1_step_length = std(accuracy_1_step_length)/10; %From mm to cm
    % std_accuracy_2_step_length = std(accuracy_2_step_length)/10; %From mm to cm
    std_accuracy_0_step_length = std(mean(reshape(accuracy_0_step_length, 9, []), 2), 'omitnan')/10; %do first the mean of the stds each 9 positions, which correspond to the same condition for all participants
    std_accuracy_1_step_length = std(mean(reshape(accuracy_1_step_length, 9, []), 2), 'omitnan')/10;
    std_accuracy_2_step_length = std(mean(reshape(accuracy_2_step_length, 9, []), 2), 'omitnan')/10;
    
    average_step_length_byaccuracy = [avg_accuracy_0_step_length, avg_accuracy_1_step_length, avg_accuracy_2_step_length];
    std_step_length_byaccuracy = [std_accuracy_0_step_length, std_accuracy_1_step_length, std_accuracy_2_step_length];
    
    figure;
    b = bar(average_step_length_byaccuracy, 'FaceColor', 'flat');
    hold on;
    
    for k = 1:length(average_speeds)
        b.CData(k, :) = colors(k, :);
    end
    
    errorbar(1:3, average_step_length_byaccuracy, std_step_length_byaccuracy, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    
    set(gca, 'xticklabel', {'None', 'Medium', 'High'});
    xlabel('Accuracy Prompt');
    ylabel('Average Step Length (cm)');
    title('Average Step Length by Accuracy Prompt - ' + bmh_text)
    ax = gca;               % Get current axis
    ax.FontSize = 14;       % Set axis number font size (e.g., 12);
    set(gca, 'FontSize', 22);
    hold off;
    
    % Step width variability, which is the std of each single trial
    avg_accuracy_0_variability = mean(accuracy_0_variability)/10; %From mm to cm
    avg_accuracy_1_variability = mean(accuracy_1_variability)/10; %From mm to cm
    avg_accuracy_2_variability = mean(accuracy_2_variability)/10; %From mm to cm
    % std_accuracy_0_variability = std(accuracy_0_variability)/10; %From mm to cm
    % std_accuracy_1_variability = std(accuracy_1_variability)/10; %From mm to cm
    % std_accuracy_2_variability = std(accuracy_2_variability)/10; %From mm to cm
    std_accuracy_0_variability = std(mean(reshape(accuracy_0_variability, 9, []), 2), 'omitnan')/10; %do first the mean of the stds each 9 positions, which correspond to the same condition for all participants
    std_accuracy_1_variability = std(mean(reshape(accuracy_1_variability, 9, []), 2), 'omitnan')/10;
    std_accuracy_2_variability = std(mean(reshape(accuracy_2_variability, 9, []), 2), 'omitnan')/10;
    
    average_step_length_variability_byaccuracy = [avg_accuracy_0_variability, avg_accuracy_1_variability, avg_accuracy_2_variability];
    std_step_length_variability_byaccuracy = [std_accuracy_0_variability, std_accuracy_1_variability, std_accuracy_2_variability];
    
    figure;
    b = bar(average_step_length_variability_byaccuracy, 'FaceColor', 'flat');
    hold on;
    
    for k = 1:length(average_speeds)
        b.CData(k, :) = colors(k, :);
    end
    
    errorbar(1:3, average_step_length_variability_byaccuracy, std_step_length_variability_byaccuracy, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    
    set(gca, 'xticklabel', {'None', 'Medium', 'High'});
    xlabel('Accuracy Prompt');
    ylabel('Average Step Length Variability (cm)');
    title('Average Step Length Variability by Accuracy Prompt - ' + bmh_text);
    ax = gca;               % Get current axis
    ax.FontSize = 14;       % Set axis number font size (e.g., 12)
    set(gca, 'FontSize', 22);
    hold off;

    %% Errors s3a2b3 across different balance conditions ---------------------- Lines plot (Plot 10) --------------------------------
    if shaded_plot == 0
        figure;
        hold on;
    
        balance_levels = {'None', 'Medium', 'High'}; 
    
        colors_2 = lines(length(bmh_list));  
        %softer_colors = colors * 0.6 + 0.4; 
        softer_colors = colors_participants * 0.3 + 0.7; % Softer colors 
        all_participants_errors = []; 
    
        for i = 1:length(bmh_list) % For each bmh
            bmh = bmh_list(i);
            s3a2b0 = [];
            s3a2b1 = [];
            s3a2b2 = [];
    
            for j = 3:29
                if double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 && double(imported_data.(bmh).prompt.balance(j,2)) == 0
                    s3a2b0 = [s3a2b0, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10]; %From mm to cm
                elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 && double(imported_data.(bmh).prompt.balance(j,2)) == 1
                    s3a2b1 = [s3a2b1, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10]; %From mm to cm
                elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 && double(imported_data.(bmh).prompt.balance(j,2)) == 2
                    s3a2b2 = [s3a2b2, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10]; %From mm to cm
                end
            end
    
            errors_std = [std(s3a2b0), std(s3a2b1), std(s3a2b2)];
            errors_participant = [mean(s3a2b0), mean(s3a2b1), mean(s3a2b2)];
            errors_participant = reshape(errors_participant, 1, []); 
            all_participants_errors = [all_participants_errors; errors_participant];
            errorbar(1:3, errors_participant, errors_std, '--o', 'LineWidth', 2.5, 'MarkerSize', 5, 'Color', softer_colors(i, :), 'MarkerFaceColor', softer_colors(i, :)); 
        end
    
        % Average calculation
        average_errors = mean(all_participants_errors, 1);
        average_std = std(all_participants_errors, 0, 1);
        red_line = errorbar(1:3, average_errors, average_std, '--o', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0]); % Red average line
        bmh_list = [bmh_list; "Average"]; % Add average for the legend
    
        %legend(red_line, 'Average', 'FontSize', 14); % Adjust the font size as needed
        legend(red_line, 'Average', 'FontSize', 14); % Only show the red average line
    
        % No baselines because here we are sorting by balance, not speed
        title("Step Error for all Participants for High Accuracy Condition")
    
        set(gca, 'FontSize', 22);
        xlabel('Visual Perturbation');
        ylabel('Step Error (cm)');
        xticks(1:3);
        xticklabels(balance_levels);
        bmh_list(strcmp(bmh_list, "Average")) = []; % Remove average to bmh_list
        grid on;
        hold off;
    else

        figure;
        hold on;
        balance_levels = {'None', 'Medium', 'High'}; 
        
        colors_2 = lines(length(bmh_list));  
        softer_colors = colors_participants * 0.3 + 0.7;
        all_participants_errors = []; 
        
        legend_handles = []; 
        
        for i = 1:length(bmh_list) 
            bmh = bmh_list(i);
            s3a2b0 = [];
            s3a2b1 = [];
            s3a2b2 = [];
        
            for j = 3:29
                if double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 && double(imported_data.(bmh).prompt.balance(j,2)) == 0
                    s3a2b0 = [s3a2b0, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10]; % De mm a cm
                elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 && double(imported_data.(bmh).prompt.balance(j,2)) == 1
                    s3a2b1 = [s3a2b1, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10]; % De mm a cm
                elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 && double(imported_data.(bmh).prompt.balance(j,2)) == 2
                    s3a2b2 = [s3a2b2, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10]; % De mm a cm
                end
            end
        
            errors_std = [std(s3a2b0), std(s3a2b1), std(s3a2b2)];
            errors_participant = [mean(s3a2b0), mean(s3a2b1), mean(s3a2b2)];
            errors_participant = reshape(errors_participant, 1, []); 
            all_participants_errors = [all_participants_errors; errors_participant];
            
            x = 1:3;
            fill([x, fliplr(x)], [errors_participant + errors_std, fliplr(errors_participant - errors_std)], ...
                 softer_colors(i, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        
            h = plot(1:3, errors_participant, '--o', 'LineWidth', 2.5, 'MarkerSize', 5, ...
                 'Color', softer_colors(i, :), 'MarkerFaceColor', softer_colors(i, :), 'DisplayName', string(bmh)); 
            legend_handles = [legend_handles, h]; 
        end
        
        average_errors = mean(all_participants_errors, 1);
        average_std = std(all_participants_errors, 0, 1);
        
        % Shaded line for the average
        %fill([x, fliplr(x)], [average_errors + average_std, fliplr(average_errors - average_std)], [1, 0, 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        
        % Mean
        h_avg = plot(1:3, average_errors, '--o', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0], 'DisplayName', 'Average');
        legend_handles = [legend_handles, h_avg]; 
        
        legend(legend_handles, 'FontSize', 14);
        title("Step Error for all Participants for High Accuracy Condition");
        set(gca, 'FontSize', 22);
        xlabel('Visual Perturbation');
        ylabel('Step Error (cm)');
        xticks(1:3);
        xticklabels(balance_levels);
        grid on;
        hold off;
    end

    %% Errors s3a1b3 across different balance conditions - PLOT 11

    figure
    hold on
    
    balance_levels = {'None', 'Medium', 'High'};  
    colors_2 = lines(length(bmh_list));  
    %softer_colors = colors * 0.6 + 0.4; 
    softer_colors = colors_participants * 0.3 + 0.7; % Softer colors
    all_participants_errors = []; 

    legend_handles = []; 
    
    for i = 1:length(bmh_list) % One iteration for each BMH
        bmh = bmh_list(i);
        s3a1b0 = [];
        s3a1b1 = [];
        s3a1b2 = [];

        for j = 3:29
            if double(imported_data.(bmh).prompt.accuracy(j,2)) == 1 && double(imported_data.(bmh).prompt.balance(j,2)) == 0
                s3a1b0 = [s3a1b0, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10]; % From mm to cm
            elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 1 && double(imported_data.(bmh).prompt.balance(j,2)) == 1
                s3a1b1 = [s3a1b1, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10]; % From mm to cm
            elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 1 && double(imported_data.(bmh).prompt.balance(j,2)) == 2
                s3a1b2 = [s3a1b2, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10]; % From mm to cm
            end
        end
        
        errors_std = [std(s3a1b0), std(s3a1b1), std(s3a1b2)]; % std calculation before subtracting the baseline(s)
        errors_participant = [mean(s3a1b0), mean(s3a1b1), mean(s3a1b2)];
        all_participants_errors = [all_participants_errors; errors_participant];        
        if shaded_plot == 0
            errorbar(1:3, errors_participant, errors_std, '--o', 'LineWidth', 2.5, 'MarkerSize', 5, 'Color', softer_colors(i, :), 'MarkerFaceColor', softer_colors(i, :)); 
        else
            x = 1:3;
            fill([x, fliplr(x)], [errors_participant + errors_std, fliplr(errors_participant - errors_std)], ...
                 softer_colors(i, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        
            h = plot(1:3, errors_participant, '--o', 'LineWidth', 2.5, 'MarkerSize', 5, ...
                 'Color', softer_colors(i, :), 'MarkerFaceColor', softer_colors(i, :), 'DisplayName', string(bmh)); 
            legend_handles = [legend_handles, h]; 
        end
    end
    
    % Average calculation
    average_errors = mean(all_participants_errors, 1);
    average_std = std(all_participants_errors, 0, 1);
    if shaded_plot == 0
        red_line = errorbar(1:3, average_errors, average_std, '--o', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0]); % Red average line
        bmh_list = [bmh_list; "Average"]; % Add average for the legend
        legend(red_line, 'Average', 'FontSize', 14); % Adjust the font size as needed
        bmh_list(strcmp(bmh_list, "Average")) = []; % Remove average to bmh_list
    else
        h_avg = plot(1:3, average_errors, '--o', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0], 'DisplayName', 'Average');
        legend_handles = [legend_handles, h_avg];
        legend(legend_handles, 'FontSize', 14);
    end

    % No baselines because here we are sorting by balance, not speed
    title("Step Error for all Participants for Medium Accuracy Condition")
    
    set(gca, 'FontSize', 22);
    xlabel('Visual Perturbation');
    ylabel('Step Error (cm)');
    xticks(1:3);
    ylim([0 32])
    xticklabels(balance_levels);
    grid on;
    hold off

    %% Errors s3a2b0 across different speed condition - PLOT 12

    figure
    hold on
    
    speed_levels = {'Slow', 'Medium', 'Fast'}; 
    colors_2 = lines(length(bmh_list));  
    %softer_colors = colors * 0.6 + 0.4; 
    softer_colors = colors_participants * 0.3 + 0.7; % Softer colors 
    all_participants_errors = []; 

    legend_handles = []; 
    
    for i = 1:length(bmh_list)
        bmh = bmh_list(i);
        s0a2b0_ = [];
        s1a2b0_ = [];
        s2a2b0_ = [];

        for j = 3:29
            if double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 && double(imported_data.(bmh).prompt.speed(j,2)) == 0 && double(imported_data.(bmh).prompt.balance(j,2)) == 0
                s0a2b0_ = [s0a2b0_, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10];
            elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 && double(imported_data.(bmh).prompt.speed(j,2)) == 1 && double(imported_data.(bmh).prompt.balance(j,2)) == 0
                s1a2b0_ = [s1a2b0_, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10];
            elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 && double(imported_data.(bmh).prompt.speed(j,2)) == 2 && double(imported_data.(bmh).prompt.balance(j,2)) == 0
                s2a2b0_ = [s2a2b0_, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10];
            end
        end
        
        errors_std = [std(s0a2b0_), std(s1a2b0_), std(s2a2b0_)]; % std calculation before subtracting the baseline(s)

        if subtract_all == 1
            val = (s0a2b0 + s1a2b0 + s2a2b0)/3;
            s0a2b0_ = abs(s0a2b0_ - val);
            s1a2b0_ = abs(s1a2b0_ - val);
            s2a2b0_ = abs(s2a2b0_ - val);
        elseif subtract_speeds == 1
            % s0a2b0_ = abs(s0a2b0_ - (baseline_errors.(bmh).speed_s0a2b0)/10);
            % s1a2b0_ = abs(s1a2b0_ - (baseline_errors.(bmh).speed_s1a2b0)/10);
            % s2a2b0_ = abs(s2a2b0_ - (baseline_errors.(bmh).speed_s2a2b0)/10);
            s0a2b0_ = (s0a2b0_ - (baseline_errors.(bmh).speed_s0a2b0)/10);
            s1a2b0_ = (s1a2b0_ - (baseline_errors.(bmh).speed_s1a2b0)/10);
            s2a2b0_ = (s2a2b0_ - (baseline_errors.(bmh).speed_s2a2b0)/10);
        end

        errors_participant = [mean(s0a2b0_), mean(s1a2b0_), mean(s2a2b0_)];
        %plot(1:3, errors_participant, '--o', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', colors(i, :), 'MarkerFaceColor', colors(i, :));
        all_participants_errors = [all_participants_errors; errors_participant];        
        if shaded_plot == 0 
            errorbar(1:3, errors_participant, errors_std, '--o', 'LineWidth', 2.5, 'MarkerSize', 5, 'Color', softer_colors(i, :), 'MarkerFaceColor', softer_colors(i, :)); 
        else
            x = 1:3;
            fill([x, fliplr(x)], [errors_participant + errors_std, fliplr(errors_participant - errors_std)], ...
                 softer_colors(i, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        
            h = plot(1:3, errors_participant, '--o', 'LineWidth', 2.5, 'MarkerSize', 5, ...
                 'Color', softer_colors(i, :), 'MarkerFaceColor', softer_colors(i, :), 'DisplayName', string(bmh)); 
            legend_handles = [legend_handles, h]; 
        end
    end
    
    if subtract_all == 1
        title("Step Error for all Participants for no Balance Disturbace and High Accuracy - Average Error Subtracted")
    elseif subtract_speeds == 1
        title("Step Error for all Participants for no Balance Disturbace and High Accuracy - Every Speed Error Subtracted")
    else
        title("Step Error for all Participants for no Balance Disturbace and High Accuracy")
    end
    
    % Average calculation
    average_errors = mean(all_participants_errors, 1); 
    average_std = std(all_participants_errors, 0, 1); 
    if shaded_plot == 0
        red_line = errorbar(1:3, average_errors, average_std, '--o', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0]); % Red average line
        bmh_list = [bmh_list; "Average"]; % Add average for the legend
        legend(red_line, 'Average', 'FontSize', 14); % Adjust the font size as needed
        bmh_list(strcmp(bmh_list, "Average")) = []; % Remove average to bmh_list
    else
        h_avg = plot(1:3, average_errors, '--o', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0], 'DisplayName', 'Average');
        legend_handles = [legend_handles, h_avg];
        legend(legend_handles, 'FontSize', 14);
    end

    set(gca, 'FontSize', 22);
    xlabel('Speed Prompt');
    ylabel('Step Error (cm)');
    xticks(1:3);
    xticklabels(speed_levels);
    grid on;
    hold off

    %% Errors s3a2b1 across different speed condition - PLOT 13

    figure
    hold on
    
    speed_levels = {'Slow', 'Medium', 'Fast'}; 
    colors_2 = lines(length(bmh_list));  
    %softer_colors = colors * 0.6 + 0.4; 
    softer_colors = colors_participants * 0.3 + 0.7; % Softer colors 
    all_participants_errors = [];   
    legend_handles = []; 
    
    for i = 1:length(bmh_list)
        bmh = bmh_list(i);
        s0a2b1 = [];
        s1a2b1 = [];
        s2a2b1 = [];

        for j = 3:29
            if double(imported_data.(bmh).prompt.balance(j,2)) == 1
                if double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 && double(imported_data.(bmh).prompt.speed(j,2)) == 0
                    s0a2b1 = [s0a2b1, imported_data.(bmh).overall_errors.straight_mean_absolute(j)];
                elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 && double(imported_data.(bmh).prompt.speed(j,2)) == 1
                    s1a2b1 = [s1a2b1, imported_data.(bmh).overall_errors.straight_mean_absolute(j)];
                elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 && double(imported_data.(bmh).prompt.speed(j,2)) == 2
                    s2a2b1 = [s2a2b1, imported_data.(bmh).overall_errors.straight_mean_absolute(j)];
                end
            end
        end
        
        errors_std = [std(s0a2b1), std(s1a2b1), std(s2a2b1)]; % std calculation before subtracting the baseline(s)

        if subtract_all == 1
            val = (s0a2b0 + s1a2b0 + s2a2b0)/3;
            s0a2b1 = abs(s0a2b1 - val);
            s1a2b1 = abs(s1a2b1 - val);
            s2a2b1 = abs(s2a2b1 - val);
        elseif subtract_speeds == 1
            % s0a2b1 = abs(s0a2b1/10 - (baseline_errors.(bmh).speed_s0a2b0)/10);
            % s1a2b1 = abs(s1a2b1/10 - (baseline_errors.(bmh).speed_s1a2b0)/10);
            % s2a2b1 = abs(s2a2b1/10 - (baseline_errors.(bmh).speed_s2a2b0)/10);
            s0a2b1 = (s0a2b1/10 - (baseline_errors.(bmh).speed_s0a2b0)/10);
            s1a2b1 = (s1a2b1/10 - (baseline_errors.(bmh).speed_s1a2b0)/10);
            s2a2b1 = (s2a2b1/10 - (baseline_errors.(bmh).speed_s2a2b0)/10);
        end

        errors_participant = [mean(s0a2b1), mean(s1a2b1), mean(s2a2b1)]/10; % From mm to cm
        
        %plot(1:3, errors_participant, '--o', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', colors(i, :), 'MarkerFaceColor', colors(i, :));
        all_participants_errors = [all_participants_errors; errors_participant];        
        if shaded_plot == 0 
            errorbar(1:3, errors_participant, errors_std, '--o', 'LineWidth', 2.5, 'MarkerSize', 5, 'Color', softer_colors(i, :), 'MarkerFaceColor', softer_colors(i, :)); 
        else
            x = 1:3;
            fill([x, fliplr(x)], [errors_participant + errors_std, fliplr(errors_participant - errors_std)], ...
                 softer_colors(i, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        
            h = plot(1:3, errors_participant, '--o', 'LineWidth', 2.5, 'MarkerSize', 5, ...
                 'Color', softer_colors(i, :), 'MarkerFaceColor', softer_colors(i, :), 'DisplayName', string(bmh)); 
            legend_handles = [legend_handles, h]; 
        end
    end
    
    if subtract_all == 1
        title("Step Error for all Participants for Medium Balance Disturbace and High Accuracy - Average Error Subtracted")
    elseif subtract_speeds == 1
        title("Step Error for all Participants for Medium Balance Disturbace and High Accuracy - Every Speed Error Subtracted")
    else
        title("Step Error for all Participants for Medium Balance Disturbace and High Accuracy")
    end
    
    % Average calculation
    average_errors = mean(all_participants_errors, 1);
    average_std = std(all_participants_errors, 0, 1);
    if shaded_plot == 0
        red_line = errorbar(1:3, average_errors, average_std, '--o', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0]); % Red average line
        bmh_list = [bmh_list; "Average"]; % Add average for the legend
        legend(red_line, 'Average', 'FontSize', 14); % Adjust the font size as needed
        bmh_list(strcmp(bmh_list, "Average")) = []; % Remove average to bmh_list
    else
        h_avg = plot(1:3, average_errors, '--o', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0], 'DisplayName', 'Average');
        legend_handles = [legend_handles, h_avg];
        legend(legend_handles, 'FontSize', 14);
    end

    set(gca, 'FontSize', 22);
    xlabel('Speed Prompt');
    ylabel('Step Error (cm)');
    xticks(1:3);
    xticklabels(speed_levels);
    grid on;
    hold off

    %% Errors s3a2b2 across different speed condition - PLOT 14

    figure
    hold on
    
    speed_levels = {'Slow', 'Medium', 'Fast'}; 
    colors_2 = lines(length(bmh_list));  
    %softer_colors = colors * 0.6 + 0.4; 
    softer_colors = colors_participants * 0.3 + 0.7; % Softer colors
    all_participants_errors = []; 
    
    for i = 1:length(bmh_list)
        bmh = bmh_list(i);
        s0a2b2 = [];
        s1a2b2 = [];
        s2a2b2 = [];

        for j = 3:29
            if double(imported_data.(bmh).prompt.balance(j,2)) == 2
                if double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 && double(imported_data.(bmh).prompt.speed(j,2)) == 0
                    s0a2b2 = [s0a2b2, imported_data.(bmh).overall_errors.straight_mean_absolute(j)];
                elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 && double(imported_data.(bmh).prompt.speed(j,2)) == 1
                    s1a2b2 = [s1a2b2, imported_data.(bmh).overall_errors.straight_mean_absolute(j)];
                elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 && double(imported_data.(bmh).prompt.speed(j,2)) == 2
                    s2a2b2 = [s2a2b2, imported_data.(bmh).overall_errors.straight_mean_absolute(j)];
                end
            end
        end
        
        errors_std = [std(s0a2b2), std(s1a2b2), std(s2a2b2)]/10; % From mm to cm; std calculation before subtracting the baseline(s)

        if subtract_all == 1
            val = (s0a2b0 + s1a2b0 + s2a2b0)/3;
            s0a2b2 = abs(s0a2b2 - val);
            s1a2b2 = abs(s1a2b2 - val);
            s2a2b2 = abs(s2a2b2 - val);
        elseif subtract_speeds == 1
            % s0a2b2 = abs(s0a2b2/10 - (baseline_errors.(bmh).speed_s0a2b0)/10);
            % s1a2b2 = abs(s1a2b2/10 - (baseline_errors.(bmh).speed_s1a2b0)/10);
            % s2a2b2 = abs(s2a2b2/10 - (baseline_errors.(bmh).speed_s2a2b0)/10);
            s0a2b2 = (s0a2b2/10 - (baseline_errors.(bmh).speed_s0a2b0)/10);
            s1a2b2 = (s1a2b2/10 - (baseline_errors.(bmh).speed_s1a2b0)/10);
            s2a2b2 = (s2a2b2/10 - (baseline_errors.(bmh).speed_s2a2b0)/10);
        end

        errors_participant = [mean(s0a2b2), mean(s1a2b2), mean(s2a2b2)]/10; % From mm to cm
        
        %plot(1:3, errors_participant, '--o', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', colors(i, :), 'MarkerFaceColor', colors(i, :));
        all_participants_errors = [all_participants_errors; errors_participant];        
        errorbar(1:3, errors_participant, errors_std, '--o', 'LineWidth', 2.5, 'MarkerSize', 5, 'Color', softer_colors(i, :), 'MarkerFaceColor', softer_colors(i, :)); 
    end
    
    if subtract_all == 1
        title("Step Error for all Participants for High Balance Disturbace and High Accuracy - Average Error Subtracted")
    elseif subtract_speeds == 1
        title("Step Error for all Participants for High Balance Disturbace and High Accuracy - Every Speed Error Subtracted")
    else
        title("Step Error for all Participants for High Balance Disturbace and High Accuracy")
    end
    
    % Average calculation
    average_errors = mean(all_participants_errors, 1);
    average_std = std(all_participants_errors, 0, 1);
    red_line = errorbar(1:3, average_errors, average_std, '--o', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0]); % Red average line

    bmh_list = [bmh_list; "Average"]; % Add average for the legend
    legend(red_line, 'Average', 'FontSize', 14); % Adjust the font size as needed
    bmh_list(strcmp(bmh_list, "Average")) = []; % Remove average to bmh_list

    set(gca, 'FontSize', 22);
    xlabel('Speed Prompt');
    ylabel('Step Error (cm)');
    xticks(1:3);
    ylim([0 28.5])
    xticklabels(speed_levels);
    grid on;
    hold off

    %% Errors s3a1b0 across different speed condition - PLOT 15

    figure
    hold on
    
    speed_levels = {'Slow', 'Medium', 'Fast'}; 
    colors_2 = lines(length(bmh_list));  
    %softer_colors = colors * 0.6 + 0.4; 
    softer_colors = colors_participants * 0.3 + 0.7; % Softer colors
    all_participants_errors = [];
    
    for i = 1:length(bmh_list)
        bmh = bmh_list(i);
        s0a1b0 = [];
        s1a1b0 = [];
        s2a1b0 = [];

        for j = 3:29
            if double(imported_data.(bmh).prompt.balance(j,2)) == 0
                if double(imported_data.(bmh).prompt.accuracy(j,2)) == 1 && double(imported_data.(bmh).prompt.speed(j,2)) == 0
                    s0a1b0 = [s0a1b0, imported_data.(bmh).overall_errors.straight_mean_absolute(j)];
                elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 1 && double(imported_data.(bmh).prompt.speed(j,2)) == 1
                    s1a1b0 = [s1a1b0, imported_data.(bmh).overall_errors.straight_mean_absolute(j)];
                elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 1 && double(imported_data.(bmh).prompt.speed(j,2)) == 2
                    s2a1b0 = [s2a1b0, imported_data.(bmh).overall_errors.straight_mean_absolute(j)];
                end
            end
        end
        
        errors_participant = [mean(s0a1b0), mean(s1a1b0), mean(s2a1b0)]/10; % From mm to cm
        errors_std = [std(s0a1b0), std(s1a1b0), std(s2a1b0)]/10; % From mm to cm
        
        %plot(1:3, errors_participant, '--o', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', colors(i, :), 'MarkerFaceColor', colors(i, :));
        all_participants_errors = [all_participants_errors; errors_participant];        
        errorbar(1:3, errors_participant, errors_std, '--o', 'LineWidth', 2.5, 'MarkerSize', 5, 'Color', softer_colors(i, :), 'MarkerFaceColor', softer_colors(i, :)); 
    end
    
    % Average calculation
    average_errors = mean(all_participants_errors, 1);
    average_std = std(all_participants_errors, 0, 1);
    red_line = errorbar(1:3, average_errors, average_std, '--o', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0]); % Red average line

    bmh_list = [bmh_list; "Average"]; % Add average for the legend
    legend(red_line, 'Average', 'FontSize', 14); % Adjust the font size as needed
    bmh_list(strcmp(bmh_list, "Average")) = []; % Remove average to bmh_list

    % No baselines because this is medium accuracy, and baselines are for high accuracy
    title("Step Error for all Participants for no Balance Disturbace and Medium Accuracy")

    set(gca, 'FontSize', 22);
    %legend(bmh_list);
    ylim([0 31.5])
    xlabel('Speed Prompt');
    ylabel('Step Error (cm)');
    xticks(1:3);
    xticklabels(speed_levels);
    grid on;
    hold off

    %% Ignore lines when a=0? - PLOT 16
    
    fields = fieldnames(data_labels); % Get all the labels
    relevant_labels = fields(contains(fields, 'a0')); % Filter labels where accuracy 'a' is 0
    
    deviations = []; % To store deviations from the closest line
    
    for i = 1:length(relevant_labels)
        label = relevant_labels{i};
        stepping_data = data_labels.(label).step_error; % Extract stepping accuracies
        deviation = abs(stepping_data); % Absolute deviation (already relative to nearest line)
        
        % Store deviations
        deviations = [deviations; deviation];
    end
    
    % Perform statistical test (e.g., one-sample t-test against zero deviation)
    [h, p, ci, stats] = ttest(deviations, 0); % Null hypothesis: mean deviation = 0
    
    % Display test results
    disp('One-sample t-test results:');
    disp(['p-value: ', num2str(p)]);
    disp(['Mean deviation: ', num2str(mean(deviations))]);
    disp(['Confidence interval: [', num2str(ci(1)), ', ', num2str(ci(2)), ']']);
    disp('');

    % Small summary for interpretation:

        %p-value: if it is small, it indicates that the null hypothesis (mean deviation = 0) is rejected with very high confidence. This means that 
        %there is a significant deviation from 0 cm, i.e., participants are not consistently stepping on the lines and are ignoring them (on average).

        %Mean Deviation: This is the average deviation of stepping accuracy (mean error) regarding the  closest line (60 cm apart).

        %Confidence Interval: This is the range within which the true mean deviation is likely to fall.
    
    figure;
    set(gca, 'FontSize', 22);
    boxplot(deviations, 'Positions', 0);
    set(gca, 'XTick', 0, 'XTickLabel', {'a = 0'});
    title('Deviation from Nearest Line (Ignoring Lines)');
    xlabel('Accuracy Condition');
    ylabel('Deviation (cm)');

    %% Walking parameters headset vs no headset s1a0b0  - PLOT 17

    s1a0b0_headset = struct ();
    s1a0b0_headset.meanwidthstraights = [];
    s1a0b0_headset.meanwidthcurves = [];
    s1a0b0_headset.straightsvariability = [];
    s1a0b0_headset.curvesvariability = [];
    s1a0b0_headset.meanlength_straights = [];
    s1a0b0_headset.meanlength_curves = [];
    s1a0b0_headset.meanlength_straights_variability = [];
    s1a0b0_headset.meanlength_curves_variability = [];

    s1a0b0_no_headset = struct ();
    s1a0b0_no_headset.meanwidthstraights = [];
    s1a0b0_no_headset.meanwidthcurves = [];
    s1a0b0_no_headset.straightsvariability = [];
    s1a0b0_no_headset.curvesvariability = [];
    s1a0b0_no_headset.meanlength_straights = [];
    s1a0b0_no_headset.meanlength_curves = [];
    s1a0b0_no_headset.meanlength_straights_variability = [];
    s1a0b0_no_headset.meanlength_curves_variability = [];

    %Trial 1 and 35: Walk at your typical speed and ignore the red lines on the ground (no headset)
    %Trial 2 and 34: Walk at your typical speed and ignore the red lines on the ground (headset)
    for i=1:size(bmh_list,1)
        bmh = bmh_list(i);

        s1a0b0_headset.meanwidthstraights = [s1a0b0_headset.meanwidthstraights, imported_data.(bmh).meanwidthstraights([headset_beginning(i), headset_end(i)])]; % Take only trials 1 and 35
        s1a0b0_headset.meanwidthcurves = [s1a0b0_headset.meanwidthcurves, imported_data.(bmh).meanwidthcurves([headset_beginning(i), headset_end(i)])];
        s1a0b0_headset.straightsvariability = [s1a0b0_headset.straightsvariability, imported_data.(bmh).straightsvariability([headset_beginning(i), headset_end(i)])];
        s1a0b0_headset.curvesvariability = [s1a0b0_headset.curvesvariability, imported_data.(bmh).curvesvariability([headset_beginning(i), headset_end(i)])];
        s1a0b0_headset.meanlength_straights = [s1a0b0_headset.meanlength_straights, imported_data.(bmh).meanlength_straights([headset_beginning(i), headset_end(i)])];
        s1a0b0_headset.meanlength_curves = [s1a0b0_headset.meanlength_curves, imported_data.(bmh).meanlength_curves([headset_beginning(i), headset_end(i)])];
        s1a0b0_headset.meanlength_straights_variability = [s1a0b0_headset.meanlength_straights_variability, imported_data.(bmh).meanlength_straights_variability([headset_beginning(i), headset_end(i)])];
        s1a0b0_headset.meanlength_curves_variability = [s1a0b0_headset.meanlength_curves_variability, imported_data.(bmh).meanlength_curves_variability([headset_beginning(i), headset_end(i)])];
    
        s1a0b0_no_headset.meanwidthstraights = [s1a0b0_no_headset.meanwidthstraights, imported_data.(bmh).meanwidthstraights([no_headset_beginning(i), no_headset_end(i)])]; % Take only trials 2 and 34
        s1a0b0_no_headset.meanwidthcurves = [s1a0b0_no_headset.meanwidthcurves, imported_data.(bmh).meanwidthcurves([no_headset_beginning(i), no_headset_end(i)])];
        s1a0b0_no_headset.straightsvariability = [s1a0b0_no_headset.straightsvariability, imported_data.(bmh).straightsvariability([no_headset_beginning(i), no_headset_end(i)])];
        s1a0b0_no_headset.curvesvariability = [s1a0b0_no_headset.curvesvariability, imported_data.(bmh).curvesvariability([no_headset_beginning(i), no_headset_end(i)])];
        s1a0b0_no_headset.meanlength_straights = [s1a0b0_no_headset.meanlength_straights, imported_data.(bmh).meanlength_straights([no_headset_beginning(i), no_headset_end(i)])];
        s1a0b0_no_headset.meanlength_curves = [s1a0b0_no_headset.meanlength_curves, imported_data.(bmh).meanlength_curves([no_headset_beginning(i), no_headset_end(i)])];
        s1a0b0_no_headset.meanlength_straights_variability = [s1a0b0_no_headset.meanlength_straights_variability, imported_data.(bmh).meanlength_straights_variability([no_headset_beginning(i), no_headset_end(i)])];
        s1a0b0_no_headset.meanlength_curves_variability = [s1a0b0_no_headset.meanlength_curves_variability, imported_data.(bmh).meanlength_curves_variability([no_headset_beginning(i), no_headset_end(i)])];

    end

    s1a0b0_headset.meanwidthstraights_avg = mean(s1a0b0_headset.meanwidthstraights)/10; % From mm to cm;
    s1a0b0_headset.meanwidthcurves_avg = mean(s1a0b0_headset.meanwidthcurves)/10; % From mm to cm
    s1a0b0_headset.straightsvariability_avg = mean(s1a0b0_headset.straightsvariability)/10; % From mm to cm
    s1a0b0_headset.curvesvariability_avg = mean(s1a0b0_headset.curvesvariability)/10; % From mm to cm
    s1a0b0_headset.meanlength_straights_avg = mean(s1a0b0_headset.meanlength_straights)/10; % From mm to cm
    s1a0b0_headset.meanlength_curves_avg = mean(s1a0b0_headset.meanlength_curves)/10; % From mm to cm
    s1a0b0_headset.meanlength_straights_variability_avg = mean(s1a0b0_headset.meanlength_straights_variability)/10; % From mm to cm
    s1a0b0_headset.meanlength_curves_variability_avg = mean(s1a0b0_headset.meanlength_curves_variability)/10; % From mm to cm
    
    s1a0b0_no_headset.meanwidthstraights_avg = mean(s1a0b0_no_headset.meanwidthstraights)/10; % From mm to cm
    s1a0b0_no_headset.meanwidthcurves_avg = mean(s1a0b0_no_headset.meanwidthcurves)/10; % From mm to cm
    s1a0b0_no_headset.straightsvariability_avg = mean(s1a0b0_no_headset.straightsvariability)/10; % From mm to cm
    s1a0b0_no_headset.curvesvariability_avg = mean(s1a0b0_no_headset.curvesvariability)/10; % From mm to cm
    s1a0b0_no_headset.meanlength_straights_avg = mean(s1a0b0_no_headset.meanlength_straights)/10; % From mm to cm
    s1a0b0_no_headset.meanlength_curves_avg = mean(s1a0b0_no_headset.meanlength_curves)/10; % From mm to cm
    s1a0b0_no_headset.meanlength_straights_variability_avg = mean(s1a0b0_no_headset.meanlength_straights_variability)/10; % From mm to cm
    s1a0b0_no_headset.meanlength_curves_variability_avg = mean(s1a0b0_no_headset.meanlength_curves_variability)/10; % From mm to cm

    s1a0b0_headset.meanwidthstraights_std = std(s1a0b0_headset.meanwidthstraights)/10; % From mm to cm
    s1a0b0_headset.meanwidthcurves_avg_std = std(s1a0b0_headset.meanwidthcurves)/10; % From mm to cm
    s1a0b0_headset.straightsvariability_avg_std = std(s1a0b0_headset.straightsvariability)/10; % From mm to cm
    s1a0b0_headset.curvesvariability_avg_std = std(s1a0b0_headset.curvesvariability)/10; % From mm to cm
    s1a0b0_headset.meanlength_straights_avg_std = std(s1a0b0_headset.meanlength_straights)/10; % From mm to cm
    s1a0b0_headset.meanlength_curves_avg_std = std(s1a0b0_headset.meanlength_curves)/10; % From mm to cm
    s1a0b0_headset.meanlength_straights_variability_avg_std = std(s1a0b0_headset.meanlength_straights_variability)/10; % From mm to cm
    s1a0b0_headset.meanlength_curves_variability_avg_std = std(s1a0b0_headset.meanlength_curves_variability)/10; % From mm to cm
    
    s1a0b0_no_headset.meanwidthstraights_std = std(s1a0b0_no_headset.meanwidthstraights)/10; % From mm to cm
    s1a0b0_no_headset.meanwidthcurves_avg_std = std(s1a0b0_no_headset.meanwidthcurves)/10; % From mm to cm
    s1a0b0_no_headset.straightsvariability_avg_std = std(s1a0b0_no_headset.straightsvariability)/10; % From mm to cm
    s1a0b0_no_headset.curvesvariability_avg_std = std(s1a0b0_no_headset.curvesvariability)/10; % From mm to cm
    s1a0b0_no_headset.meanlength_straights_avg_std = std(s1a0b0_no_headset.meanlength_straights)/10; % From mm to cm
    s1a0b0_no_headset.meanlength_curves_avg_std = std(s1a0b0_no_headset.meanlength_curves)/10; % From mm to cm
    s1a0b0_no_headset.meanlength_straights_variability_avg_std = std(s1a0b0_no_headset.meanlength_straights_variability)/10; % From mm to cm
    s1a0b0_no_headset.meanlength_curves_variability_avg_std = std(s1a0b0_no_headset.meanlength_curves_variability)/10; % From mm to cm

    headset_avgs = [
        s1a0b0_headset.meanwidthstraights_avg, ...
        s1a0b0_headset.straightsvariability_avg, ...
        s1a0b0_headset.meanwidthcurves_avg, ...
        s1a0b0_headset.curvesvariability_avg, ...
        s1a0b0_headset.meanlength_straights_avg, ...
        s1a0b0_headset.meanlength_straights_variability_avg, ...
        s1a0b0_headset.meanlength_curves_avg, ...
        s1a0b0_headset.meanlength_curves_variability_avg
    ];
    
    no_headset_avgs = [
        s1a0b0_no_headset.meanwidthstraights_avg, ...
        s1a0b0_no_headset.straightsvariability_avg, ...
        s1a0b0_no_headset.meanwidthcurves_avg, ...
        s1a0b0_no_headset.curvesvariability_avg, ...
        s1a0b0_no_headset.meanlength_straights_avg, ...
        s1a0b0_no_headset.meanlength_straights_variability_avg, ...
        s1a0b0_no_headset.meanlength_curves_avg, ...
        s1a0b0_no_headset.meanlength_curves_variability_avg
    ];
    
    headset_stds = [
        s1a0b0_headset.meanwidthstraights_std, ...
        s1a0b0_headset.straightsvariability_avg_std, ...
        s1a0b0_headset.meanwidthcurves_avg_std, ...
        s1a0b0_headset.curvesvariability_avg_std, ...
        s1a0b0_headset.meanlength_straights_avg_std, ...
        s1a0b0_headset.meanlength_straights_variability_avg_std, ...
        s1a0b0_headset.meanlength_curves_avg_std, ...
        s1a0b0_headset.meanlength_curves_variability_avg_std
    ];
    
    no_headset_stds = [
        s1a0b0_no_headset.meanwidthstraights_std, ...
        s1a0b0_no_headset.straightsvariability_avg_std, ...
        s1a0b0_no_headset.meanwidthcurves_avg_std, ...
        s1a0b0_no_headset.curvesvariability_avg_std, ...
        s1a0b0_no_headset.meanlength_straights_avg_std, ...
        s1a0b0_no_headset.meanlength_straights_variability_avg_std, ...
        s1a0b0_no_headset.meanlength_curves_avg_std, ...
        s1a0b0_no_headset.meanlength_curves_variability_avg_std
    ];
   
    variable_names = {
    'SW Straights', ...
    'SWV Straights', ...
    'SW Curves', ...
    'SWV Curves', ...
    'SL Straights', ...
    'SLV Straights', ...
    'SL Curves', ...
    'SLV Curves'
    };

    left_indices = [5, 7]; % Step Length Straights and Curves
    right_indices = [6, 8, 1, 2, 3, 4]; % Step Length Variability, Step Width, and Width Variability
    
    left_data = [headset_avgs(left_indices); no_headset_avgs(left_indices)]';
    left_stds = [headset_stds(left_indices); no_headset_stds(left_indices)]';
    right_data = [headset_avgs(right_indices); no_headset_avgs(right_indices)]';
    right_stds = [headset_stds(right_indices); no_headset_stds(right_indices)]';
    
    figure;

    % --- Left subplot: Step Length (Straights and Curves) ---
    subplot('Position', [0.05, 0.3, 0.25, 0.6]); % Ajustar posición del subplot izquierdo
    bar_handle_left = bar(left_data, 'grouped');
    hold on;
    
    % Error bars for left plot
    num_groups_left = size(left_data, 1);
    num_bars_left = size(left_data, 2);
    group_width_left = min(0.8, num_bars_left / (num_bars_left + 1.5));
    
    x_left = zeros(num_groups_left, num_bars_left);
    for i = 1:num_bars_left
        x_left(:, i) = (1:num_groups_left) - group_width_left / 2 + ...
                       (2 * i - 1) * group_width_left / (2 * num_bars_left);
    end
    
    % Add error bars
    errorbar(x_left(:, 1), left_data(:, 1), left_stds(:, 1), 'k.', 'LineWidth', 1.5); % Headset
    errorbar(x_left(:, 2), left_data(:, 2), left_stds(:, 2), 'k.', 'LineWidth', 1.5); % No Headset
    hold off;
    
    % Customize left plot
    set(gca, 'FontSize', 22);
    xticks(1:num_groups_left);
    xticklabels(variable_names(left_indices));
    xtickangle(45);
    ylabel('Mean Value (cm)');
    %title('Step Length Metrics');
    title('Step Length');
    grid on;
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.YTick = min(ylim):10:82; % Grid horizontal cada 0.2
    bar_handle_left(1).FaceColor = [0.2, 0.6, 0.8]; % Headset color
    bar_handle_left(2).FaceColor = [0.9, 0.4, 0.4]; % No Headset color
    
    % --- Right subplot: Step Length Variability, Step Width, and Width Variability ---
    subplot('Position', [0.4, 0.3, 0.55, 0.6]); % Ajustar posición del subplot derecho
    bar_handle_right = bar(right_data, 'grouped');
    hold on;
    
    % Error bars for right plot
    num_groups_right = size(right_data, 1);
    num_bars_right = size(right_data, 2);
    group_width_right = min(0.8, num_bars_right / (num_bars_right + 1.5));
    
    x_right = zeros(num_groups_right, num_bars_right);
    for i = 1:num_bars_right
        x_right(:, i) = (1:num_groups_right) - group_width_right / 2 + ...
                        (2 * i - 1) * group_width_right / (2 * num_bars_right);
    end
    
    % Add error bars
    errorbar(x_right(:, 1), right_data(:, 1), right_stds(:, 1), 'k.', 'LineWidth', 1.5); % Headset
    errorbar(x_right(:, 2), right_data(:, 2), right_stds(:, 2), 'k.', 'LineWidth', 1.5); % No Headset
    hold off;
    
    % Customize right plot
    set(gca, 'FontSize', 22);
    xticks(1:num_groups_right);
    xticklabels(variable_names(right_indices));
    xtickangle(45);
    ylabel('Mean Value (cm)');
    %title('Step Length Variability and Step Width Metrics');
    title('Various Parameters');
    grid on;
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.YTick = min(ylim):2:max(ylim); % Grid horizontal cada 0.2
    legend({'Headset', 'No Headset'}, 'Location', 'northwest');
    bar_handle_right(1).FaceColor = [0.2, 0.6, 0.8]; % Headset color
    bar_handle_right(2).FaceColor = [0.9, 0.4, 0.4]; % No Headset color
    
    % Resize figure for better visualization
    set(gcf, 'Position', [100, 100, 1400, 700]); % Ampliar el tamaño general de la figura
    set(gca, 'LooseInset', [0.1, 0.1, 0.1, 0.2]); % Ajustar márgenes para evitar texto cortado
    
    % Calculate and display percentage variation
    percentage_variation = ((no_headset_avgs - headset_avgs) ./ headset_avgs) * 100;
    disp('Percentage Variation Between Headset and No Headset for Each Variable:');
    for i = 1:length(variable_names)
        fprintf('%s: %.2f%%\n', variable_names{i}, percentage_variation(i));
    end

    %% Errors speed vs balance when a=2 - PLOT 18

    figure
    hold on
    
    balance_levels = {'None', 'Medium', 'High'}; 
    colors_2 = lines(length(bmh_list));  
    softer_colors = colors_participants * 0.3 + 0.7; % Softer colors 
    all_participants_errors = []; 
    legend_handles = []; 
    
    for i = 1:length(bmh_list) % One iteration for each BMH
        bmh = bmh_list(i);
        s3a2b0 = [];
        s3a2b1 = [];
        s3a2b2 = [];

        for j = 3:29
            if double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 && double(imported_data.(bmh).prompt.balance(j,2)) == 0
                s3a2b0 = [s3a2b0, imported_data.(bmh).walkingspeed(j)];
            elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 && double(imported_data.(bmh).prompt.balance(j,2)) == 1
                s3a2b1 = [s3a2b1, imported_data.(bmh).walkingspeed(j)];
            elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 && double(imported_data.(bmh).prompt.balance(j,2)) == 2
                s3a2b2 = [s3a2b2, imported_data.(bmh).walkingspeed(j)];
            end
        end
        
        speed_std = [std(s3a2b0), std(s3a2b1), std(s3a2b2)]; % std calculation before subtracting the baseline(s)
        speed_participant = [mean(s3a2b0), mean(s3a2b1), mean(s3a2b2)];
                
        %plot(1:3, errors_participant, '--o', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', colors(i, :), 'MarkerFaceColor', colors(i, :));
        all_participants_errors = [all_participants_errors; speed_participant];        
        if shaded_plot == 0 
            errorbar(1:3, speed_participant, speed_std, '--o', 'LineWidth', 2.5, 'MarkerSize', 5, 'Color', softer_colors(i, :), 'MarkerFaceColor', softer_colors(i, :)); 
        else
            x = 1:3;
            fill([x, fliplr(x)], [speed_participant + speed_std, fliplr(speed_participant - speed_std)], ...
                 softer_colors(i, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        
            h = plot(1:3, speed_participant, '--o', 'LineWidth', 2.5, 'MarkerSize', 5, ...
                 'Color', softer_colors(i, :), 'MarkerFaceColor', softer_colors(i, :), 'DisplayName', string(bmh)); 
            legend_handles = [legend_handles, h]; 
        end
    end
    
    % Average calculation
    average_errors = mean(all_participants_errors, 1);
    average_std = std(all_participants_errors, 0, 1);
    if shaded_plot == 0
        red_line = errorbar(1:3, average_errors, average_std, '--o', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0]); % Red average line
        bmh_list = [bmh_list; "Average"]; % Add average for the legend
        legend(red_line, 'Average', 'FontSize', 14); % Adjust the font size as needed
        bmh_list(strcmp(bmh_list, "Average")) = []; % Remove average to bmh_list
    else
        h_avg = plot(1:3, average_errors, '--o', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0], 'DisplayName', 'Average');
        legend_handles = [legend_handles, h_avg];
        legend(legend_handles, 'FontSize', 14);
    end

    set(gca, 'FontSize', 22);
    title("Walking Speed for all Participants for High Accuracy Condition")
    xlabel('Visual Perturbation');
    ylabel('Walking Speed (m/s)');
    xticks(1:3);
    xticklabels(balance_levels);
    grid on;
    hold off

    %% Stepping errors for each label (across all participants) - PLOT 19

    list_labels = fieldnames(data_labels); % labels stored in the struct

    figure
    hold on
    
    for k = 1:length(list_labels)
        current_label = string(list_labels{k});
        current_mean = mean(data_labels.(current_label).step_error); 
        current_std = std(data_labels.(current_label).step_error); 
        
        bar(k, current_mean, 'FaceColor', [0, 0.5, 1], 'EdgeColor', 'k');
        errorbar(k, current_mean, current_std, 'k', 'LineWidth', 1.5);
    end
    
    xticks(1:length(list_labels));
    xticklabels(list_labels); 
    xlabel('Conditions');
    ylabel('Mean Step Error (cm)');
    title('Step Error for Each Label Condition');
    set(gca, 'FontSize', 22);
    hold off;

    %% Plots 2,3,5,6 together - Plot 20

    avg_step_length = average_step_length_bybalance; 
    std_step_length = std_step_length_bybalance;
    
    avg_other_data = [average_step_length_variability_bybalance; average_step_width_bybalance; average_step_width_variability_bybalance];
    std_other_data = [std_step_length_variability_bybalance; std_step_width_bybalance; std_step_width_variability_bybalance];
    
    figure;
    
    % Adjust layout for the left plot to be 1/3 smaller
    subplot('Position', [0.05, 0.2, 0.25, 0.6]); % [left, bottom, width, height]
    bar_handle_left = bar(avg_step_length); % Create bars for step length
    hold on;
    
    num_groups_left = size(avg_step_length, 1); % Number of conditions
    num_bars_left = size(avg_step_length, 2); % Number of levels
    group_width_left = min(0.8, num_bars_left / (num_bars_left + 1.5));
    
    colors = lines(num_bars_left); % Use distinct colors
    
    for i = 1:num_bars_left
        b = bar(i, avg_step_length(:, i), 'FaceColor', colors(i, :)); % Set color for each bar
        errorbar(i, avg_step_length(:, i), std_step_length(:, i), 'k.', 'LineWidth', 1.2); % Add std error bars
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:num_groups_left);
    xticklabels({}); % Remove x-axis labels for the left plot
    xlabel('Step Length');
    ylabel('Distance (cm)');
    title('Step Length by Balance Condition');
    grid on;
    hold off;
    
    % Adjust layout for the right plot (2/3 larger)
    subplot('Position', [0.35, 0.2, 0.6, 0.6]); % [left, bottom, width, height]
    bar_handle_right = bar(avg_other_data); % Create bars for other variables
    hold on;
    
    num_groups_right = size(avg_other_data, 1); % Number of conditions
    num_bars_right = size(avg_other_data, 2); % Number of levels
    group_width_right = min(0.8, num_bars_right / (num_bars_right + 1.5));
    
    for i = 1:num_bars_right
        x = (1:num_groups_right) - group_width_right / 2 + (2 * i - 1) * group_width_right / (2 * num_bars_right);
        errorbar(x, avg_other_data(:, i), std_other_data(:, i), 'k.', 'LineWidth', 1.2);
    end
    
    for i = 1:numel(bar_handle_right)
        bar_handle_right(i).FaceColor = colors(i, :); % Use the same colors as left plot
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:num_groups_right);
    xticklabels({'Step Length Variability', 'Step Width', 'Step Width Variability'});
    ylabel('Distance (cm)');
    title('Other Variables by Balance Condition');
    grid on;
    hold off;
    
    legend({'No Perturbation', 'Medium Perturbation', 'High Perturbation'});

    %% Accuracy vs step width and step width variability - PLOT 21,22

    accuracy_0_step_width = [];
    accuracy_1_step_width = [];
    accuracy_2_step_width = [];
    
    accuracy_0_variability = [];
    accuracy_1_variability = [];
    accuracy_2_variability = [];
    
    for i=1:length(fieldnames(imported_data))
        bmh = bmh_list(i);
        for j=1:size(imported_data.(bmh).prompt.accuracy,1)
            if j>=3 && j<=29
                if double(imported_data.(bmh).prompt.accuracy(j,2)) == 0
                    accuracy_0_step_width = [accuracy_0_step_width, imported_data.(bmh).meanwidthstraights(1,j)];
                    accuracy_0_variability = [accuracy_0_variability, imported_data.(bmh).straightsvariability(1,j)];
                elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 1
                    accuracy_1_step_width = [accuracy_1_step_width, imported_data.(bmh).meanwidthstraights(1,j)];    
                    accuracy_1_variability = [accuracy_1_variability, imported_data.(bmh).straightsvariability(1,j)];
                elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 2
                    accuracy_2_step_width = [accuracy_2_step_width, imported_data.(bmh).meanwidthstraights(1,j)];
                    accuracy_2_variability = [accuracy_2_variability, imported_data.(bmh).straightsvariability(1,j)];
                end
            end
        end
    end
    
    % Step width and its std (which is the std across all step widths). 
    avg_accuracy_0_step_width = mean(accuracy_0_step_width)/10; % From mm to cm
    avg_accuracy_1_step_width = mean(accuracy_1_step_width)/10; % From mm to cm
    avg_accuracy_2_step_width = mean(accuracy_2_step_width)/10; % From mm to cm
    % std_accuracy_0_step_width = std(accuracy_0_step_width)/10; % From mm to cm
    % std_accuracy_1_step_width = std(accuracy_1_step_width)/10; % From mm to cm
    % std_accuracy_2_step_width = std(accuracy_2_step_width)/10; % From mm to cm
    std_accuracy_0_step_width = std(mean(reshape(accuracy_0_step_width, 9, []), 2), 'omitnan')/10; %do first the mean of the stds each 9 positions, which correspond to the same condition for all participants
    std_accuracy_1_step_width = std(mean(reshape(accuracy_1_step_width, 9, []), 2), 'omitnan')/10;
    std_accuracy_2_step_width = std(mean(reshape(accuracy_2_step_width, 9, []), 2), 'omitnan')/10;
    
    average_step_width_byaccuracy = [avg_accuracy_0_step_width, avg_accuracy_1_step_width, avg_accuracy_2_step_width];
    std_step_width_byaccuracy = [std_accuracy_0_step_width, std_accuracy_1_step_width, std_accuracy_2_step_width];
    
    figure;
    b = bar(average_step_width_byaccuracy, 'FaceColor', 'flat');
    hold on;
    
    for k = 1:length(average_speeds)
        b.CData(k, :) = colors(k, :);
    end
    
    errorbar(1:3, average_step_width_byaccuracy, std_step_width_byaccuracy, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    
    set(gca, 'xticklabel', {'None', 'Medium', 'High'});
    xlabel('Accuracy Prompt');
    ylabel('Average Step Width (cm)');
    title('Average Step Width by Accuracy Prompt - ' + bmh_text)
    ax = gca;               % Get current axis
    ax.FontSize = 14;       % Set axis number font size (e.g., 12);
    hold off;
    
    % Step width variability, which is the std of each single trial
    avg_accuracy_0_variability = mean(accuracy_0_variability)/10; % From mm to cm
    avg_accuracy_1_variability = mean(accuracy_1_variability)/10; % From mm to cm
    avg_accuracy_2_variability = mean(accuracy_2_variability)/10; % From mm to cm
    % std_accuracy_0_variability = std(accuracy_0_variability)/10; % From mm to cm
    % std_accuracy_1_variability = std(accuracy_1_variability)/10; % From mm to cm
    % std_accuracy_2_variability = std(accuracy_2_variability)/10; % From mm to cm
    std_accuracy_0_variability = std(mean(reshape(accuracy_0_variability, 9, []), 2), 'omitnan')/10; %do first the mean of the stds each 9 positions, which correspond to the same condition for all participants
    std_accuracy_1_variability = std(mean(reshape(accuracy_1_variability, 9, []), 2), 'omitnan')/10;
    std_accuracy_2_variability = std(mean(reshape(accuracy_2_variability, 9, []), 2), 'omitnan')/10;
    
    average_step_width_variability_byaccuracy = [avg_accuracy_0_variability, avg_accuracy_1_variability, avg_accuracy_2_variability];
    std_step_width_variability_byaccuracy = [std_accuracy_0_variability, std_accuracy_1_variability, std_accuracy_2_variability];
    
    figure;
    b = bar(average_step_width_variability_byaccuracy, 'FaceColor', 'flat');
    hold on;
    
    for k = 1:length(average_speeds)
        b.CData(k, :) = colors(k, :);
    end
    
    errorbar(1:3, average_step_width_variability_byaccuracy, std_step_width_variability_byaccuracy, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    
    set(gca, 'xticklabel', {'None', 'Medium', 'High'});
    xlabel('Accuracy Prompt');
    ylabel('Average Step Width Variability (cm)');
    title('Average Step Width Variability by Accuracy Prompt - ' + bmh_text);
    ax = gca;               % Get current axis
    ax.FontSize = 14;       % Set axis number font size (e.g., 12)
    hold off;

    %% Step length, width, and variability vs accuracy - PLOT 23

    avg_step_length = average_step_length_byaccuracy;
    std_step_length = std_step_length_byaccuracy;
    
    avg_other_data = [average_step_length_variability_byaccuracy; average_step_width_byaccuracy; average_step_width_variability_byaccuracy]; 
    std_other_data = [std_step_length_variability_byaccuracy; std_step_width_byaccuracy; std_step_width_variability_byaccuracy]; 
    
    figure;
    
    % Adjust layout for the left plot to be 1/3 smaller
    subplot('Position', [0.05, 0.2, 0.25, 0.6]); % [left, bottom, width, height]
    bar_handle_left = bar(avg_step_length); % Create bars for step length
    hold on;
    
    num_groups_left = size(avg_step_length, 1); % Number of conditions
    num_bars_left = size(avg_step_length, 2); % Number of levels
    group_width_left = min(0.8, num_bars_left / (num_bars_left + 1.5));
    
    colors = lines(num_bars_left); % Use distinct colors
    
    for i = 1:num_bars_left
        b = bar(i, avg_step_length(:, i), 'FaceColor', colors(i, :)); % Set color for each bar
        errorbar(i, avg_step_length(:, i), std_step_length(:, i), 'k.', 'LineWidth', 1.2); % Add std error bars
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:num_groups_left);
    xticklabels({}); % Remove x-axis labels for the left plot
    xlabel('Step Length');
    ylabel('Distance (cm)');
    title('Step Length by Accuracy Condition');
    grid on;
    hold off;
    
    % Adjust layout for the right plot (2/3 larger)
    subplot('Position', [0.35, 0.2, 0.6, 0.6]); % [left, bottom, width, height]
    bar_handle_right = bar(avg_other_data); % Create bars for other variables
    hold on;
    
    num_groups_right = size(avg_other_data, 1); % Number of conditions
    num_bars_right = size(avg_other_data, 2); % Number of levels
    group_width_right = min(0.8, num_bars_right / (num_bars_right + 1.5));
    
    for i = 1:num_bars_right
        x = (1:num_groups_right) - group_width_right / 2 + (2 * i - 1) * group_width_right / (2 * num_bars_right);
        errorbar(x, avg_other_data(:, i), std_other_data(:, i), 'k.', 'LineWidth', 1.2);
    end
    
    for i = 1:numel(bar_handle_right)
        bar_handle_right(i).FaceColor = colors(i, :); % Use the same colors as left plot
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:num_groups_right);
    xticklabels({'Step Length Variability', 'Step Width', 'Step Width Variability'}); % X-axis labels for right plot
    ylabel('Distance (cm)');
    title('Other Variables by Accuracy Condition');
    grid on;
    hold off;
    
    legend({'Null Accuracy', 'Medium Accuracy', 'High Accuracy'});

    %% Speed vs step length and step length variability - PLOT 24, 25

    speed_0_step_length = [];
    speed_1_step_length = [];
    speed_2_step_length = [];
    
    speed_0_variability = [];
    speed_1_variability = [];
    speed_2_variability = [];
    
    for i=1:length(fieldnames(imported_data))
        bmh = bmh_list(i);
        for j=1:size(imported_data.(bmh).prompt.speed,1)
            if j>=3 && j<=29
                if double(imported_data.(bmh).prompt.speed(j,2)) == 0
                    speed_0_step_length = [speed_0_step_length, imported_data.(bmh).meanlength_straights(1,j)];
                    speed_0_variability = [speed_0_variability, imported_data.(bmh).meanlength_straights_variability(1,j)];
                elseif double(imported_data.(bmh).prompt.speed(j,2)) == 1
                    speed_1_step_length = [speed_1_step_length, imported_data.(bmh).meanlength_straights(1,j)];    
                    speed_1_variability = [speed_1_variability, imported_data.(bmh).meanlength_straights_variability(1,j)];
                elseif double(imported_data.(bmh).prompt.speed(j,2)) == 2
                    speed_2_step_length = [speed_2_step_length, imported_data.(bmh).meanlength_straights(1,j)];
                    speed_2_variability = [speed_2_variability, imported_data.(bmh).meanlength_straights_variability(1,j)];
                end
            end
        end
    end
    
    % Step length and its std (which is the std across all step lengths). 
    avg_speed_0_step_length = mean(speed_0_step_length)/10; % From mm to cm
    avg_speed_1_step_length = mean(speed_1_step_length)/10; % From mm to cm
    avg_speed_2_step_length = mean(speed_2_step_length)/10; % From mm to cm
    % std_speed_0_step_length = std(speed_0_step_length)/10; % From mm to cm
    % std_speed_1_step_length = std(speed_1_step_length)/10; % From mm to cm
    % std_speed_2_step_length = std(speed_2_step_length)/10; % From mm to cm
    std_speed_0_step_length = std(mean(reshape(speed_0_step_length, 9, []), 2), 'omitnan')/10; %do first the mean of the stds each 9 positions, which correspond to the same condition for all participants
    std_speed_1_step_length = std(mean(reshape(speed_1_step_length, 9, []), 2), 'omitnan')/10;
    std_speed_2_step_length = std(mean(reshape(speed_2_step_length, 9, []), 2), 'omitnan')/10;
    
    average_step_length_byspeed = [avg_speed_0_step_length, avg_speed_1_step_length, avg_speed_2_step_length];
    std_step_length_byspeed = [std_speed_0_step_length, std_speed_1_step_length, std_speed_2_step_length];
    
    figure;
    b = bar(average_step_length_byspeed, 'FaceColor', 'flat');
    hold on;
    
    for k = 1:length(average_speeds)
        b.CData(k, :) = colors(k, :);
    end
    
    errorbar(1:3, average_step_length_byspeed, std_step_length_byspeed, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    
    set(gca, 'xticklabel', {'Slow', 'Medium', 'Fast'});
    xlabel('Speed Prompt');
    ylabel('Average Step Length (cm)');
    title('Average Step Length by Speed Prompt - ' + bmh_text)
    ax = gca;               % Get current axis
    ax.FontSize = 14;       % Set axis number font size (e.g., 12);
    hold off;
    
    % Step width variability, which is the std of each single trial
    avg_speed_0_variability = mean(speed_0_variability)/10; % From mm to cm
    avg_speed_1_variability = mean(speed_1_variability)/10; % From mm to cm
    avg_speed_2_variability = mean(speed_2_variability)/10; % From mm to cm
    % std_speed_0_variability = std(speed_0_variability)/10; % From mm to cm
    % std_speed_1_variability = std(speed_1_variability)/10; % From mm to cm
    % std_speed_2_variability = std(speed_2_variability)/10; % From mm to cm
    std_speed_0_variability = std(mean(reshape(speed_0_variability, 9, []), 2), 'omitnan')/10; %do first the mean of the stds each 9 positions, which correspond to the same condition for all participants
    std_speed_1_variability = std(mean(reshape(speed_1_variability, 9, []), 2), 'omitnan')/10;
    std_speed_2_variability = std(mean(reshape(speed_2_variability, 9, []), 2), 'omitnan')/10;
    
    average_step_length_variability_byspeed = [avg_speed_0_variability, avg_speed_1_variability, avg_speed_2_variability];
    std_step_length_variability_byspeed = [std_speed_0_variability, std_speed_1_variability, std_speed_2_variability];
    
    figure;
    b = bar(average_step_length_variability_byspeed, 'FaceColor', 'flat');
    hold on;
    
    for k = 1:length(average_speeds)
        b.CData(k, :) = colors(k, :);
    end
    
    errorbar(1:3, average_step_length_variability_byspeed, std_step_length_variability_byspeed, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    
    set(gca, 'xticklabel', {'None', 'Medium', 'High'});
    xlabel('Speed Prompt');
    ylabel('Average Step Length Variability (cm)');
    title('Average Step Length Variability by Speed Prompt - ' + bmh_text);
    ax = gca;               % Get current axis
    ax.FontSize = 14;       % Set axis number font size (e.g., 12)
    hold off;
    
    %% Speed vs step width and step width variability - PLOT 26, 27
    
    speed_0_step_width = [];
    speed_1_step_width = [];
    speed_2_step_width = [];
    
    speed_0_variability = [];
    speed_1_variability = [];
    speed_2_variability = [];
    
    for i=1:length(fieldnames(imported_data))
        bmh = bmh_list(i);
        for j=1:size(imported_data.(bmh).prompt.speed,1)
            if j>=3 && j<=29
                if double(imported_data.(bmh).prompt.speed(j,2)) == 0
                    speed_0_step_width = [speed_0_step_width, imported_data.(bmh).meanwidthstraights(1,j)];
                    speed_0_variability = [speed_0_variability, imported_data.(bmh).straightsvariability(1,j)];
                elseif double(imported_data.(bmh).prompt.speed(j,2)) == 1
                    speed_1_step_width = [speed_1_step_width, imported_data.(bmh).meanwidthstraights(1,j)];    
                    speed_1_variability = [speed_1_variability, imported_data.(bmh).straightsvariability(1,j)];
                elseif double(imported_data.(bmh).prompt.speed(j,2)) == 2
                    speed_2_step_width = [speed_2_step_width, imported_data.(bmh).meanwidthstraights(1,j)];
                    speed_2_variability = [speed_2_variability, imported_data.(bmh).straightsvariability(1,j)];
                end
            end
        end
    end
    
    % Step width and its std (which is the std across all step widths). 
    avg_speed_0_step_width = mean(speed_0_step_width)/10; % From mm to cm
    avg_speed_1_step_width = mean(speed_1_step_width)/10; % From mm to cm
    avg_speed_2_step_width = mean(speed_2_step_width)/10; % From mm to cm
    % std_speed_0_step_width = std(speed_0_step_width)/10; % From mm to cm
    % std_speed_1_step_width = std(speed_1_step_width)/10; % From mm to cm
    % std_speed_2_step_width = std(speed_2_step_width)/10; % From mm to cm
    std_speed_0_step_width = std(mean(reshape(speed_0_step_width, 9, []), 2), 'omitnan')/10; %do first the mean of the stds each 9 positions, which correspond to the same condition for all participants
    std_speed_1_step_width = std(mean(reshape(speed_1_step_width, 9, []), 2), 'omitnan')/10;
    std_speed_2_step_width = std(mean(reshape(speed_2_step_width, 9, []), 2), 'omitnan')/10;
    
    average_step_width_byspeed = [avg_speed_0_step_width, avg_speed_1_step_width, avg_speed_2_step_width];
    std_step_width_byspeed = [std_speed_0_step_width, std_speed_1_step_width, std_speed_2_step_width];
    
    figure;
    b = bar(average_step_width_byspeed, 'FaceColor', 'flat');
    hold on;
    
    for k = 1:length(average_speeds)
        b.CData(k, :) = colors(k, :);
    end
    
    errorbar(1:3, average_step_width_byspeed, std_step_width_byspeed, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    
    set(gca, 'xticklabel', {'None', 'Medium', 'High'});
    xlabel('Speed Prompt');
    ylabel('Average Step Width (cm)');
    title('Average Step Width by Speed Prompt - ' + bmh_text)
    ax = gca;               % Get current axis
    ax.FontSize = 14;       % Set axis number font size (e.g., 12);
    hold off;
    
    % Step width variability, which is the std of each single trial
    avg_speed_0_variability = mean(speed_0_variability)/10; % From mm to cm
    avg_speed_1_variability = mean(speed_1_variability)/10; % From mm to cm
    avg_speed_2_variability = mean(speed_2_variability)/10; % From mm to cm
    % std_speed_0_variability = std(speed_0_variability)/10; % From mm to cm
    % std_speed_1_variability = std(speed_1_variability)/10; % From mm to cm
    % std_speed_2_variability = std(speed_2_variability)/10; % From mm to cm
    std_speed_0_variability = std(mean(reshape(speed_0_variability, 9, []), 2), 'omitnan')/10; %do first the mean of the stds each 9 positions, which correspond to the same condition for all participants
    std_speed_1_variability = std(mean(reshape(speed_1_variability, 9, []), 2), 'omitnan')/10;
    std_speed_2_variability = std(mean(reshape(speed_2_variability, 9, []), 2), 'omitnan')/10;
    
    average_step_width_variability_byspeed = [avg_speed_0_variability, avg_speed_1_variability, avg_speed_2_variability];
    std_step_width_variability_byspeed = [std_speed_0_variability, std_speed_1_variability, std_speed_2_variability];
    
    figure;
    b = bar(average_step_width_variability_byspeed, 'FaceColor', 'flat');
    hold on;
    
    for k = 1:length(average_speeds)
        b.CData(k, :) = colors(k, :);
    end
    
    errorbar(1:3, average_step_width_variability_byspeed, std_step_width_variability_byspeed, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
    
    set(gca, 'xticklabel', {'Slow', 'Medium', 'Fast'});
    xlabel('Speed Prompt');
    ylabel('Average Step Width Variability (cm)');
    title('Average Step Width Variability by Accuracy Prompt - ' + bmh_text);
    ax = gca;               % Get current axis
    ax.FontSize = 14;       % Set axis number font size (e.g., 12)
    hold off;

    %% Step length, width, and variability vs speed - PLOT 28

    avg_step_length = average_step_length_byspeed; % Step Length
    std_step_length = std_step_length_byspeed;
    
    avg_other_data = [average_step_length_variability_byspeed; average_step_width_byspeed; average_step_width_variability_byspeed];
    std_other_data = [std_step_length_variability_byspeed; std_step_width_byspeed; std_step_width_variability_byspeed];
    
    figure;
    
    % Adjust layout for the left plot to be 1/3 smaller
    subplot('Position', [0.05, 0.2, 0.25, 0.6]); % [left, bottom, width, height]
    bar_handle_left = bar(avg_step_length); % Create bars for step length
    hold on;
    
    num_groups_left = size(avg_step_length, 1); % Number of conditions
    num_bars_left = size(avg_step_length, 2); % Number of levels
    group_width_left = min(0.8, num_bars_left / (num_bars_left + 1.5));
    
    colors = lines(num_bars_left); % Use distinct colors
    
    for i = 1:num_bars_left
        b = bar(i, avg_step_length(:, i), 'FaceColor', colors(i, :)); % Set color for each bar
        errorbar(i, avg_step_length(:, i), std_step_length(:, i), 'k.', 'LineWidth', 1.2); % Add std error bars
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:num_groups_left);
    xticklabels({}); % Remove x-axis labels for the left plot
    xlabel('Step Length');
    ylabel('Distance (cm)');
    title('Step Length by Speed Condition');
    grid on;
    hold off;
    
    % Adjust layout for the right plot (2/3 larger)
    subplot('Position', [0.35, 0.2, 0.6, 0.6]); % [left, bottom, width, height]
    bar_handle_right = bar(avg_other_data); % Create bars for other variables
    hold on;
    
    num_groups_right = size(avg_other_data, 1); % Number of conditions
    num_bars_right = size(avg_other_data, 2); % Number of levels
    group_width_right = min(0.8, num_bars_right / (num_bars_right + 1.5));
    
    for i = 1:num_bars_right
        x = (1:num_groups_right) - group_width_right / 2 + (2 * i - 1) * group_width_right / (2 * num_bars_right);
        errorbar(x, avg_other_data(:, i), std_other_data(:, i), 'k.', 'LineWidth', 1.2);
    end
    
    for i = 1:numel(bar_handle_right)
        bar_handle_right(i).FaceColor = colors(i, :); % Use the same colors as left plot
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:num_groups_right);
    xticklabels({'Step Length Variability', 'Step Width', 'Step Width Variability'}); % X-axis labels for right plot
    ylabel('Distance (cm)');
    title('Other Variables by Speed Condition');
    grid on;
    hold off;
    
    % Add a single legend for both plots, placed on the top-right
    legend({'Slow', 'Medium', 'Fast'});

    %% Walking speed by condition - PLOT 29

    speed_balance0 = [];
    speed_balance1 = [];
    speed_balance2 = [];
    speed_speed0 = [];
    speed_speed1 = [];
    speed_speed2 = [];
    speed_accuracy0 = [];
    speed_accuracy1 = [];
    speed_accuracy2 = [];

    for i=1:length(list_labels)
        current_label = list_labels{i};

        if contains(current_label, 's0')
            speed_speed0 = [speed_speed0, data_labels.(current_label).walking_speed(:,1)];
        end
        if contains(current_label, 's1')
            speed_speed1 = [speed_speed1, data_labels.(current_label).walking_speed(:,1)];
        end
        if contains(current_label, 's2')
            speed_speed2 = [speed_speed2, data_labels.(current_label).walking_speed(:,1)];
        end

        if contains(current_label, 'a0')
            speed_accuracy0 = [speed_accuracy0, data_labels.(current_label).walking_speed(:,1)];
        end
        if contains(current_label, 'a1')
            speed_accuracy1 = [speed_accuracy1, data_labels.(current_label).walking_speed(:,1)];
        end
        if contains(current_label, 'a2')
            speed_accuracy2 = [speed_accuracy2, data_labels.(current_label).walking_speed(:,1)];
        end

        if contains(current_label, 'b0')
            speed_balance0 = [speed_balance0, data_labels.(current_label).walking_speed(:,1)];
        end
        if contains(current_label, 'b1')
            speed_balance1 = [speed_balance1, data_labels.(current_label).walking_speed(:,1)];
        end
        if contains(current_label, 'b2')
            speed_balance2 = [speed_balance2, data_labels.(current_label).walking_speed(:,1)];
        end
    end

    speed_balance0_std = std(mean(speed_balance0,2));
    speed_balance1_std = std(mean(speed_balance1,2));
    speed_balance2_std = std(mean(speed_balance2,2));
    speed_speed0_std = std(mean(speed_speed0,2));
    speed_speed1_std = std(mean(speed_speed1,2));
    speed_speed2_std = std(mean(speed_speed2,2));
    speed_accuracy0_std = std(mean(speed_accuracy0,2));
    speed_accuracy1_std = std(mean(speed_accuracy1,2));
    speed_accuracy2_std = std(mean(speed_accuracy2,2));

    speed_balance0_avg = mean(speed_balance0(:));
    speed_balance1_avg = mean(speed_balance1(:));
    speed_balance2_avg = mean(speed_balance2(:));
    speed_speed0_avg = mean(speed_speed0(:));
    speed_speed1_avg = mean(speed_speed1(:));
    speed_speed2_avg = mean(speed_speed2(:));
    speed_accuracy0_avg = mean(speed_accuracy0(:));
    speed_accuracy1_avg = mean(speed_accuracy1(:));
    speed_accuracy2_avg = mean(speed_accuracy2(:));

    means = [speed_balance0_avg, speed_balance1_avg, speed_balance2_avg; 
             speed_speed0_avg, speed_speed1_avg, speed_speed2_avg; 
             speed_accuracy0_avg, speed_accuracy1_avg, speed_accuracy2_avg];
    
    stds = [speed_balance0_std, speed_balance1_std, speed_balance2_std; 
            speed_speed0_std, speed_speed1_std, speed_speed2_std; 
            speed_accuracy0_std, speed_accuracy1_std, speed_accuracy2_std];
        
    figure;
    hold on;
    categories = {'Visual Perturbation', 'Speed Prompt', 'Accuracy Prompt'};
    x = 1:3; % Groups
    bar_handles = bar(x, means, 'grouped'); 

    if plot_eachparticipant == 1
        points_means = {mean(speed_balance0,2), mean(speed_balance1,2), mean(speed_balance2,2); 
                    mean(speed_speed0,2), mean(speed_speed1,2), mean(speed_speed2,2); 
                    mean(speed_accuracy0,2), mean(speed_accuracy1,2), mean(speed_accuracy2,2)};
        
        num_trials = size(points_means{1,1}, 1); % Assuming all groups have the same number of trials
        %colors = lines(num_trials); % Use 'lines' colormap to generate distinct colors
        %colors = min(colors * 1.4, 1); % Increase intensity, limit to 1
        for i = 1:size(means, 2) % Loop through balance, speed, accuracy columns
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars
            
            % Scatter plot for individual points with unique colors per row
            for j = 1:3 % Loop through categories (Balance, Speed, Accuracy)
                data_points = points_means{j, i}; % Extract row-wise mean values
                for k = 1:num_trials % Loop through trials (rows)
                    scatter(x_offset(j), data_points(k), point_size, colors_participants(k, :), 'filled', 'MarkerFaceAlpha', 0.6);
                end
            end
        end
    else
        for i = 1:size(means, 2) 
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            errorbar(x_offset, means(:,i), stds(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1.5);
        end
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:3); 
    xticklabels(categories);
    %ylabel('Walking Speed (m/s)');
    ylabel('Mean Value (m/s)')
    legend({'Null/Low', 'Medium', 'High'});
    title('Walking Speed');
    grid on;
    hold off;

    %% Step length by condition - PLOT 30

    step_length_balance0 = [];
    step_length_balance1 = [];
    step_length_balance2 = [];
    step_length_speed0 = [];
    step_length_speed1 = [];
    step_length_speed2 = [];
    step_length_accuracy0 = [];
    step_length_accuracy1 = [];
    step_length_accuracy2 = [];

    for i=1:length(list_labels)
        current_label = list_labels{i};

        if contains(current_label, 's0')
            step_length_speed0 = [step_length_speed0, data_labels.(current_label).step_length(:,1)];
        end
        if contains(current_label, 's1')
            step_length_speed1 = [step_length_speed1, data_labels.(current_label).step_length(:,1)];
        end
        if contains(current_label, 's2')
            step_length_speed2 = [step_length_speed2, data_labels.(current_label).step_length(:,1)];
        end

        if contains(current_label, 'a0')
            step_length_accuracy0 = [step_length_accuracy0, data_labels.(current_label).step_length(:,1)];
        end
        if contains(current_label, 'a1')
            step_length_accuracy1 = [step_length_accuracy1, data_labels.(current_label).step_length(:,1)];
        end
        if contains(current_label, 'a2')
            step_length_accuracy2 = [step_length_accuracy2, data_labels.(current_label).step_length(:,1)];
        end

        if contains(current_label, 'b0')
            step_length_balance0 = [step_length_balance0, data_labels.(current_label).step_length(:,1)];
        end
        if contains(current_label, 'b1')
            step_length_balance1 = [step_length_balance1, data_labels.(current_label).step_length(:,1)];
        end
        if contains(current_label, 'b2')
            step_length_balance2 = [step_length_balance2, data_labels.(current_label).step_length(:,1)];
        end
    end

    step_length_balance0_std = std(mean(step_length_balance0,2));
    step_length_balance1_std = std(mean(step_length_balance1,2));
    step_length_balance2_std = std(mean(step_length_balance2,2));
    step_length_speed0_std = std(mean(step_length_speed0,2));
    step_length_speed1_std = std(mean(step_length_speed1,2));
    step_length_speed2_std = std(mean(step_length_speed2,2));
    step_length_accuracy0_std = std(mean(step_length_accuracy0,2));
    step_length_accuracy1_std = std(mean(step_length_accuracy1,2));
    step_length_accuracy2_std = std(mean(step_length_accuracy2,2));

    step_length_balance0_avg = mean(step_length_balance0(:));
    step_length_balance1_avg = mean(step_length_balance1(:));
    step_length_balance2_avg = mean(step_length_balance2(:));
    step_length_speed0_avg = mean(step_length_speed0(:));
    step_length_speed1_avg = mean(step_length_speed1(:));
    step_length_speed2_avg = mean(step_length_speed2(:));
    step_length_accuracy0_avg = mean(step_length_accuracy0(:));
    step_length_accuracy1_avg = mean(step_length_accuracy1(:));
    step_length_accuracy2_avg = mean(step_length_accuracy2(:));

    means = [step_length_balance0_avg, step_length_balance1_avg, step_length_balance2_avg; 
             step_length_speed0_avg, step_length_speed1_avg, step_length_speed2_avg; 
             step_length_accuracy0_avg, step_length_accuracy1_avg, step_length_accuracy2_avg];
    
    stds = [step_length_balance0_std, step_length_balance1_std, step_length_balance2_std; 
            step_length_speed0_std, step_length_speed1_std, step_length_speed2_std; 
            step_length_accuracy0_std, step_length_accuracy1_std, step_length_accuracy2_std];
    
    categories = {'Visual Perturbation', 'Speed Prompt', 'Accuracy Prompt'};
    x = 1:3; % Groups
    
    figure;
    hold on;
    bar_handles = bar(x, means, 'grouped'); % Bar plot

    if plot_eachparticipant == 1
        points_means = {mean(step_length_balance0,2), mean(step_length_balance1,2), mean(step_length_balance2,2); 
                    mean(step_length_speed0,2), mean(step_length_speed1,2), mean(step_length_speed2,2); 
                    mean(step_length_accuracy0,2), mean(step_length_accuracy1,2), mean(step_length_accuracy2,2)};
        
        num_trials = size(points_means{1,1}, 1); % Assuming all groups have the same number of trials
        %colors = lines(num_trials); % Use 'lines' colormap to generate distinct colors
        %colors = min(colors * 1.4, 1); % Increase intensity, limit to 1
        for i = 1:size(means, 2) % Loop through balance, speed, accuracy columns
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars
            
            % Scatter plot for individual points with unique colors per row
            for j = 1:3 % Loop through categories (Balance, Speed, Accuracy)
                data_points = points_means{j, i}; % Extract row-wise mean values
                for k = 1:num_trials % Loop through trials (rows)
                    scatter(x_offset(j), data_points(k), point_size, colors_participants(k, :), 'filled', 'MarkerFaceAlpha', 0.6);
                end
            end
        end
    else
        for i = 1:size(means, 2) % Loop through each group
            % Error bars
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            errorbar(x_offset, means(:,i), stds(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1.5);
        end
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:3); % Set x-ticks
    ylim([0 97])
    xticklabels(categories); % Label groups
    %ylabel('Average Step Length (cm)');
    ylabel('Mean Value (cm)')
    legend({'Null/Low', 'Medium', 'High'});
    title('Step Length');
    grid on;
    hold off;

    %% Step length variability by condition - PLOT 31

    step_length_var_balance0 = [];
    step_length_var_balance1 = [];
    step_length_var_balance2 = [];
    step_length_var_speed0 = [];
    step_length_var_speed1 = [];
    step_length_var_speed2 = [];
    step_length_var_accuracy0 = [];
    step_length_var_accuracy1 = [];
    step_length_var_accuracy2 = [];

    for i=1:length(list_labels)
        current_label = list_labels{i};

        if contains(current_label, 's0')
            step_length_var_speed0 = [step_length_var_speed0, data_labels.(current_label).step_length_var(:,1)];
        end
        if contains(current_label, 's1')
            step_length_var_speed1 = [step_length_var_speed1, data_labels.(current_label).step_length_var(:,1)];
        end
        if contains(current_label, 's2')
            step_length_var_speed2 = [step_length_var_speed2, data_labels.(current_label).step_length_var(:,1)];
        end

        if contains(current_label, 'a0')
            step_length_var_accuracy0 = [step_length_var_accuracy0, data_labels.(current_label).step_length_var(:,1)];
        end
        if contains(current_label, 'a1')
            step_length_var_accuracy1 = [step_length_var_accuracy1, data_labels.(current_label).step_length_var(:,1)];
        end
        if contains(current_label, 'a2')
            step_length_var_accuracy2 = [step_length_var_accuracy2, data_labels.(current_label).step_length_var(:,1)];
        end

        if contains(current_label, 'b0')
            step_length_var_balance0 = [step_length_var_balance0, data_labels.(current_label).step_length_var(:,1)];
        end
        if contains(current_label, 'b1')
            step_length_var_balance1 = [step_length_var_balance1, data_labels.(current_label).step_length_var(:,1)];
        end
        if contains(current_label, 'b2')
            step_length_var_balance2 = [step_length_var_balance2, data_labels.(current_label).step_length_var(:,1)];
        end
    end

    step_length_var_balance0_std = std(mean(step_length_var_balance0,2));
    step_length_var_balance1_std = std(mean(step_length_var_balance1,2));
    step_length_var_balance2_std = std(mean(step_length_var_balance2,2));
    step_length_var_speed0_std = std(mean(step_length_var_speed0,2));
    step_length_var_speed1_std = std(mean(step_length_var_speed1,2));
    step_length_var_speed2_std = std(mean(step_length_var_speed2,2));
    step_length_var_accuracy0_std = std(mean(step_length_var_accuracy0,2));
    step_length_var_accuracy1_std = std(mean(step_length_var_accuracy1,2));
    step_length_var_accuracy2_std = std(mean(step_length_var_accuracy2,2));

    step_length_var_balance0_avg = mean(step_length_var_balance0(:));
    step_length_var_balance1_avg = mean(step_length_var_balance1(:));
    step_length_var_balance2_avg = mean(step_length_var_balance2(:));
    step_length_var_speed0_avg = mean(step_length_var_speed0(:));
    step_length_var_speed1_avg = mean(step_length_var_speed1(:));
    step_length_var_speed2_avg = mean(step_length_var_speed2(:));
    step_length_var_accuracy0_avg = mean(step_length_var_accuracy0(:));
    step_length_var_accuracy1_avg = mean(step_length_var_accuracy1(:));
    step_length_var_accuracy2_avg = mean(step_length_var_accuracy2(:));

    means = [step_length_var_balance0_avg, step_length_var_balance1_avg, step_length_var_balance2_avg; 
             step_length_var_speed0_avg, step_length_var_speed1_avg, step_length_var_speed2_avg; 
             step_length_var_accuracy0_avg, step_length_var_accuracy1_avg, step_length_var_accuracy2_avg];
    
    stds = [step_length_var_balance0_std, step_length_var_balance1_std, step_length_var_balance2_std; 
            step_length_var_speed0_std, step_length_var_speed1_std, step_length_var_speed2_std; 
            step_length_var_accuracy0_std, step_length_var_accuracy1_std, step_length_var_accuracy2_std];

    categories = {'Visual Perturbation', 'Speed Prompt', 'Accuracy Prompt'};
    x = 1:3; % Groups
    
    figure;
    hold on;
    bar_handles = bar(x, means, 'grouped'); % Bar plot

    if plot_eachparticipant == 1
        points_means = {mean(step_length_var_balance0,2), mean(step_length_var_balance1,2), mean(step_length_var_balance2,2); 
                    mean(step_length_var_speed0,2), mean(step_length_var_speed1,2), mean(step_length_var_speed2,2); 
                    mean(step_length_var_accuracy0,2), mean(step_length_var_accuracy1,2), mean(step_length_var_accuracy2,2)};
        
        num_trials = size(points_means{1,1}, 1); % Assuming all groups have the same number of trials
        %colors = lines(num_trials); % Use 'lines' colormap to generate distinct colors
        %colors = min(colors * 1.4, 1); % Increase intensity, limit to 1
        for i = 1:size(means, 2) % Loop through balance, speed, accuracy columns
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars
            
            % Scatter plot for individual points with unique colors per row
            for j = 1:3 % Loop through categories (Balance, Speed, Accuracy)
                data_points = points_means{j, i}; % Extract row-wise mean values
                for k = 1:num_trials % Loop through trials (rows)
                    scatter(x_offset(j), data_points(k), point_size, colors_participants(k, :), 'filled', 'MarkerFaceAlpha', 0.6);
                end
            end
        end
    else
        for i = 1:size(means, 2) % Loop through each group
            % Error bars
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            errorbar(x_offset, means(:,i), stds(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1.5);
        end
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:3); % Set x-ticks
    ylim([0 21])
    xticklabels(categories); % Label groups
    %ylabel('Step Length Variability (cm)');
    ylabel('Mean Value (cm)')
    legend({'Null/Low', 'Medium', 'High'});
    title('Step Length Variability');
    grid on;
    hold off;

    %% Step width by condition - PLOT 32

    step_width_balance0 = [];
    step_width_balance1 = [];
    step_width_balance2 = [];
    step_width_speed0 = [];
    step_width_speed1 = [];
    step_width_speed2 = [];
    step_width_accuracy0 = [];
    step_width_accuracy1 = [];
    step_width_accuracy2 = [];

    for i=1:length(list_labels)
        current_label = list_labels{i};

        if contains(current_label, 's0')
            step_width_speed0 = [step_width_speed0, data_labels.(current_label).step_width(:,1)];
        end
        if contains(current_label, 's1')
            step_width_speed1 = [step_width_speed1, data_labels.(current_label).step_width(:,1)];
        end
        if contains(current_label, 's2')
            step_width_speed2 = [step_width_speed2, data_labels.(current_label).step_width(:,1)];
        end

        if contains(current_label, 'a0')
            step_width_accuracy0 = [step_width_accuracy0, data_labels.(current_label).step_width(:,1)];
        end
        if contains(current_label, 'a1')
            step_width_accuracy1 = [step_width_accuracy1, data_labels.(current_label).step_width(:,1)];
        end
        if contains(current_label, 'a2')
            step_width_accuracy2 = [step_width_accuracy2, data_labels.(current_label).step_width(:,1)];
        end

        if contains(current_label, 'b0')
            step_width_balance0 = [step_width_balance0, data_labels.(current_label).step_width(:,1)];
        end
        if contains(current_label, 'b1')
            step_width_balance1 = [step_width_balance1, data_labels.(current_label).step_width(:,1)];
        end
        if contains(current_label, 'b2')
            step_width_balance2 = [step_width_balance2, data_labels.(current_label).step_width(:,1)];
        end
    end

    step_width_balance0_std = std(mean(step_width_balance0,2));
    step_width_balance1_std = std(mean(step_width_balance1,2));
    step_width_balance2_std = std(mean(step_width_balance2,2));
    step_width_speed0_std = std(mean(step_width_speed0,2));
    step_width_speed1_std = std(mean(step_width_speed1,2));
    step_width_speed2_std = std(mean(step_width_speed2,2));
    step_width_accuracy0_std = std(mean(step_width_accuracy0,2));
    step_width_accuracy1_std = std(mean(step_width_accuracy1,2));
    step_width_accuracy2_std = std(mean(step_width_accuracy2,2));

    step_width_balance0_avg = mean(step_width_balance0(:));
    step_width_balance1_avg = mean(step_width_balance1(:));
    step_width_balance2_avg = mean(step_width_balance2(:));
    step_width_speed0_avg = mean(step_width_speed0(:));
    step_width_speed1_avg = mean(step_width_speed1(:));
    step_width_speed2_avg = mean(step_width_speed2(:));
    step_width_accuracy0_avg = mean(step_width_accuracy0(:));
    step_width_accuracy1_avg = mean(step_width_accuracy1(:));
    step_width_accuracy2_avg = mean(step_width_accuracy2(:));

    means = [step_width_balance0_avg, step_width_balance1_avg, step_width_balance2_avg; 
             step_width_speed0_avg, step_width_speed1_avg, step_width_speed2_avg; 
             step_width_accuracy0_avg, step_width_accuracy1_avg, step_width_accuracy2_avg];
    
    stds = [step_width_balance0_std, step_width_balance1_std, step_width_balance2_std; 
            step_width_speed0_std, step_width_speed1_std, step_width_speed2_std; 
            step_width_accuracy0_std, step_width_accuracy1_std, step_width_accuracy2_std];

    categories = {'Visual Perturbation', 'Speed Prompt', 'Accuracy Prompt'};
    x = 1:3; % Groups
    
    figure;
    hold on;
    bar_handles = bar(x, means, 'grouped'); % Bar plot

    if plot_eachparticipant == 1
        points_means = {mean(step_width_balance0,2), mean(step_width_balance1,2), mean(step_width_balance2,2); 
                    mean(step_width_speed0,2), mean(step_width_speed1,2), mean(step_width_speed2,2); 
                    mean(step_width_accuracy0,2), mean(step_width_accuracy1,2), mean(step_width_accuracy2,2)};
        
        num_trials = size(points_means{1,1}, 1); % Assuming all groups have the same number of trials
        %colors = lines(num_trials); % Use 'lines' colormap to generate distinct colors
        %colors = min(colors * 1.4, 1); % Increase intensity, limit to 1
        for i = 1:size(means, 2) % Loop through balance, speed, accuracy columns
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars
            
            % Scatter plot for individual points with unique colors per row
            for j = 1:3 % Loop through categories (Balance, Speed, Accuracy)
                data_points = points_means{j, i}; % Extract row-wise mean values
                for k = 1:num_trials % Loop through trials (rows)
                    scatter(x_offset(j), data_points(k), point_size, colors_participants(k, :), 'filled', 'MarkerFaceAlpha', 0.6);
                end
            end
        end
    else
        for i = 1:size(means, 2) % Loop through each group
            % Error bars
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            errorbar(x_offset, means(:,i), stds(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1.5);
        end
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:3); % Set x-ticks
    ylim([0 25])
    xticklabels(categories); % Label groups
    %ylabel('Step Width (cm)');
    ylabel('Mean Value (cm)')
    legend({'Null/Low', 'Medium', 'High'});
    title('Step Width');
    grid on;
    hold off;

    %% Step width variability by condition - PLOT 33

    step_width_var_balance0 = [];
    step_width_var_balance1 = [];
    step_width_var_balance2 = [];
    step_width_var_speed0 = [];
    step_width_var_speed1 = [];
    step_width_var_speed2 = [];
    step_width_var_accuracy0 = [];
    step_width_var_accuracy1 = [];
    step_width_var_accuracy2 = [];

    for i=1:length(list_labels)
        current_label = list_labels{i};

        if contains(current_label, 's0')
            step_width_var_speed0 = [step_width_var_speed0, data_labels.(current_label).step_width_var(:,1)];
        end
        if contains(current_label, 's1')
            step_width_var_speed1 = [step_width_var_speed1, data_labels.(current_label).step_width_var(:,1)];
        end
        if contains(current_label, 's2')
            step_width_var_speed2 = [step_width_var_speed2, data_labels.(current_label).step_width_var(:,1)];
        end

        if contains(current_label, 'a0')
            step_width_var_accuracy0 = [step_width_var_accuracy0, data_labels.(current_label).step_width_var(:,1)];
        end
        if contains(current_label, 'a1')
            step_width_var_accuracy1 = [step_width_var_accuracy1, data_labels.(current_label).step_width_var(:,1)];
        end
        if contains(current_label, 'a2')
            step_width_var_accuracy2 = [step_width_var_accuracy2, data_labels.(current_label).step_width_var(:,1)];
        end

        if contains(current_label, 'b0')
            step_width_var_balance0 = [step_width_var_balance0, data_labels.(current_label).step_width_var(:,1)];
        end
        if contains(current_label, 'b1')
            step_width_var_balance1 = [step_width_var_balance1, data_labels.(current_label).step_width_var(:,1)];
        end
        if contains(current_label, 'b2')
            step_width_var_balance2 = [step_width_var_balance2, data_labels.(current_label).step_width_var(:,1)];
        end
    end

    step_width_var_balance0_std = std(mean(step_width_var_balance0,2));
    step_width_var_balance1_std = std(mean(step_width_var_balance1,2));
    step_width_var_balance2_std = std(mean(step_width_var_balance2,2));
    step_width_var_speed0_std = std(mean(step_width_var_speed0,2));
    step_width_var_speed1_std = std(mean(step_width_var_speed1,2));
    step_width_var_speed2_std = std(mean(step_width_var_speed2,2));
    step_width_var_accuracy0_std = std(mean(step_width_var_accuracy0,2));
    step_width_var_accuracy1_std = std(mean(step_width_var_accuracy1,2));
    step_width_var_accuracy2_std = std(mean(step_width_var_accuracy2,2));

    step_width_var_balance0_avg = mean(step_width_var_balance0(:));
    step_width_var_balance1_avg = mean(step_width_var_balance1(:));
    step_width_var_balance2_avg = mean(step_width_var_balance2(:));
    step_width_var_speed0_avg = mean(step_width_var_speed0(:));
    step_width_var_speed1_avg = mean(step_width_var_speed1(:));
    step_width_var_speed2_avg = mean(step_width_var_speed2(:));
    step_width_var_accuracy0_avg = mean(step_width_var_accuracy0(:));
    step_width_var_accuracy1_avg = mean(step_width_var_accuracy1(:));
    step_width_var_accuracy2_avg = mean(step_width_var_accuracy2(:));

    means = [step_width_var_balance0_avg, step_width_var_balance1_avg, step_width_var_balance2_avg; 
             step_width_var_speed0_avg, step_width_var_speed1_avg, step_width_var_speed2_avg; 
             step_width_var_accuracy0_avg, step_width_var_accuracy1_avg, step_width_var_accuracy2_avg];
    
    stds = [step_width_var_balance0_std, step_width_var_balance1_std, step_width_var_balance2_std; 
            step_width_var_speed0_std, step_width_var_speed1_std, step_width_var_speed2_std; 
            step_width_var_accuracy0_std, step_width_var_accuracy1_std, step_width_var_accuracy2_std];

    categories = {'Visual Perturbation', 'Speed Prompt', 'Accuracy Prompt'};
    x = 1:3; % Groups
    
    figure;
    hold on;
    bar_handles = bar(x, means, 'grouped'); % Bar plot

    if plot_eachparticipant == 1
        points_means = {mean(step_width_var_balance0,2), mean(step_width_var_balance1,2), mean(step_width_var_balance2,2); 
                    mean(step_width_var_speed0,2), mean(step_width_var_speed1,2), mean(step_width_var_speed2,2); 
                    mean(step_width_var_accuracy0,2), mean(step_width_var_accuracy1,2), mean(step_width_var_accuracy2,2)};
        
        num_trials = size(points_means{1,1}, 1); % Assuming all groups have the same number of trials
        %colors = lines(num_trials); % Use 'lines' colormap to generate distinct colors
        %colors = min(colors * 1.4, 1); % Increase intensity, limit to 1
        for i = 1:size(means, 2) % Loop through balance, speed, accuracy columns
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars
            
            % Scatter plot for individual points with unique colors per row
            for j = 1:3 % Loop through categories (Balance, Speed, Accuracy)
                data_points = points_means{j, i}; % Extract row-wise mean values
                for k = 1:num_trials % Loop through trials (rows)
                    scatter(x_offset(j), data_points(k), point_size, colors_participants(k, :), 'filled', 'MarkerFaceAlpha', 0.6);
                end
            end
        end
    else
        for i = 1:size(means, 2) % Loop through each group
            % Error bars
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            errorbar(x_offset, means(:,i), stds(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1.5);
        end
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:3); % Set x-ticks
    xticklabels(categories); % Label groups
    %ylabel('Step Width Variability (cm)');
    ylabel('Mean Value (cm)')
    legend({'Null/Low', 'Medium', 'High'});
    title('Step Width Variability');
    grid on;
    hold off;

    %% Head angle by condition - PLOT 34

    head_angle_balance0 = [];
    head_angle_balance1 = [];
    head_angle_balance2 = [];
    head_angle_speed0 = [];
    head_angle_speed1 = [];
    head_angle_speed2 = [];
    head_angle_accuracy0 = [];
    head_angle_accuracy1 = [];
    head_angle_accuracy2 = [];

    for i=1:length(list_labels)
        current_label = list_labels{i};

        if contains(current_label, 's0')
            head_angle_speed0 = [head_angle_speed0, data_labels.(current_label).head_angle(:,1)];
        end
        if contains(current_label, 's1')
            head_angle_speed1 = [head_angle_speed1, data_labels.(current_label).head_angle(:,1)];
        end
        if contains(current_label, 's2')
            head_angle_speed2 = [head_angle_speed2, data_labels.(current_label).head_angle(:,1)];
        end

        if contains(current_label, 'a0')
            head_angle_accuracy0 = [head_angle_accuracy0, data_labels.(current_label).head_angle(:,1)];
        end
        if contains(current_label, 'a1')
            head_angle_accuracy1 = [head_angle_accuracy1, data_labels.(current_label).head_angle(:,1)];
        end
        if contains(current_label, 'a2')
            head_angle_accuracy2 = [head_angle_accuracy2, data_labels.(current_label).head_angle(:,1)];
        end

        if contains(current_label, 'b0')
            head_angle_balance0 = [head_angle_balance0, data_labels.(current_label).head_angle(:,1)];
        end
        if contains(current_label, 'b1')
            head_angle_balance1 = [head_angle_balance1, data_labels.(current_label).head_angle(:,1)];
        end
        if contains(current_label, 'b2')
            head_angle_balance2 = [head_angle_balance2, data_labels.(current_label).head_angle(:,1)];
        end
    end

    head_angle_balance0_std = std(mean(head_angle_balance0,2));
    head_angle_balance1_std = std(mean(head_angle_balance1,2));
    head_angle_balance2_std = std(mean(head_angle_balance2,2));
    head_angle_speed0_std = std(mean(head_angle_speed0,2));
    head_angle_speed1_std = std(mean(head_angle_speed1,2));
    head_angle_speed2_std = std(mean(head_angle_speed2,2));
    head_angle_accuracy0_std = std(mean(head_angle_accuracy0,2));
    head_angle_accuracy1_std = std(mean(head_angle_accuracy1,2));
    head_angle_accuracy2_std = std(mean(head_angle_accuracy2,2));

    head_angle_balance0_avg = mean(head_angle_balance0(:));
    head_angle_balance1_avg = mean(head_angle_balance1(:));
    head_angle_balance2_avg = mean(head_angle_balance2(:));
    head_angle_speed0_avg = mean(head_angle_speed0(:));
    head_angle_speed1_avg = mean(head_angle_speed1(:));
    head_angle_speed2_avg = mean(head_angle_speed2(:));
    head_angle_accuracy0_avg = mean(head_angle_accuracy0(:));
    head_angle_accuracy1_avg = mean(head_angle_accuracy1(:));
    head_angle_accuracy2_avg = mean(head_angle_accuracy2(:));

    means = [head_angle_balance0_avg, head_angle_balance1_avg, head_angle_balance2_avg; 
             head_angle_speed0_avg, head_angle_speed1_avg, head_angle_speed2_avg; 
             head_angle_accuracy0_avg, head_angle_accuracy1_avg, head_angle_accuracy2_avg];
    
    stds = [head_angle_balance0_std, head_angle_balance1_std, head_angle_balance2_std; 
            head_angle_speed0_std, head_angle_speed1_std, head_angle_speed2_std; 
            head_angle_accuracy0_std, head_angle_accuracy1_std, head_angle_accuracy2_std];
        
    figure;
    hold on;
    categories = {'Visual Perturbation', 'Speed Prompt', 'Accuracy Prompt'};
    x = 1:3; % Groups
    bar_handles = bar(x, means, 'grouped'); 

    if plot_eachparticipant == 1
        points_means = {mean(head_angle_balance0,2), mean(head_angle_balance1,2), mean(head_angle_balance2,2); 
                    mean(head_angle_speed0,2), mean(head_angle_speed1,2), mean(head_angle_speed2,2); 
                    mean(head_angle_accuracy0,2), mean(head_angle_accuracy1,2), mean(head_angle_accuracy2,2)};
        
        num_trials = size(points_means{1,1}, 1); % Assuming all groups have the same number of trials
        %colors = lines(num_trials); % Use 'lines' colormap to generate distinct colors
        %colors = min(colors * 1.4, 1); % Increase intensity, limit to 1
        for i = 1:size(means, 2) % Loop through balance, speed, accuracy columns
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars
            
            % Scatter plot for individual points with unique colors per row
            for j = 1:3 % Loop through categories (Balance, Speed, Accuracy)
                data_points = points_means{j, i}; % Extract row-wise mean values
                for k = 1:num_trials % Loop through trials (rows)
                    scatter(x_offset(j), data_points(k), point_size, colors_participants(k, :), 'filled', 'MarkerFaceAlpha', 0.6);
                end
            end
        end
    else
        for i = 1:size(means, 2) % Loop through each group
            % Error bars
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            errorbar(x_offset, means(:,i), stds(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1.5);
        end
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:3); 
    ylim([0 74])
    xticklabels(categories);
    %ylabel('Head Angle (deg)');
    ylabel('Mean Value (deg)')
    legend({'Null/Low', 'Medium', 'High'});
    title('Head Angle');
    grid on;
    hold off;

    %% Trunk angle by condition - PLOT 35

    trunk_angle_balance0 = [];
    trunk_angle_balance1 = [];
    trunk_angle_balance2 = [];
    trunk_angle_speed0 = [];
    trunk_angle_speed1 = [];
    trunk_angle_speed2 = [];
    trunk_angle_accuracy0 = [];
    trunk_angle_accuracy1 = [];
    trunk_angle_accuracy2 = [];

    for i=1:length(list_labels)
        current_label = list_labels{i};

        if contains(current_label, 's0')
            trunk_angle_speed0 = [trunk_angle_speed0, data_labels.(current_label).trunk_angle(:,1)];
        end
        if contains(current_label, 's1')
            trunk_angle_speed1 = [trunk_angle_speed1, data_labels.(current_label).trunk_angle(:,1)];
        end
        if contains(current_label, 's2')
            trunk_angle_speed2 = [trunk_angle_speed2, data_labels.(current_label).trunk_angle(:,1)];
        end

        if contains(current_label, 'a0')
            trunk_angle_accuracy0 = [trunk_angle_accuracy0, data_labels.(current_label).trunk_angle(:,1)];
        end
        if contains(current_label, 'a1')
            trunk_angle_accuracy1 = [trunk_angle_accuracy1, data_labels.(current_label).trunk_angle(:,1)];
        end
        if contains(current_label, 'a2')
            trunk_angle_accuracy2 = [trunk_angle_accuracy2, data_labels.(current_label).trunk_angle(:,1)];
        end

        if contains(current_label, 'b0')
            trunk_angle_balance0 = [trunk_angle_balance0, data_labels.(current_label).trunk_angle(:,1)];
        end
        if contains(current_label, 'b1')
            trunk_angle_balance1 = [trunk_angle_balance1, data_labels.(current_label).trunk_angle(:,1)];
        end
        if contains(current_label, 'b2')
            trunk_angle_balance2 = [trunk_angle_balance2, data_labels.(current_label).trunk_angle(:,1)];
        end
    end

    trunk_angle_balance0_std = std(mean(trunk_angle_balance0,2));
    trunk_angle_balance1_std = std(mean(trunk_angle_balance1,2));
    trunk_angle_balance2_std = std(mean(trunk_angle_balance2,2));
    trunk_angle_speed0_std = std(mean(trunk_angle_speed0,2));
    trunk_angle_speed1_std = std(mean(trunk_angle_speed1,2));
    trunk_angle_speed2_std = std(mean(trunk_angle_speed2,2));
    trunk_angle_accuracy0_std = std(mean(trunk_angle_accuracy0,2));
    trunk_angle_accuracy1_std = std(mean(trunk_angle_accuracy1,2));
    trunk_angle_accuracy2_std = std(mean(trunk_angle_accuracy2,2));

    trunk_angle_balance0_avg = mean(trunk_angle_balance0(:));
    trunk_angle_balance1_avg = mean(trunk_angle_balance1(:));
    trunk_angle_balance2_avg = mean(trunk_angle_balance2(:));
    trunk_angle_speed0_avg = mean(trunk_angle_speed0(:));
    trunk_angle_speed1_avg = mean(trunk_angle_speed1(:));
    trunk_angle_speed2_avg = mean(trunk_angle_speed2(:));
    trunk_angle_accuracy0_avg = mean(trunk_angle_accuracy0(:));
    trunk_angle_accuracy1_avg = mean(trunk_angle_accuracy1(:));
    trunk_angle_accuracy2_avg = mean(trunk_angle_accuracy2(:));

    means = [trunk_angle_balance0_avg, trunk_angle_balance1_avg, trunk_angle_balance2_avg; 
             trunk_angle_speed0_avg, trunk_angle_speed1_avg, trunk_angle_speed2_avg; 
             trunk_angle_accuracy0_avg, trunk_angle_accuracy1_avg, trunk_angle_accuracy2_avg];
    
    stds = [trunk_angle_balance0_std, trunk_angle_balance1_std, trunk_angle_balance2_std; 
            trunk_angle_speed0_std, trunk_angle_speed1_std, trunk_angle_speed2_std; 
            trunk_angle_accuracy0_std, trunk_angle_accuracy1_std, trunk_angle_accuracy2_std];
        
    figure;
    hold on;
    categories = {'Visual Perturbation', 'Speed Prompt', 'Accuracy Prompt'};
    x = 1:3; % Groups
    bar_handles = bar(x, means, 'grouped'); 

    if plot_eachparticipant == 1
        points_means = {mean(trunk_angle_balance0,2), mean(trunk_angle_balance1,2), mean(trunk_angle_balance2,2); 
                    mean(trunk_angle_speed0,2), mean(trunk_angle_speed1,2), mean(trunk_angle_speed2,2); 
                    mean(trunk_angle_accuracy0,2), mean(trunk_angle_accuracy1,2), mean(trunk_angle_accuracy2,2)};
        
        num_trials = size(points_means{1,1}, 1); % Assuming all groups have the same number of trials
        %colors = lines(num_trials); % Use 'lines' colormap to generate distinct colors
        %colors = min(colors * 1.4, 1); % Increase intensity, limit to 1
        for i = 1:size(means, 2) % Loop through balance, speed, accuracy columns
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars
            
            % Scatter plot for individual points with unique colors per row
            for j = 1:3 % Loop through categories (Balance, Speed, Accuracy)
                data_points = points_means{j, i}; % Extract row-wise mean values
                for k = 1:num_trials % Loop through trials (rows)
                    scatter(x_offset(j), data_points(k), point_size, colors_participants(k, :), 'filled', 'MarkerFaceAlpha', 0.6);
                end
            end
        end
    else
        for i = 1:size(means, 2) % Loop through each group
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            errorbar(x_offset, means(:,i), stds(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1.5);
        end
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:3); 
    ylim([0 26])
    xticklabels(categories);
    %ylabel('Trunk Angle (deg)');
    ylabel('Mean Value (deg)')
    legend({'Null/Low', 'Medium', 'High'});
    title('Trunk Angle');
    grid on;
    hold off;

    %% Head + Trunk angle by condition - PLOT 36

    means = [trunk_angle_balance0_avg + head_angle_balance0_avg, trunk_angle_balance1_avg + head_angle_balance1_avg, trunk_angle_balance2_avg + head_angle_balance2_avg; 
             trunk_angle_speed0_avg + head_angle_speed0_avg, trunk_angle_speed1_avg + head_angle_speed1_avg, trunk_angle_speed2_avg + head_angle_speed2_avg; 
             trunk_angle_accuracy0_avg + head_angle_accuracy0_avg, trunk_angle_accuracy1_avg + head_angle_accuracy1_avg, trunk_angle_accuracy2_avg + head_angle_accuracy2_avg];
    
    stds = [head_angle_balance0_std, head_angle_balance1_std, head_angle_balance2_std; 
            head_angle_speed0_std, head_angle_speed1_std, head_angle_speed2_std; 
            head_angle_accuracy0_std, head_angle_accuracy1_std, head_angle_accuracy2_std];
        
    figure;
    hold on;
    categories = {'Visual Perturbation', 'Speed Prompt', 'Accuracy Prompt'};
    x = 1:3; % Groups
    bar_handles = bar(x, means, 'grouped'); 
    for i = 1:size(means, 2) 
        x_offset = bar_handles(i).XEndPoints; % Get bar positions
        errorbar(x_offset, means(:,i), stds(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1.5);
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:3); 
    xticklabels(categories);
    ylabel('Head + Trunk Angle (deg)');
    legend({'Null/Low', 'Medium', 'High'});
    title('Head + Trunk Angle');
    grid on;
    hold off;

    %% Walking speed variability by condition - PLOT 37

    speed_std_balance0 = [];
    speed_std_balance1 = [];
    speed_std_balance2 = [];
    speed_std_speed0 = [];
    speed_std_speed1 = [];
    speed_std_speed2 = [];
    speed_std_accuracy0 = [];
    speed_std_accuracy1 = [];
    speed_std_accuracy2 = [];

    for i=1:length(list_labels)
        current_label = list_labels{i};

        if contains(current_label, 's0')
            speed_std_speed0 = [speed_std_speed0, data_labels.(current_label).walking_speed_std(:,1)];
        end
        if contains(current_label, 's1')
            speed_std_speed1 = [speed_std_speed1, data_labels.(current_label).walking_speed_std(:,1)];
        end
        if contains(current_label, 's2')
            speed_std_speed2 = [speed_std_speed2, data_labels.(current_label).walking_speed_std(:,1)];
        end

        if contains(current_label, 'a0')
            speed_std_accuracy0 = [speed_std_accuracy0, data_labels.(current_label).walking_speed_std(:,1)];
        end
        if contains(current_label, 'a1')
            speed_std_accuracy1 = [speed_std_accuracy1, data_labels.(current_label).walking_speed_std(:,1)];
        end
        if contains(current_label, 'a2')
            speed_std_accuracy2 = [speed_std_accuracy2, data_labels.(current_label).walking_speed_std(:,1)];
        end

        if contains(current_label, 'b0')
            speed_std_balance0 = [speed_std_balance0, data_labels.(current_label).walking_speed_std(:,1)];
        end
        if contains(current_label, 'b1')
            speed_std_balance1 = [speed_std_balance1, data_labels.(current_label).walking_speed_std(:,1)];
        end
        if contains(current_label, 'b2')
            speed_std_balance2 = [speed_std_balance2, data_labels.(current_label).walking_speed_std(:,1)];
        end
    end

    speed_std_balance0_std = nanstd(nanmean(speed_std_balance0, 2));
    speed_std_balance1_std = nanstd(nanmean(speed_std_balance1, 2));
    speed_std_balance2_std = nanstd(nanmean(speed_std_balance2, 2));
    speed_std_speed0_std = nanstd(nanmean(speed_std_speed0, 2));
    speed_std_speed1_std = nanstd(nanmean(speed_std_speed1, 2));
    speed_std_speed2_std = nanstd(nanmean(speed_std_speed2, 2));
    speed_std_accuracy0_std = nanstd(nanmean(speed_std_accuracy0, 2));
    speed_std_accuracy1_std = nanstd(nanmean(speed_std_accuracy1, 2));
    speed_std_accuracy2_std = nanstd(nanmean(speed_std_accuracy2, 2));
    
    speed_std_balance0_avg = nanmean(speed_std_balance0(:));
    speed_std_balance1_avg = nanmean(speed_std_balance1(:));
    speed_std_balance2_avg = nanmean(speed_std_balance2(:));
    speed_std_speed0_avg = nanmean(speed_std_speed0(:));
    speed_std_speed1_avg = nanmean(speed_std_speed1(:));
    speed_std_speed2_avg = nanmean(speed_std_speed2(:));
    speed_std_accuracy0_avg = nanmean(speed_std_accuracy0(:));
    speed_std_accuracy1_avg = nanmean(speed_std_accuracy1(:));
    speed_std_accuracy2_avg = nanmean(speed_std_accuracy2(:));


    means = [speed_std_balance0_avg, speed_std_balance1_avg, speed_std_balance2_avg; 
             speed_std_speed0_avg, speed_std_speed1_avg, speed_std_speed2_avg; 
             speed_std_accuracy0_avg, speed_std_accuracy1_avg, speed_std_accuracy2_avg];
    
    stds = [speed_std_balance0_std, speed_std_balance1_std, speed_std_balance2_std; 
            speed_std_speed0_std, speed_std_speed1_std, speed_std_speed2_std; 
            speed_std_accuracy0_std, speed_std_accuracy1_std, speed_std_accuracy2_std];
        
    figure;
    hold on;
    categories = {'Visual Perturbation', 'Speed Prompt', 'Accuracy Prompt'};
    x = 1:3; % Groups
    bar_handles = bar(x, means, 'grouped'); 

    if plot_eachparticipant == 1
        points_means = {mean(speed_std_balance0,2), mean(speed_std_balance1,2), mean(speed_std_balance2,2); 
                    mean(speed_std_speed0,2), mean(speed_std_speed1,2), mean(speed_std_speed2,2); 
                    mean(speed_std_accuracy0,2), mean(speed_std_accuracy1,2), mean(speed_std_accuracy2,2)};
        
        num_trials = size(points_means{1,1}, 1); % Assuming all groups have the same number of trials
        %colors = lines(num_trials); % Use 'lines' colormap to generate distinct colors
        %colors = min(colors * 1.4, 1); % Increase intensity, limit to 1
        for i = 1:size(means, 2) % Loop through balance, speed, accuracy columns
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars
            
            % Scatter plot for individual points with unique colors per row
            for j = 1:3 % Loop through categories (Balance, Speed, Accuracy)
                data_points = points_means{j, i}; % Extract row-wise mean values
                for k = 1:num_trials % Loop through trials (rows)
                    scatter(x_offset(j), data_points(k), point_size, colors_participants(k, :), 'filled', 'MarkerFaceAlpha', 0.6);
                end
            end
        end
    else
        for i = 1:size(means, 2) % Loop through each group
            % Error bars
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            errorbar(x_offset, means(:,i), stds(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1.5);
        end
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:3); 
    ylim([0 0.27])
    xticklabels(categories);
    %ylabel('Walking Speed Variability (m/s)');
    ylabel('Mean Value (m/s)')
    legend({'Null/Low', 'Medium', 'High'});
    title('Walking Speed Variability');
    grid on;
    hold off;

    %% Step error by all 27 conditions - PLOT 38-41

    for i=1:length(list_labels)
        current_label = list_labels{i};
        step_error_avg_array(i) = mean(data_labels.(current_label).step_error); 
        step_error_avg_std(i) = std(data_labels.(current_label).step_error);  
    end
    
    figure;
    hold on;

    % All trials - Sort by accuracy
    ordered_indexes = [1,2,3,10,11,12,19,20,21,4,5,6,13,14,15,22,23,24,7,8,9,16,17,18,25,26,27];
    bar_handle = bar(1:27, step_error_avg_array(ordered_indexes)); 
    
    hold on; % Keep the bars while adding error bars
    errorbar(1:27, step_error_avg_array(ordered_indexes), step_error_avg_std(ordered_indexes), 'k.', 'LineWidth', 1.5);
    hold off;
    
    for i = 1:27
        if i <= 9
            bar_handle.FaceColor = 'flat';
            bar_handle.CData(i, :) = colors(1, :);
        elseif i <= 18
            bar_handle.FaceColor = 'flat';
            bar_handle.CData(i, :) = colors(2, :);
        else
            bar_handle.FaceColor = 'flat';
            bar_handle.CData(i, :) = colors(3, :);
        end
    end

    set(gca, 'FontSize', 22);
    set(gca, 'XTick', 1:27);
    set(gca, 'XTickLabel', list_labels(ordered_indexes));
    xtickangle(45);
    xlabel('Trials');
    ylabel('Mean Step Error (cm)');
    title('Mean Step Error by Trial');
    
    hold off;

    % Accuracy 1,2
    selected_indices = [4, 5, 6, 7, 8, 9, 13, 14, 15, 16, 17, 18, 22, 23, 24, 25, 26, 27];

    selected_means = step_error_avg_array(selected_indices);
    selected_stds = step_error_avg_std(selected_indices);
    selected_labels = list_labels(selected_indices);
    
    figure;
    hold on;
    
    bar_handle = bar(1:length(selected_indices), selected_means);
    errorbar(1:length(selected_indices), selected_means, selected_stds, 'k.', 'LineWidth', 1.5);
    
    set(gca, 'FontSize', 22);
    set(gca, 'XTick', 1:length(selected_indices));
    set(gca, 'XTickLabel', selected_labels);
    xtickangle(45);
    xlabel('Medium and High Accuracy Trials');
    ylabel('Mean Step Error (cm)');
    title('Mean Step Error - Medium and High Accuracy Trials');
    hold off;  

    % Accuracy 1
    selected_indices = [4, 5, 6, 13, 14, 15, 22, 23, 24];

    selected_means = step_error_avg_array(selected_indices);
    selected_stds = step_error_avg_std(selected_indices);
    selected_labels = list_labels(selected_indices);
    
    figure;
    hold on;
    
    bar_handle = bar(1:length(selected_indices), selected_means);
    errorbar(1:length(selected_indices), selected_means, selected_stds, 'k.', 'LineWidth', 1.5);
    
    set(gca, 'FontSize', 22);
    set(gca, 'XTick', 1:length(selected_indices));
    set(gca, 'XTickLabel', selected_labels);
    xtickangle(45);
    xlabel('Medium Accuracy Trials');
    ylabel('Mean Step Error (cm)');
    title('Mean Step Error - Medium Accuracy Trials');
    hold off;  

    % Accuracy 2
    selected_indices = [7, 8, 9, 16, 17, 18, 25, 26, 27];

    selected_means = step_error_avg_array(selected_indices);
    selected_stds = step_error_avg_std(selected_indices);
    selected_labels = list_labels(selected_indices);
    
    figure;
    hold on;
    
    bar_handle = bar(1:length(selected_indices), selected_means);
    errorbar(1:length(selected_indices), selected_means, selected_stds, 'k.', 'LineWidth', 1.5);
    
    set(gca, 'FontSize', 22);
    set(gca, 'XTick', 1:length(selected_indices));
    set(gca, 'XTickLabel', selected_labels);
    xtickangle(45);
    xlabel('High Accuracy Trials');
    ylabel('Mean Step Error (cm)');
    title('Mean Step Error - High Accuracy Trials');
    hold off; 

    %% Stride Duration - PLOT 42

    stride_time_balance0 = [];
    stride_time_balance1 = [];
    stride_time_balance2 = [];
    stride_time_speed0 = [];
    stride_time_speed1 = [];
    stride_time_speed2 = [];
    stride_time_accuracy0 = [];
    stride_time_accuracy1 = [];
    stride_time_accuracy2 = [];

    for i=1:length(list_labels)
        current_label = list_labels{i};

        if contains(current_label, 's0')
            stride_time_speed0 = [stride_time_speed0, data_labels.(current_label).stride_time(:,1)];
        end
        if contains(current_label, 's1')
            stride_time_speed1 = [stride_time_speed1, data_labels.(current_label).stride_time(:,1)];
        end
        if contains(current_label, 's2')
            stride_time_speed2 = [speed_speed2, data_labels.(current_label).stride_time(:,1)];
        end

        if contains(current_label, 'a0')
            stride_time_accuracy0 = [stride_time_accuracy0, data_labels.(current_label).stride_time(:,1)];
        end
        if contains(current_label, 'a1')
            stride_time_accuracy1 = [stride_time_accuracy1, data_labels.(current_label).stride_time(:,1)];
        end
        if contains(current_label, 'a2')
            stride_time_accuracy2 = [stride_time_accuracy2, data_labels.(current_label).stride_time(:,1)];
        end

        if contains(current_label, 'b0')
            stride_time_balance0 = [stride_time_balance0, data_labels.(current_label).stride_time(:,1)];
        end
        if contains(current_label, 'b1')
            stride_time_balance1 = [stride_time_balance1, data_labels.(current_label).stride_time(:,1)];
        end
        if contains(current_label, 'b2')
            stride_time_balance2 = [stride_time_balance2, data_labels.(current_label).stride_time(:,1)];
        end
    end

    stride_time_balance0_std = std(mean(stride_time_balance0,2));
    stride_time_balance1_std = std(mean(stride_time_balance1,2));
    stride_time_balance2_std = std(mean(stride_time_balance2,2));
    stride_time_speed0_std = std(mean(stride_time_speed0,2));
    stride_time_speed1_std = std(mean(stride_time_speed1,2));
    stride_time_speed2_std = std(mean(stride_time_speed2,2));
    stride_time_accuracy0_std = std(mean(stride_time_accuracy0,2));
    stride_time_accuracy1_std = std(mean(stride_time_accuracy1,2));
    stride_time_accuracy2_std = std(mean(stride_time_accuracy2,2));

    stride_time_balance0_avg = mean(mean(stride_time_balance0,2));
    stride_time_balance1_avg = mean(mean(stride_time_balance1,2));
    stride_time_balance2_avg = mean(mean(stride_time_balance2,2));
    stride_time_speed0_avg = mean(mean(stride_time_speed0,2));
    stride_time_speed1_avg = mean(mean(stride_time_speed1,2));
    stride_time_speed2_avg = mean(mean(stride_time_speed2,2));
    stride_time_accuracy0_avg = mean(mean(stride_time_accuracy0,2));
    stride_time_accuracy1_avg = mean(mean(stride_time_accuracy1,2));
    stride_time_accuracy2_avg = mean(mean(stride_time_accuracy2,2));

    means = [stride_time_balance0_avg, stride_time_balance1_avg, stride_time_balance2_avg; 
             stride_time_speed0_avg, stride_time_speed1_avg, stride_time_speed2_avg; 
             stride_time_accuracy0_avg, stride_time_accuracy1_avg, stride_time_accuracy2_avg];
    
    stds = [stride_time_balance0_std, stride_time_balance1_std, stride_time_balance2_std; 
            stride_time_speed0_std, stride_time_speed1_std, stride_time_speed2_std; 
            stride_time_accuracy0_std, stride_time_accuracy1_std, stride_time_accuracy2_std];
        
    figure;
    hold on;
    categories = {'Visual Perturbation', 'Speed Prompt', 'Accuracy Prompt'};
    x = 1:3; % Groups
    bar_handles = bar(x, means, 'grouped'); 

    if plot_eachparticipant == 1
        points_means = {mean(stride_time_balance0,2), mean(stride_time_balance1,2), mean(stride_time_balance2,2); 
                    mean(stride_time_speed0,2), mean(stride_time_speed1,2), mean(stride_time_speed2,2); 
                    mean(stride_time_accuracy0,2), mean(stride_time_accuracy1,2), mean(stride_time_accuracy2,2)};
        
        num_trials = size(points_means{1,1}, 1); % Assuming all groups have the same number of trials
        %colors = lines(num_trials); % Use 'lines' colormap to generate distinct colors
        %colors = min(colors * 1.4, 1); % Increase intensity, limit to 1
        for i = 1:size(means, 2) % Loop through balance, speed, accuracy columns
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars
            
            % Scatter plot for individual points with unique colors per row
            for j = 1:3 % Loop through categories (Balance, Speed, Accuracy)
                data_points = points_means{j, i}; % Extract row-wise mean values
                for k = 1:num_trials % Loop through trials (rows)
                    scatter(x_offset(j), data_points(k), point_size, colors_participants(k, :), 'filled', 'MarkerFaceAlpha', 0.6);
                end
            end
        end
    else
        for i = 1:size(means, 2) % Loop through each group
            % Error bars
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            errorbar(x_offset, means(:,i), stds(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1.5);
        end
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:3); 
    ylim([0 2.1])
    xticklabels(categories);
    %ylabel('Stride Duration (s)');
    ylabel('Mean Value (s)')
    legend({'Null/Low', 'Medium', 'High'});
    title('Stride Duration');
    grid on;
    hold off;
    
    %% Stride Duration variability - PLOT 43

    stride_time_std_balance0 = [];
    stride_time_std_balance1 = [];
    stride_time_std_balance2 = [];
    stride_time_std_speed0 = [];
    stride_time_std_speed1 = [];
    stride_time_std_speed2 = [];
    stride_time_std_accuracy0 = [];
    stride_time_std_accuracy1 = [];
    stride_time_std_accuracy2 = [];

    for i=1:length(list_labels)
        current_label = list_labels{i};

        if contains(current_label, 's0')
            stride_time_std_speed0 = [stride_time_std_speed0, data_labels.(current_label).stride_time_var(:,1)];
        end
        if contains(current_label, 's1')
            stride_time_std_speed1 = [stride_time_std_speed1, data_labels.(current_label).stride_time_var(:,1)];
        end
        if contains(current_label, 's2')
            stride_time_std_speed2 = [stride_time_std_speed2, data_labels.(current_label).stride_time_var(:,1)];
        end

        if contains(current_label, 'a0')
            stride_time_std_accuracy0 = [stride_time_std_accuracy0, data_labels.(current_label).stride_time_var(:,1)];
        end
        if contains(current_label, 'a1')
            stride_time_std_accuracy1 = [stride_time_std_accuracy1, data_labels.(current_label).stride_time_var(:,1)];
        end
        if contains(current_label, 'a2')
            stride_time_std_accuracy2 = [stride_time_std_accuracy2, data_labels.(current_label).stride_time_var(:,1)];
        end

        if contains(current_label, 'b0')
            stride_time_std_balance0 = [stride_time_std_balance0, data_labels.(current_label).stride_time_var(:,1)];
        end
        if contains(current_label, 'b1')
            stride_time_std_balance1 = [stride_time_std_balance1, data_labels.(current_label).stride_time_var(:,1)];
        end
        if contains(current_label, 'b2')
            stride_time_std_balance2 = [stride_time_std_balance2, data_labels.(current_label).stride_time_var(:,1)];
        end
    end

    % Remove this value manually because it is messing up the stds (it is correct though, I checked the qtm BMH02 Trial5)
    idx = abs(stride_time_std_accuracy2 - 1.4525) <= 0.01;
    stride_time_std_accuracy2(idx) = NaN;
    idx = abs(stride_time_std_balance2 - 1.4525) <= 0.01;
    stride_time_std_balance2(idx) = NaN;
    idx = abs(stride_time_std_speed0 - 1.4525) <= 0.01;
    stride_time_std_speed0(idx) = NaN;

    stride_time_std_balance0_std = std(mean(stride_time_std_balance0,2), 'omitnan');
    stride_time_std_balance1_std = std(mean(stride_time_std_balance1,2), 'omitnan');
    stride_time_std_balance2_std = std(mean(stride_time_std_balance2,2), 'omitnan');
    stride_time_std_speed0_std = std(mean(stride_time_std_speed0,2), 'omitnan');
    stride_time_std_speed1_std = std(mean(stride_time_std_speed1,2), 'omitnan');
    stride_time_std_speed2_std = std(mean(stride_time_std_speed2,2), 'omitnan');
    stride_time_std_accuracy0_std = std(mean(stride_time_std_accuracy0,2), 'omitnan');
    stride_time_std_accuracy1_std = std(mean(stride_time_std_accuracy1,2), 'omitnan');
    stride_time_std_accuracy2_std = std(mean(stride_time_std_accuracy2,2), 'omitnan');

    stride_time_std_balance0_avg = mean(mean(stride_time_std_balance0,2), 'omitnan');
    stride_time_std_balance1_avg = mean(mean(stride_time_std_balance1,2), 'omitnan');
    stride_time_std_balance2_avg = mean(mean(stride_time_std_balance2,2), 'omitnan');
    stride_time_std_speed0_avg = mean(mean(stride_time_std_speed0,2), 'omitnan');
    stride_time_std_speed1_avg = mean(mean(stride_time_std_speed1,2), 'omitnan');
    stride_time_std_speed2_avg = mean(mean(stride_time_std_speed2,2), 'omitnan');
    stride_time_std_accuracy0_avg = mean(mean(stride_time_std_accuracy0,2), 'omitnan');
    stride_time_std_accuracy1_avg = mean(mean(stride_time_std_accuracy1,2), 'omitnan');
    stride_time_std_accuracy2_avg = mean(mean(stride_time_std_accuracy2,2), 'omitnan');

    means = [stride_time_std_balance0_avg, stride_time_std_balance1_avg, stride_time_std_balance2_avg; 
             stride_time_std_speed0_avg, stride_time_std_speed1_avg, stride_time_std_speed2_avg; 
             stride_time_std_accuracy0_avg, stride_time_std_accuracy1_avg, stride_time_std_accuracy2_avg];
    
    stds = [stride_time_std_balance0_std, stride_time_std_balance1_std, stride_time_std_balance2_std; 
            stride_time_std_speed0_std, stride_time_std_speed1_std, stride_time_std_speed2_std; 
            stride_time_std_accuracy0_std, stride_time_std_accuracy1_std, stride_time_std_accuracy2_std];
        
    figure;
    hold on;
    categories = {'Visual Perturbation', 'Speed Prompt', 'Accuracy Prompt'};
    x = 1:3; % Groups
    bar_handles = bar(x, means, 'grouped'); 
    
    if plot_eachparticipant == 1
        points_means = {mean(stride_time_std_balance0,2), mean(stride_time_std_balance1,2), mean(stride_time_std_balance2,2); 
                    mean(stride_time_std_speed0,2), mean(stride_time_std_speed1,2), mean(stride_time_std_speed2,2); 
                    mean(stride_time_std_accuracy0,2), mean(stride_time_std_accuracy1,2), mean(stride_time_std_accuracy2,2)};
        
        num_trials = size(points_means{1,1}, 1); % Assuming all groups have the same number of trials
        %colors = lines(num_trials); % Use 'lines' colormap to generate distinct colors
        %colors = min(colors * 1.4, 1); % Increase intensity, limit to 1
        for i = 1:size(means, 2) % Loop through balance, speed, accuracy columns
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars
            
            % Scatter plot for individual points with unique colors per row
            for j = 1:3 % Loop through categories (Balance, Speed, Accuracy)
                data_points = points_means{j, i}; % Extract row-wise mean values
                for k = 1:num_trials % Loop through trials (rows)
                    scatter(x_offset(j), data_points(k), point_size, colors_participants(k, :), 'filled', 'MarkerFaceAlpha', 0.6);
                end
            end
        end
    else
        for i = 1:size(means, 2) % Loop through each group
            % Error bars
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            errorbar(x_offset, means(:,i), stds(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1.5);
        end
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:3); 
    ylim([0 0.37])
    xticklabels(categories);
    %ylabel('Stride Duration Variability (s)');
    ylabel('Mean Value (s)')
    legend({'Null/Low', 'Medium', 'High'});
    title('Stride Duration Variability');
    grid on;
    hold off;

    %% Step error variability - PLOT 44

    step_error_std_balance0 = [];
    step_error_std_balance1 = [];
    step_error_std_balance2 = [];
    step_error_std_speed0 = [];
    step_error_std_speed1 = [];
    step_error_std_speed2 = [];
    step_error_std_accuracy0 = [];
    step_error_std_accuracy1 = [];
    step_error_std_accuracy2 = [];

    for i=1:length(list_labels)
        current_label = list_labels{i};

        if contains(current_label, 's0')
            step_error_std_speed0 = [step_error_std_speed0, data_labels.(current_label).step_error_std(:,1)];
        end
        if contains(current_label, 's1')
            step_error_std_speed1 = [step_error_std_speed1, data_labels.(current_label).step_error_std(:,1)];
        end
        if contains(current_label, 's2')
            step_error_std_speed2 = [step_error_std_speed2, data_labels.(current_label).step_error_std(:,1)];
        end

        if contains(current_label, 'a0')
            step_error_std_accuracy0 = [step_error_std_accuracy0, data_labels.(current_label).step_error_std(:,1)];
        end
        if contains(current_label, 'a1')
            step_error_std_accuracy1 = [step_error_std_accuracy1, data_labels.(current_label).step_error_std(:,1)];
        end
        if contains(current_label, 'a2')
            step_error_std_accuracy2 = [step_error_std_accuracy2, data_labels.(current_label).step_error_std(:,1)];
        end

        if contains(current_label, 'b0')
            step_error_std_balance0 = [step_error_std_balance0, data_labels.(current_label).step_error_std(:,1)];
        end
        if contains(current_label, 'b1')
            step_error_std_balance1 = [step_error_std_balance1, data_labels.(current_label).step_error_std(:,1)];
        end
        if contains(current_label, 'b2')
            step_error_std_balance2 = [step_error_std_balance2, data_labels.(current_label).step_error_std(:,1)];
        end
    end

    step_error_std_balance0_std = std(mean(step_error_std_balance0,2));
    step_error_std_balance1_std = std(mean(step_error_std_balance1,2));
    step_error_std_balance2_std = std(mean(step_error_std_balance2,2));
    step_error_std_speed0_std = std(mean(step_error_std_speed0,2));
    step_error_std_speed1_std = std(mean(step_error_std_speed1,2));
    step_error_std_speed2_std = std(mean(step_error_std_speed2,2));
    step_error_std_accuracy0_std = std(mean(step_error_std_accuracy0,2));
    step_error_std_accuracy1_std = std(mean(step_error_std_accuracy1,2));
    step_error_std_accuracy2_std = std(mean(step_error_std_accuracy2,2));

    step_error_std_balance0_avg = mean(step_error_std_balance0(:));
    step_error_std_balance1_avg = mean(step_error_std_balance1(:));
    step_error_std_balance2_avg = mean(step_error_std_balance2(:));
    step_error_std_speed0_avg = mean(step_error_std_speed0(:));
    step_error_std_speed1_avg = mean(step_error_std_speed1(:));
    step_error_std_speed2_avg = mean(step_error_std_speed2(:));
    step_error_std_accuracy0_avg = mean(step_error_std_accuracy0(:));
    step_error_std_accuracy1_avg = mean(step_error_std_accuracy1(:));
    step_error_std_accuracy2_avg = mean(step_error_std_accuracy2(:));

    means = [step_error_std_balance0_avg, step_error_std_balance1_avg, step_error_std_balance2_avg; 
             step_error_std_speed0_avg, step_error_std_speed1_avg, step_error_std_speed2_avg; 
             step_error_std_accuracy0_avg, step_error_std_accuracy1_avg, step_error_std_accuracy2_avg];
    
    stds = [step_error_std_balance0_std, step_error_std_balance1_std, step_error_std_balance2_std; 
            step_error_std_speed0_std, step_error_std_speed1_std, step_error_std_speed2_std; 
            step_error_std_accuracy0_std, step_error_std_accuracy1_std, step_error_std_accuracy2_std];
        
    figure;
    hold on;
    categories = {'Visual Perturbation', 'Speed Prompt', 'Accuracy Prompt'};
    x = 1:3; % Groups
    bar_handles = bar(x, means, 'grouped'); 

    if plot_eachparticipant == 1
        points_means = {mean(step_error_std_balance0,2), mean(step_error_std_balance1,2), mean(step_error_std_balance2,2); 
                    mean(step_error_std_speed0,2), mean(step_error_std_speed1,2), mean(step_error_std_speed2,2); 
                    mean(step_error_std_accuracy0,2), mean(step_error_std_accuracy1,2), mean(step_error_std_accuracy2,2)};
        
        num_trials = size(points_means{1,1}, 1); % Assuming all groups have the same number of trials
        %colors = lines(num_trials); % Use 'lines' colormap to generate distinct colors
        %colors = min(colors * 1.4, 1); % Increase intensity, limit to 1
        for i = 1:size(means, 2) % Loop through balance, speed, accuracy columns
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars
            
            % Scatter plot for individual points with unique colors per row
            for j = 1:3 % Loop through categories (Balance, Speed, Accuracy)
                data_points = points_means{j, i}; % Extract row-wise mean values
                for k = 1:num_trials % Loop through trials (rows)
                    scatter(x_offset(j), data_points(k), point_size, colors_participants(k, :), 'filled', 'MarkerFaceAlpha', 0.6);
                end
            end
        end
    else
        for i = 1:size(means, 2) % Loop through each group
            % Error bars
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            errorbar(x_offset, means(:,i), stds(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1.5);
        end
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:3); 
    xticklabels(categories);
    %ylabel('Step Error Variability (cm)');
    ylabel('Mean Value (cm)')
    legend({'Null/Low', 'Medium', 'High'});
    title('Step Error Variability');
    grid on;
    hold off;

    %% STJ-headset angle - PLOT 45

    stj_headset_angle_balance0 = [];
    stj_headset_angle_balance1 = [];
    stj_headset_angle_balance2 = [];
    stj_headset_angle_speed0 = [];
    stj_headset_angle_speed1 = [];
    stj_headset_angle_speed2 = [];
    stj_headset_angle_accuracy0 = [];
    stj_headset_angle_accuracy1 = [];
    stj_headset_angle_accuracy2 = [];

    for i=1:length(list_labels)
        current_label = list_labels{i};

        if contains(current_label, 's0')
            stj_headset_angle_speed0 = [stj_headset_angle_speed0, data_labels.(current_label).stj_headset_angle(:,1)];
        end
        if contains(current_label, 's1')
            stj_headset_angle_speed1 = [stj_headset_angle_speed1, data_labels.(current_label).stj_headset_angle(:,1)];
        end
        if contains(current_label, 's2')
            stj_headset_angle_speed2 = [stj_headset_angle_speed2, data_labels.(current_label).stj_headset_angle(:,1)];
        end

        if contains(current_label, 'a0')
            stj_headset_angle_accuracy0 = [stj_headset_angle_accuracy0, data_labels.(current_label).stj_headset_angle(:,1)];
        end
        if contains(current_label, 'a1')
            stj_headset_angle_accuracy1 = [stj_headset_angle_accuracy1, data_labels.(current_label).stj_headset_angle(:,1)];
        end
        if contains(current_label, 'a2')
            stj_headset_angle_accuracy2 = [stj_headset_angle_accuracy2, data_labels.(current_label).stj_headset_angle(:,1)];
        end

        if contains(current_label, 'b0')
            stj_headset_angle_balance0 = [stj_headset_angle_balance0, data_labels.(current_label).stj_headset_angle(:,1)];
        end
        if contains(current_label, 'b1')
            stj_headset_angle_balance1 = [stj_headset_angle_balance1, data_labels.(current_label).stj_headset_angle(:,1)];
        end
        if contains(current_label, 'b2')
            stj_headset_angle_balance2 = [stj_headset_angle_balance2, data_labels.(current_label).stj_headset_angle(:,1)];
        end
    end

    stj_headset_angle_balance0_std = std(mean(stj_headset_angle_balance0,2), 'omitnan');
    stj_headset_angle_balance1_std = std(mean(stj_headset_angle_balance1,2), 'omitnan');
    stj_headset_angle_balance2_std = std(mean(stj_headset_angle_balance2,2), 'omitnan');
    stj_headset_angle_speed0_std = std(mean(stj_headset_angle_speed0,2), 'omitnan');
    stj_headset_angle_speed1_std = std(mean(stj_headset_angle_speed1,2), 'omitnan');
    stj_headset_angle_speed2_std = std(mean(stj_headset_angle_speed2,2), 'omitnan');
    stj_headset_angle_accuracy0_std = std(mean(stj_headset_angle_accuracy0,2), 'omitnan');
    stj_headset_angle_accuracy1_std = std(mean(stj_headset_angle_accuracy1,2), 'omitnan');
    stj_headset_angle_accuracy2_std = std(mean(stj_headset_angle_accuracy2,2), 'omitnan');

    stj_headset_angle_balance0 = mean(stj_headset_angle_balance0(:), 'omitnan');
    stj_headset_angle_balance1 = mean(stj_headset_angle_balance1(:), 'omitnan');
    stj_headset_angle_balance2 = mean(stj_headset_angle_balance2(:), 'omitnan');
    stj_headset_angle_speed0 = mean(stj_headset_angle_speed0(:), 'omitnan');
    stj_headset_angle_speed1 = mean(stj_headset_angle_speed1(:), 'omitnan');
    stj_headset_angle_speed2 = mean(stj_headset_angle_speed2(:), 'omitnan');
    stj_headset_angle_accuracy0 = mean(stj_headset_angle_accuracy0(:), 'omitnan');
    stj_headset_angle_accuracy1 = mean(stj_headset_angle_accuracy1(:), 'omitnan');
    stj_headset_angle_accuracy2 = mean(stj_headset_angle_accuracy2(:), 'omitnan');


    means = [stj_headset_angle_balance0, stj_headset_angle_balance1, stj_headset_angle_balance2; 
             stj_headset_angle_speed0, stj_headset_angle_speed1, stj_headset_angle_speed2; 
             stj_headset_angle_accuracy0, stj_headset_angle_accuracy1, stj_headset_angle_accuracy2];
    
    stds = [stj_headset_angle_balance0_std, stj_headset_angle_balance1_std, stj_headset_angle_balance2_std; 
            stj_headset_angle_speed0_std, stj_headset_angle_speed1_std, stj_headset_angle_speed2_std; 
            stj_headset_angle_accuracy0_std, stj_headset_angle_accuracy1_std, stj_headset_angle_accuracy2_std];
        
    figure;
    hold on;
    categories = {'Visual Perturbation', 'Speed Prompt', 'Accuracy Prompt'};
    x = 1:3; % Groups
    bar_handles = bar(x, means, 'grouped'); 
    for i = 1:size(means, 2) 
        x_offset = bar_handles(i).XEndPoints; % Get bar positions
        errorbar(x_offset, means(:,i), stds(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1.5);
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:3); 
    xticklabels(categories);
    ylabel('STJ-Headset Angle (deg)');
    legend({'Null/Low', 'Medium', 'High'});
    title('Mean STJ-Headset Angle');
    grid on;
    hold off;

    %% Where are they looking at (dist_gaze)? - PLOT 46

    dist_gaze_balance0 = [];
    dist_gaze_balance1 = [];
    dist_gaze_balance2 = [];
    dist_gaze_speed0 = [];
    dist_gaze_speed1 = [];
    dist_gaze_speed2 = [];
    dist_gaze_accuracy0 = [];
    dist_gaze_accuracy1 = [];
    dist_gaze_accuracy2 = [];

    for i=1:length(list_labels)
        current_label = list_labels{i};
        
        distances_label(:,1) = (min(data_labels.(current_label).dist_gaze(:,1), 400))/100; % if any value is bigger than 400 cm (4m), change it to that

        if contains(current_label, 's0')
            dist_gaze_speed0 = [dist_gaze_speed0, distances_label(:,1)];
        end
        if contains(current_label, 's1')
            dist_gaze_speed1 = [dist_gaze_speed1, distances_label(:,1)];
        end
        if contains(current_label, 's2')
            dist_gaze_speed2 = [dist_gaze_speed2, distances_label(:,1)];
        end

        if contains(current_label, 'a0')
            dist_gaze_accuracy0 = [dist_gaze_accuracy0, distances_label(:,1)];
        end
        if contains(current_label, 'a1')
            dist_gaze_accuracy1 = [dist_gaze_accuracy1, distances_label(:,1)];
        end
        if contains(current_label, 'a2')
            dist_gaze_accuracy2 = [dist_gaze_accuracy2, distances_label(:,1)];
        end

        if contains(current_label, 'b0')
            dist_gaze_balance0 = [dist_gaze_balance0, distances_label(:,1)];
        end
        if contains(current_label, 'b1')
            dist_gaze_balance1 = [dist_gaze_balance1, distances_label(:,1)];
        end
        if contains(current_label, 'b2')
            dist_gaze_balance2 = [dist_gaze_balance2, distances_label(:,1)];
        end
    end

    dist_gaze_balance0_std = std(mean(dist_gaze_balance0,2));
    dist_gaze_balance1_std = std(mean(dist_gaze_balance1,2));
    dist_gaze_balance2_std = std(mean(dist_gaze_balance2,2));
    dist_gaze_speed0_std = std(mean(dist_gaze_speed0,2));
    dist_gaze_speed1_std = std(mean(dist_gaze_speed1,2));
    dist_gaze_speed2_std = std(mean(dist_gaze_speed2,2));
    dist_gaze_accuracy0_std = std(mean(dist_gaze_accuracy0,2));
    dist_gaze_accuracy1_std = std(mean(dist_gaze_accuracy1,2));
    dist_gaze_accuracy2_std = std(mean(dist_gaze_accuracy2,2));

    dist_gaze_balance0_avg = mean(dist_gaze_balance0(:));
    dist_gaze_balance1_avg = mean(dist_gaze_balance1(:));
    dist_gaze_balance2_avg = mean(dist_gaze_balance2(:));
    dist_gaze_speed0_avg = mean(dist_gaze_speed0(:));
    dist_gaze_speed1_avg = mean(dist_gaze_speed1(:));
    dist_gaze_speed2_avg = mean(dist_gaze_speed2(:));
    dist_gaze_accuracy0_avg = mean(dist_gaze_accuracy0(:));
    dist_gaze_accuracy1_avg = mean(dist_gaze_accuracy1(:));
    dist_gaze_accuracy2_avg = mean(dist_gaze_accuracy2(:));

    means = [dist_gaze_balance0_avg, dist_gaze_balance1_avg, dist_gaze_balance2_avg; 
             dist_gaze_speed0_avg, dist_gaze_speed1_avg, dist_gaze_speed2_avg; 
             dist_gaze_accuracy0_avg, dist_gaze_accuracy1_avg, dist_gaze_accuracy2_avg];
    
    stds = [dist_gaze_balance0_std, dist_gaze_balance1_std, dist_gaze_balance2_std; 
            dist_gaze_speed0_std, dist_gaze_speed1_std, dist_gaze_speed2_std; 
            dist_gaze_accuracy0_std, dist_gaze_accuracy1_std, dist_gaze_accuracy2_std];
        
    figure;
    hold on;
    categories = {'Visual Perturbation', 'Speed Prompt', 'Accuracy Prompt'};
    x = 1:3; % Groups
    bar_handles = bar(x, means, 'grouped'); 

    if plot_eachparticipant == 1
        points_means = {mean(dist_gaze_balance0,2), mean(dist_gaze_balance1,2), mean(dist_gaze_balance2,2); 
                    mean(dist_gaze_speed0,2), mean(dist_gaze_speed1,2), mean(dist_gaze_speed2,2); 
                    mean(dist_gaze_accuracy0,2), mean(dist_gaze_accuracy1,2), mean(dist_gaze_accuracy2,2)};
        
        num_trials = size(points_means{1,1}, 1); % Assuming all groups have the same number of trials
        %colors = lines(num_trials); % Use 'lines' colormap to generate distinct colors
        %colors = min(colors * 1.4, 1); % Increase intensity, limit to 1
        for i = 1:size(means, 2) % Loop through balance, speed, accuracy columns
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars
            
            % Scatter plot for individual points with unique colors per row
            for j = 1:3 % Loop through categories (Balance, Speed, Accuracy)
                data_points = points_means{j, i}; % Extract row-wise mean values
                for k = 1:num_trials % Loop through trials (rows)
                    scatter(x_offset(j), data_points(k), point_size, colors_participants(k, :), 'filled', 'MarkerFaceAlpha', 0.6);
                end
            end
        end
    else
        for i = 1:size(means, 2) % Loop through each group
            % Error bars
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            errorbar(x_offset, means(:,i), stds(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1.5);
        end
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:3);
    ylim([0 5])
    xticklabels(categories);
    %ylabel('Look-Ahead Distance (m)');
    ylabel('Mean Value (cm)')
    legend({'Null/Low', 'Medium', 'High'});
    title('Look-Ahead Distance');
    grid on;
    hold off;

    %% Step error by condition - Plot 47
    
    step_error_balance0 = [];
    step_error_balance1 = [];
    step_error_balance2 = [];
    step_error_speed0 = [];
    step_error_speed1 = [];
    step_error_speed2 = [];
    step_error_accuracy0 = [];
    step_error_accuracy1 = [];
    step_error_accuracy2 = [];

    for i=1:length(list_labels)
        current_label = list_labels{i};

        if contains(current_label, 's0')
            step_error_speed0 = [step_error_speed0, data_labels.(current_label).step_error(:,1)];
        end
        if contains(current_label, 's1')
            step_error_speed1 = [step_error_speed1, data_labels.(current_label).step_error(:,1)];
        end
        if contains(current_label, 's2')
            step_error_speed2 = [step_error_speed2, data_labels.(current_label).step_error(:,1)];
        end

        if contains(current_label, 'a0')
            step_error_accuracy0 = [step_error_accuracy0, data_labels.(current_label).step_error(:,1)];
        end
        if contains(current_label, 'a1')
            step_error_accuracy1 = [step_error_accuracy1, data_labels.(current_label).step_error(:,1)];
        end
        if contains(current_label, 'a2')
            step_error_accuracy2 = [step_error_accuracy2, data_labels.(current_label).step_error(:,1)];
        end

        if contains(current_label, 'b0')
            step_error_balance0 = [step_error_balance0, data_labels.(current_label).step_error(:,1)];
        end
        if contains(current_label, 'b1')
            step_error_balance1 = [step_error_balance1, data_labels.(current_label).step_error(:,1)];
        end
        if contains(current_label, 'b2')
            step_error_balance2 = [step_error_balance2, data_labels.(current_label).step_error(:,1)];
        end
    end

    step_error_balance0_std = std(mean(step_error_balance0,2));
    step_error_balance1_std = std(mean(step_error_balance1,2));
    step_error_balance2_std = std(mean(step_error_balance2,2));
    step_error_speed0_std = std(mean(step_error_speed0,2));
    step_error_speed1_std = std(mean(step_error_speed1,2));
    step_error_speed2_std = std(mean(step_error_speed2,2));
    step_error_accuracy0_std = std(mean(step_error_accuracy0,2));
    step_error_accuracy1_std = std(mean(step_error_accuracy1,2));
    step_error_accuracy2_std = std(mean(step_error_accuracy2,2));

    step_error_balance0_avg = mean(step_error_balance0(:));
    step_error_balance1_avg = mean(step_error_balance1(:));
    step_error_balance2_avg = mean(step_error_balance2(:));
    step_error_speed0_avg = mean(step_error_speed0(:));
    step_error_speed1_avg = mean(step_error_speed1(:));
    step_error_speed2_avg = mean(step_error_speed2(:));
    step_error_accuracy0_avg = mean(step_error_accuracy0(:));
    step_error_accuracy1_avg = mean(step_error_accuracy1(:));
    step_error_accuracy2_avg = mean(step_error_accuracy2(:));

    means = [step_error_balance0_avg, step_error_balance1_avg, step_error_balance2_avg; 
             step_error_speed0_avg, step_error_speed1_avg, step_error_speed2_avg; 
             step_error_accuracy0_avg, step_error_accuracy1_avg, step_error_accuracy2_avg];

    stds = [step_error_balance0_std, step_error_balance1_std, step_error_balance2_std; 
            step_error_speed0_std, step_error_speed1_std, step_error_speed2_std; 
            step_error_accuracy0_std, step_error_accuracy1_std, step_error_accuracy2_std];

    figure;
    hold on;
    categories = {'Visual Perturbation', 'Speed Prompt', 'Accuracy Prompt'};
    x = 1:3; % Groups
    bar_handles = bar(x, means, 'grouped'); 

    if plot_eachparticipant == 1
        step_error_means = {mean(step_error_balance0,2), mean(step_error_balance1,2), mean(step_error_balance2,2); 
                    mean(step_error_speed0,2), mean(step_error_speed1,2), mean(step_error_speed2,2); 
                    mean(step_error_accuracy0,2), mean(step_error_accuracy1,2), mean(step_error_accuracy2,2)};
        
        num_trials = size(step_error_means{1,1}, 1); % Assuming all groups have the same number of trials
        %colors = lines(num_trials); % Use 'lines' colormap to generate distinct colors
        %colors = min(colors * 1.4, 1); % Increase intensity, limit to 1
        for i = 1:size(means, 2) % Loop through balance, speed, accuracy columns
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars
            
            % Scatter plot for individual points with unique colors per row
            for j = 1:3 % Loop through categories (Balance, Speed, Accuracy)
                data_points = step_error_means{j, i}; % Extract row-wise mean values
                for k = 1:num_trials % Loop through trials (rows)
                    scatter(x_offset(j), data_points(k), point_size, colors_participants(k, :), 'filled', 'MarkerFaceAlpha', 0.6);
                end
            end
        end
    else
        for i = 1:size(means, 2)
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars with separate lower and upper limits
        end
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:3); 
    ylim([0 32])
    xticklabels(categories);
    %ylabel('Step Error (cm)');
    ylabel('Mean Value (cm)')
    legend({'Null/Low', 'Medium', 'High'});
    title('Step Error');
    grid on;
    hold off;

    %% Head angle different participants in s3a2b3 by balance and accuracy - Plot 48

    figure;
    hold on;

    balance_levels = {'None', 'Medium', 'High'}; 

    colors_2 = lines(length(bmh_list));  
    %softer_colors = colors * 0.6 + 0.4; 
    softer_colors = colors_participants * 0.3 + 0.7; % Softer colors 
    all_participants_errors = []; 

    for i = 1:length(bmh_list) % For each bmh
        bmh = bmh_list(i);
        s3a2b0 = [];
        s3a2b1 = [];
        s3a2b2 = [];

        labels_head_angles = string(fieldnames(imported_data.(bmh).head_angles));
        for j = 3:29
            if double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 && double(imported_data.(bmh).prompt.balance(j,2)) == 0
                s3a2b0 = [s3a2b0, mean(real(imported_data.(bmh).head_angles.(labels_head_angles(j-2))))];
            elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 && double(imported_data.(bmh).prompt.balance(j,2)) == 1
                s3a2b1 = [s3a2b1, mean(real(imported_data.(bmh).head_angles.(labels_head_angles(j-2))))];
            elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 && double(imported_data.(bmh).prompt.balance(j,2)) == 2
                s3a2b2 = [s3a2b2, mean(real(imported_data.(bmh).head_angles.(labels_head_angles(j-2))))];
            end
        end

        errors_std = [std(s3a2b0), std(s3a2b1), std(s3a2b2)];
        errors_participant = [mean(s3a2b0), mean(s3a2b1), mean(s3a2b2)];
        errors_participant = reshape(errors_participant, 1, []); 
        all_participants_errors = [all_participants_errors; errors_participant];
        errorbar(1:3, errors_participant, errors_std, '--o', 'LineWidth', 2.5, 'MarkerSize', 5, 'Color', softer_colors(i, :), 'MarkerFaceColor', softer_colors(i, :)); 
    end

    % Average calculation
    average_errors = mean(all_participants_errors, 1);
    average_std = std(all_participants_errors, 0, 1);
    red_line = errorbar(1:3, average_errors, average_std, '--o', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0]); % Red average line
    bmh_list = [bmh_list; "Average"]; % Add average for the legend

    legend(red_line, 'Average', 'FontSize', 14); % Adjust the font size as needed

    % No baselines because here we are sorting by balance, not speed
    title("Head Angle for all Participants for High Accuracy Condition")

    set(gca, 'FontSize', 22);
    xlabel('Visual Perturbation');
    ylabel('Head Angle (deg)');
    xticks(1:3);
    xticklabels(balance_levels);
    bmh_list(strcmp(bmh_list, "Average")) = []; % Remove average to bmh_list
    grid on;
    hold off;

    %% Head angle by accuracy for all participants - Plot 49

    figure;
    hold on;

    accuracy_levels = {'Null', 'Medium', 'High'}; 

    colors_2 = lines(length(bmh_list));  
    %softer_colors = colors * 0.6 + 0.4; 
    softer_colors = colors_participants * 0.3 + 0.7; % Softer colors 
    all_participants_errors = []; 

    for i = 1:length(bmh_list) % For each bmh
        bmh = bmh_list(i);
        s3a0b3 = [];
        s3a1b3 = [];
        s3a2b3 = [];

        labels_head_angles = string(fieldnames(imported_data.(bmh).head_angles));
        for j = 3:29
            if double(imported_data.(bmh).prompt.accuracy(j,2)) == 0 
                s3a0b3 = [s3a0b3, mean(real(imported_data.(bmh).head_angles.(labels_head_angles(j-2))))];
            elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 1 
                s3a1b3 = [s3a1b3, mean(real(imported_data.(bmh).head_angles.(labels_head_angles(j-2))))];
            elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 
                s3a2b3 = [s3a2b3, mean(real(imported_data.(bmh).head_angles.(labels_head_angles(j-2))))];
            end
        end

        errors_std = [std(s3a0b3), std(s3a1b3), std(s3a2b3)];
        errors_participant = [mean(s3a0b3), mean(s3a1b3), mean(s3a2b3)];
        errors_participant = reshape(errors_participant, 1, []); 
        all_participants_errors = [all_participants_errors; errors_participant];
        errorbar(1:3, errors_participant, errors_std, '--o', 'LineWidth', 2.5, 'MarkerSize', 5, 'Color', softer_colors(i, :), 'MarkerFaceColor', softer_colors(i, :)); 
    end

    % Average calculation
    average_errors = mean(all_participants_errors, 1);
    average_std = std(all_participants_errors, 0, 1);
    red_line = errorbar(1:3, average_errors, average_std, '--o', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0]); % Red average line
    bmh_list = [bmh_list; "Average"]; % Add average for the legend

    legend(red_line, 'Average', 'FontSize', 14); % Adjust the font size as needed

    % No baselines because here we are sorting by balance, not speed
    title("Head Angle for all Participants")

    set(gca, 'FontSize', 22);
    xlabel('Accuracy Prompt');
    ylabel('Head Angle (deg)');
    xticks(1:3);
    xticklabels(accuracy_levels);
    bmh_list(strcmp(bmh_list, "Average")) = []; % Remove average to bmh_list
    grid on;
    hold off;

    %% Look-ahead distance by accuracy for all participants - Plot 50

    figure;
    hold on;

    accuracy_levels = {'Null', 'Medium', 'High'}; 

    colors_2 = lines(length(bmh_list));  
    %softer_colors = colors * 0.6 + 0.4; 
    softer_colors = colors_participants * 0.3 + 0.7; % Softer colors 
    all_participants_errors = []; 

    for i = 1:length(bmh_list) % For each bmh
        bmh = bmh_list(i);
        s3a0b3 = [];
        s3a1b3 = [];
        s3a2b3 = [];

        labels_gaze = string(fieldnames(imported_data.(bmh).head_angles));
        for j = 3:29
            value = imported_data.(bmh).dist_gaze.(labels_gaze(j-2))/100;
            if value > 4
                value = 4;
            end
            if double(imported_data.(bmh).prompt.accuracy(j,2)) == 0 
                s3a0b3 = [s3a0b3, value];
            elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 1 
                s3a1b3 = [s3a1b3, value];
            elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 
                s3a2b3 = [s3a2b3, value];
            end
        end

        errors_std = [std(s3a0b3), std(s3a1b3), std(s3a2b3)];
        errors_participant = [mean(s3a0b3), mean(s3a1b3), mean(s3a2b3)];
        errors_participant = reshape(errors_participant, 1, []); 
        all_participants_errors = [all_participants_errors; errors_participant];
        errorbar(1:3, errors_participant, errors_std, '--o', 'LineWidth', 2.5, 'MarkerSize', 5, 'Color', softer_colors(i, :), 'MarkerFaceColor', softer_colors(i, :)); 
    end

    % Average calculation
    average_errors = mean(all_participants_errors, 1);
    average_std = std(all_participants_errors, 0, 1);
    red_line = errorbar(1:3, average_errors, average_std, '--o', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0]); % Red average line
    bmh_list = [bmh_list; "Average"]; % Add average for the legend

    legend(red_line, 'Average', 'FontSize', 14); % Adjust the font size as needed

    % No baselines because here we are sorting by balance, not speed
    title("Look-Ahead Distance")

    set(gca, 'FontSize', 22);
    xlabel('Accuracy Prompt');
    ylabel('Look-Ahead Distance (m)');
    xticks(1:3);
    xticklabels(accuracy_levels);
    bmh_list(strcmp(bmh_list, "Average")) = []; % Remove average to bmh_list
    grid on;
    hold off;

    %% Trunk angle different participants in s3a2b3 by balance and accuracy - Plot 51

    figure;
    hold on;

    balance_levels = {'None', 'Medium', 'High'}; 

    colors_2 = lines(length(bmh_list));  
    %softer_colors = colors * 0.6 + 0.4; 
    softer_colors = colors_participants * 0.3 + 0.7; % Softer colors 
    all_participants_errors = []; 

    for i = 1:length(bmh_list) % For each bmh
        bmh = bmh_list(i);
        s3a2b0 = [];
        s3a2b1 = [];
        s3a2b2 = [];

        labels_trunk_angles = string(fieldnames(imported_data.(bmh).trunk_angles));
        for j = 3:29
            if double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 && double(imported_data.(bmh).prompt.balance(j,2)) == 0
                s3a2b0 = [s3a2b0, mean(real(imported_data.(bmh).trunk_angles.(labels_trunk_angles(j-2))))];
            elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 && double(imported_data.(bmh).prompt.balance(j,2)) == 1
                s3a2b1 = [s3a2b1, mean(real(imported_data.(bmh).trunk_angles.(labels_trunk_angles(j-2))))];
            elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 && double(imported_data.(bmh).prompt.balance(j,2)) == 2
                s3a2b2 = [s3a2b2, mean(real(imported_data.(bmh).trunk_angles.(labels_trunk_angles(j-2))))];
            end
        end

        errors_std = [std(s3a2b0), std(s3a2b1), std(s3a2b2)];
        errors_participant = [mean(s3a2b0), mean(s3a2b1), mean(s3a2b2)];
        errors_participant = reshape(errors_participant, 1, []); 
        all_participants_errors = [all_participants_errors; errors_participant];
        errorbar(1:3, errors_participant, errors_std, '--o', 'LineWidth', 2.5, 'MarkerSize', 5, 'Color', softer_colors(i, :), 'MarkerFaceColor', softer_colors(i, :)); 
    end

    % Average calculation
    average_errors = mean(all_participants_errors, 1);
    average_std = std(all_participants_errors, 0, 1);
    red_line = errorbar(1:3, average_errors, average_std, '--o', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0]); % Red average line
    bmh_list = [bmh_list; "Average"]; % Add average for the legend

    legend(red_line, 'Average', 'FontSize', 14); % Adjust the font size as needed

    % No baselines because here we are sorting by balance, not speed
    title("Trunk Angle for all Participants for High Accuracy Condition")

    set(gca, 'FontSize', 22);
    xlabel('Visual Perturbation');
    ylabel('Head Angle (deg)');
    xticks(1:3);
    xticklabels(balance_levels);
    bmh_list(strcmp(bmh_list, "Average")) = []; % Remove average to bmh_list
    grid on;
    hold off;

    %% Trunk angle by accuracy for all participants - Plot 52

    figure;
    hold on;

    accuracy_levels = {'Null', 'Medium', 'High'}; 

    colors_2 = lines(length(bmh_list));  
    %softer_colors = colors * 0.6 + 0.4; 
    softer_colors = colors_participants * 0.3 + 0.7; % Softer colors 
    all_participants_errors = []; 

    for i = 1:length(bmh_list) % For each bmh
        bmh = bmh_list(i);
        s3a0b3 = [];
        s3a1b3 = [];
        s3a2b3 = [];

        labels_trunk_angles = string(fieldnames(imported_data.(bmh).trunk_angles));
        for j = 3:29
            if double(imported_data.(bmh).prompt.accuracy(j,2)) == 0 
                s3a0b3 = [s3a0b3, mean(real(imported_data.(bmh).trunk_angles.(labels_trunk_angles(j-2))))];
            elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 1 
                s3a1b3 = [s3a1b3, mean(real(imported_data.(bmh).trunk_angles.(labels_trunk_angles(j-2))))];
            elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 
                s3a2b3 = [s3a2b3, mean(real(imported_data.(bmh).trunk_angles.(labels_trunk_angles(j-2))))];
            end
        end

        errors_std = [std(s3a0b3), std(s3a1b3), std(s3a2b3)];
        errors_participant = [mean(s3a0b3), mean(s3a1b3), mean(s3a2b3)];
        errors_participant = reshape(errors_participant, 1, []); 
        all_participants_errors = [all_participants_errors; errors_participant];
        errorbar(1:3, errors_participant, errors_std, '--o', 'LineWidth', 2.5, 'MarkerSize', 5, 'Color', softer_colors(i, :), 'MarkerFaceColor', softer_colors(i, :)); 
    end

    % Average calculation
    average_errors = mean(all_participants_errors, 1);
    average_std = std(all_participants_errors, 0, 1);
    red_line = errorbar(1:3, average_errors, average_std, '--o', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0]); % Red average line
    bmh_list = [bmh_list; "Average"]; % Add average for the legend

    legend(red_line, 'Average', 'FontSize', 14); % Adjust the font size as needed

    % No baselines because here we are sorting by balance, not speed
    title("Trunk Angle for all Participants")

    set(gca, 'FontSize', 22);
    xlabel('Accuracy Prompt');
    ylabel('Trunk Angle (deg)');
    xticks(1:3);
    xticklabels(accuracy_levels);
    bmh_list(strcmp(bmh_list, "Average")) = []; % Remove average to bmh_list
    grid on;
    hold off;

    %% Survey balance by condition - Plot 53
    
    survey_balance_balance0 = [];
    survey_balance_balance1 = [];
    survey_balance_balance2 = [];
    survey_balance_speed0 = [];
    survey_balance_speed1 = [];
    survey_balance_speed2 = [];
    survey_balance_accuracy0 = [];
    survey_balance_accuracy1 = [];
    survey_balance_accuracy2 = [];

    for i=1:length(list_labels)
        current_label = list_labels{i};

        if contains(current_label, 's0')
            survey_balance_speed0 = [survey_balance_speed0, data_labels.(current_label).survey_balance(:,1)];
        end
        if contains(current_label, 's1')
            survey_balance_speed1 = [survey_balance_speed1, data_labels.(current_label).survey_balance(:,1)];
        end
        if contains(current_label, 's2')
            survey_balance_speed2 = [survey_balance_speed2, data_labels.(current_label).survey_balance(:,1)];
        end

        if contains(current_label, 'a0')
            survey_balance_accuracy0 = [survey_balance_accuracy0, data_labels.(current_label).survey_balance(:,1)];
        end
        if contains(current_label, 'a1')
            survey_balance_accuracy1 = [survey_balance_accuracy1, data_labels.(current_label).survey_balance(:,1)];
        end
        if contains(current_label, 'a2')
            survey_balance_accuracy2 = [survey_balance_accuracy2, data_labels.(current_label).survey_balance(:,1)];
        end

        if contains(current_label, 'b0')
            survey_balance_balance0 = [survey_balance_balance0, data_labels.(current_label).survey_balance(:,1)];
        end
        if contains(current_label, 'b1')
            survey_balance_balance1 = [survey_balance_balance1, data_labels.(current_label).survey_balance(:,1)];
        end
        if contains(current_label, 'b2')
            survey_balance_balance2 = [survey_balance_balance2, data_labels.(current_label).survey_balance(:,1)];
        end
    end

    survey_balance_balance0_std = std(mean(survey_balance_balance0,2));
    survey_balance_balance1_std = std(mean(survey_balance_balance1,2));
    survey_balance_balance2_std = std(mean(survey_balance_balance2,2));
    survey_balance_speed0_std = std(mean(survey_balance_speed0,2));
    survey_balance_speed1_std = std(mean(survey_balance_speed1,2));
    survey_balance_speed2_std = std(mean(survey_balance_speed2,2));
    survey_balance_accuracy0_std = std(mean(survey_balance_accuracy0,2));
    survey_balance_accuracy1_std = std(mean(survey_balance_accuracy1,2));
    survey_balance_accuracy2_std = std(mean(survey_balance_accuracy2,2));

    survey_balance_balance0_avg = mean(survey_balance_balance0(:));
    survey_balance_balance1_avg = mean(survey_balance_balance1(:));
    survey_balance_balance2_avg = mean(survey_balance_balance2(:));
    survey_balance_speed0_avg = mean(survey_balance_speed0(:));
    survey_balance_speed1_avg = mean(survey_balance_speed1(:));
    survey_balance_speed2_avg = mean(survey_balance_speed2(:));
    survey_balance_accuracy0_avg = mean(survey_balance_accuracy0(:));
    survey_balance_accuracy1_avg = mean(survey_balance_accuracy1(:));
    survey_balance_accuracy2_avg = mean(survey_balance_accuracy2(:));

    means = [survey_balance_balance0_avg, survey_balance_balance1_avg, survey_balance_balance2_avg; 
             survey_balance_speed0_avg, survey_balance_speed1_avg, survey_balance_speed2_avg; 
             survey_balance_accuracy0_avg, survey_balance_accuracy1_avg, survey_balance_accuracy2_avg];
    
    stds = [survey_balance_balance0_std, survey_balance_balance1_std, survey_balance_balance2_std; 
            survey_balance_speed0_std, survey_balance_speed1_std, survey_balance_speed2_std; 
            survey_balance_accuracy0_std, survey_balance_accuracy1_std, survey_balance_accuracy2_std];
        
    figure;
    hold on;
    categories = {'Visual Perturbation', 'Speed Prompt', 'Accuracy Prompt'};
    x = 1:3; % Groups
    bar_handles = bar(x, means, 'grouped'); 

    if plot_eachparticipant == 1
        points_means = {mean(survey_balance_balance0,2), mean(survey_balance_balance1,2), mean(survey_balance_balance2,2); 
                    mean(survey_balance_speed0,2), mean(survey_balance_speed1,2), mean(survey_balance_speed2,2); 
                    mean(survey_balance_accuracy0,2), mean(survey_balance_accuracy1,2), mean(survey_balance_accuracy2,2)};
        
        num_trials = size(points_means{1,1}, 1); % Assuming all groups have the same number of trials
        %colors = lines(num_trials); % Use 'lines' colormap to generate distinct colors
        %colors = min(colors * 1.4, 1); % Increase intensity, limit to 1
        for i = 1:size(means, 2) % Loop through balance, speed, accuracy columns
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars
            
            % Scatter plot for individual points with unique colors per row
            for j = 1:3 % Loop through categories (Balance, Speed, Accuracy)
                data_points = points_means{j, i}; % Extract row-wise mean values
                for k = 1:num_trials % Loop through trials (rows)
                    scatter(x_offset(j), data_points(k), point_size, colors_participants(k, :), 'filled', 'MarkerFaceAlpha', 0.6);
                end
            end
        end
    else
        for i = 1:size(means, 2) % Loop through each group
            % Error bars
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            errorbar(x_offset, means(:,i), stds(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1.5);
        end
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:3); 
    ylim([0 0.46])
    xticklabels(categories);
    ylabel('Normalized Priority');
    legend({'Null/Low', 'Medium', 'High'});
    title('Perceived Balance Priority');
    grid on;
    hold off;

    %% Survey foot placement by condition - Plot 54
    
    survey_foot_placement_balance0 = [];
    survey_foot_placement_balance1 = [];
    survey_foot_placement_balance2 = [];
    survey_foot_placement_speed0 = [];
    survey_foot_placement_speed1 = [];
    survey_foot_placement_speed2 = [];
    survey_foot_placement_accuracy0 = [];
    survey_foot_placement_accuracy1 = [];
    survey_foot_placement_accuracy2 = [];

    for i=1:length(list_labels)
        current_label = list_labels{i};

        if contains(current_label, 's0')
            survey_foot_placement_speed0 = [survey_foot_placement_speed0, data_labels.(current_label).survey_foot_placement(:,1)];
        end
        if contains(current_label, 's1')
            survey_foot_placement_speed1 = [survey_foot_placement_speed1, data_labels.(current_label).survey_foot_placement(:,1)];
        end
        if contains(current_label, 's2')
            survey_foot_placement_speed2 = [survey_foot_placement_speed2, data_labels.(current_label).survey_foot_placement(:,1)];
        end

        if contains(current_label, 'a0')
            survey_foot_placement_accuracy0 = [survey_foot_placement_accuracy0, data_labels.(current_label).survey_foot_placement(:,1)];
        end
        if contains(current_label, 'a1')
            survey_foot_placement_accuracy1 = [survey_foot_placement_accuracy1, data_labels.(current_label).survey_foot_placement(:,1)];
        end
        if contains(current_label, 'a2')
            survey_foot_placement_accuracy2 = [survey_foot_placement_accuracy2, data_labels.(current_label).survey_foot_placement(:,1)];
        end

        if contains(current_label, 'b0')
            survey_foot_placement_balance0 = [survey_foot_placement_balance0, data_labels.(current_label).survey_foot_placement(:,1)];
        end
        if contains(current_label, 'b1')
            survey_foot_placement_balance1 = [survey_foot_placement_balance1, data_labels.(current_label).survey_foot_placement(:,1)];
        end
        if contains(current_label, 'b2')
            survey_foot_placement_balance2 = [survey_foot_placement_balance2, data_labels.(current_label).survey_foot_placement(:,1)];
        end
    end

    survey_foot_placement_balance0_std = std(mean(survey_foot_placement_balance0,2));
    survey_foot_placement_balance1_std = std(mean(survey_foot_placement_balance1,2));
    survey_foot_placement_balance2_std = std(mean(survey_foot_placement_balance2,2));
    survey_foot_placement_speed0_std = std(mean(survey_foot_placement_speed0,2));
    survey_foot_placement_speed1_std = std(mean(survey_foot_placement_speed1,2));
    survey_foot_placement_speed2_std = std(mean(survey_foot_placement_speed2,2));
    survey_foot_placement_accuracy0_std = std(mean(survey_foot_placement_accuracy0,2));
    survey_foot_placement_accuracy1_std = std(mean(survey_foot_placement_accuracy1,2));
    survey_foot_placement_accuracy2_std = std(mean(survey_foot_placement_accuracy2,2));

    survey_foot_placement_balance0_avg = mean(survey_foot_placement_balance0(:));
    survey_foot_placement_balance1_avg = mean(survey_foot_placement_balance1(:));
    survey_foot_placement_balance2_avg = mean(survey_foot_placement_balance2(:));
    survey_foot_placement_speed0_avg = mean(survey_foot_placement_speed0(:));
    survey_foot_placement_speed1_avg = mean(survey_foot_placement_speed1(:));
    survey_foot_placement_speed2_avg = mean(survey_foot_placement_speed2(:));
    survey_foot_placement_accuracy0_avg = mean(survey_foot_placement_accuracy0(:));
    survey_foot_placement_accuracy1_avg = mean(survey_foot_placement_accuracy1(:));
    survey_foot_placement_accuracy2_avg = mean(survey_foot_placement_accuracy2(:));

    means = [survey_foot_placement_balance0_avg, survey_foot_placement_balance1_avg, survey_foot_placement_balance2_avg; 
             survey_foot_placement_speed0_avg, survey_foot_placement_speed1_avg, survey_foot_placement_speed2_avg; 
             survey_foot_placement_accuracy0_avg, survey_foot_placement_accuracy1_avg, survey_foot_placement_accuracy2_avg];
    
    stds = [survey_foot_placement_balance0_std, survey_foot_placement_balance1_std, survey_foot_placement_balance2_std; 
            survey_foot_placement_speed0_std, survey_foot_placement_speed1_std, survey_foot_placement_speed2_std; 
            survey_foot_placement_accuracy0_std, survey_foot_placement_accuracy1_std, survey_foot_placement_accuracy2_std];
        
    figure;
    hold on;
    categories = {'Visual Perturbation', 'Speed Prompt', 'Accuracy Prompt'};
    x = 1:3; % Groups
    bar_handles = bar(x, means, 'grouped'); 

    if plot_eachparticipant == 1
        points_means = {mean(survey_foot_placement_balance0,2), mean(survey_foot_placement_balance1,2), mean(survey_foot_placement_balance2,2); 
                    mean(survey_foot_placement_speed0,2), mean(survey_foot_placement_speed1,2), mean(survey_foot_placement_speed2,2); 
                    mean(survey_foot_placement_accuracy0,2), mean(survey_foot_placement_accuracy1,2), mean(survey_foot_placement_accuracy2,2)};
        
        num_trials = size(points_means{1,1}, 1); % Assuming all groups have the same number of trials
        %colors = lines(num_trials); % Use 'lines' colormap to generate distinct colors
        %colors = min(colors * 1.4, 1); % Increase intensity, limit to 1
        for i = 1:size(means, 2) % Loop through balance, speed, accuracy columns
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars
            
            % Scatter plot for individual points with unique colors per row
            for j = 1:3 % Loop through categories (Balance, Speed, Accuracy)
                data_points = points_means{j, i}; % Extract row-wise mean values
                for k = 1:num_trials % Loop through trials (rows)
                    scatter(x_offset(j), data_points(k), point_size, colors_participants(k, :), 'filled', 'MarkerFaceAlpha', 0.6);
                end
            end
        end
    else
        for i = 1:size(means, 2) % Loop through each group
            % Error bars
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            errorbar(x_offset, means(:,i), stds(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1.5);
        end
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:3); 
    ylim([0 0.52])
    xticklabels(categories);
    ylabel('Normalized Priority');
    legend({'Null/Low', 'Medium', 'High'});
    title('Perceived Foot Placement Priority');
    grid on;
    hold off;

    %% Survey walking speed by condition - Plot 55
    
    survey_walking_speed_balance0 = [];
    survey_walking_speed_balance1 = [];
    survey_walking_speed_balance2 = [];
    survey_walking_speed_speed0 = [];
    survey_walking_speed_speed1 = [];
    survey_walking_speed_speed2 = [];
    survey_walking_speed_accuracy0 = [];
    survey_walking_speed_accuracy1 = [];
    survey_walking_speed_accuracy2 = [];

    for i=1:length(list_labels)
        current_label = list_labels{i};

        if contains(current_label, 's0')
            survey_walking_speed_speed0 = [survey_walking_speed_speed0, data_labels.(current_label).survey_walking_speed(:,1)];
        end
        if contains(current_label, 's1')
            survey_walking_speed_speed1 = [survey_walking_speed_speed1, data_labels.(current_label).survey_walking_speed(:,1)];
        end
        if contains(current_label, 's2')
            survey_walking_speed_speed2 = [survey_walking_speed_speed2, data_labels.(current_label).survey_walking_speed(:,1)];
        end

        if contains(current_label, 'a0')
            survey_walking_speed_accuracy0 = [survey_walking_speed_accuracy0, data_labels.(current_label).survey_walking_speed(:,1)];
        end
        if contains(current_label, 'a1')
            survey_walking_speed_accuracy1 = [survey_walking_speed_accuracy1, data_labels.(current_label).survey_walking_speed(:,1)];
        end
        if contains(current_label, 'a2')
            survey_walking_speed_accuracy2 = [survey_walking_speed_accuracy2, data_labels.(current_label).survey_walking_speed(:,1)];
        end

        if contains(current_label, 'b0')
            survey_walking_speed_balance0 = [survey_walking_speed_balance0, data_labels.(current_label).survey_walking_speed(:,1)];
        end
        if contains(current_label, 'b1')
            survey_walking_speed_balance1 = [survey_walking_speed_balance1, data_labels.(current_label).survey_walking_speed(:,1)];
        end
        if contains(current_label, 'b2')
            survey_walking_speed_balance2 = [survey_walking_speed_balance2, data_labels.(current_label).survey_walking_speed(:,1)];
        end
    end

    survey_walking_speed_balance0_std = std(mean(survey_walking_speed_balance0,2));
    survey_walking_speed_balance1_std = std(mean(survey_walking_speed_balance1,2));
    survey_walking_speed_balance2_std = std(mean(survey_walking_speed_balance2,2));
    survey_walking_speed_speed0_std = std(mean(survey_walking_speed_speed0,2));
    survey_walking_speed_speed1_std = std(mean(survey_walking_speed_speed1,2));
    survey_walking_speed_speed2_std = std(mean(survey_walking_speed_speed2,2));
    survey_walking_speed_accuracy0_std = std(mean(survey_walking_speed_accuracy0,2));
    survey_walking_speed_accuracy1_std = std(mean(survey_walking_speed_accuracy1,2));
    survey_walking_speed_accuracy2_std = std(mean(survey_walking_speed_accuracy2,2));

    survey_walking_speed_balance0_avg = mean(survey_walking_speed_balance0(:));
    survey_walking_speed_balance1_avg = mean(survey_walking_speed_balance1(:));
    survey_walking_speed_balance2_avg = mean(survey_walking_speed_balance2(:));
    survey_walking_speed_speed0_avg = mean(survey_walking_speed_speed0(:));
    survey_walking_speed_speed1_avg = mean(survey_walking_speed_speed1(:));
    survey_walking_speed_speed2_avg = mean(survey_walking_speed_speed2(:));
    survey_walking_speed_accuracy0_avg = mean(survey_walking_speed_accuracy0(:));
    survey_walking_speed_accuracy1_avg = mean(survey_walking_speed_accuracy1(:));
    survey_walking_speed_accuracy2_avg = mean(survey_walking_speed_accuracy2(:));

    means = [survey_walking_speed_balance0_avg, survey_walking_speed_balance1_avg, survey_walking_speed_balance2_avg; 
             survey_walking_speed_speed0_avg, survey_walking_speed_speed1_avg, survey_walking_speed_speed2_avg; 
             survey_walking_speed_accuracy0_avg, survey_walking_speed_accuracy1_avg, survey_walking_speed_accuracy2_avg];
    
    stds = [survey_walking_speed_balance0_std, survey_walking_speed_balance1_std, survey_walking_speed_balance2_std; 
            survey_walking_speed_speed0_std, survey_walking_speed_speed1_std, survey_walking_speed_speed2_std; 
            survey_walking_speed_accuracy0_std, survey_walking_speed_accuracy1_std, survey_walking_speed_accuracy2_std];
        
    figure;
    hold on;
    categories = {'Visual Perturbation', 'Speed Prompt', 'Accuracy Prompt'};
    x = 1:3; % Groups
    bar_handles = bar(x, means, 'grouped'); 

    if plot_eachparticipant == 1
        points_means = {mean(survey_walking_speed_balance0,2), mean(survey_walking_speed_balance1,2), mean(survey_walking_speed_balance2,2); 
                    mean(survey_walking_speed_speed0,2), mean(survey_walking_speed_speed1,2), mean(survey_walking_speed_speed2,2); 
                    mean(survey_walking_speed_accuracy0,2), mean(survey_walking_speed_accuracy1,2), mean(survey_walking_speed_accuracy2,2)};
        
        num_trials = size(points_means{1,1}, 1); % Assuming all groups have the same number of trials
        %colors = lines(num_trials); % Use 'lines' colormap to generate distinct colors
        %colors = min(colors * 1.4, 1); % Increase intensity, limit to 1
        for i = 1:size(means, 2) % Loop through balance, speed, accuracy columns
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars
            
            % Scatter plot for individual points with unique colors per row
            for j = 1:3 % Loop through categories (Balance, Speed, Accuracy)
                data_points = points_means{j, i}; % Extract row-wise mean values
                for k = 1:num_trials % Loop through trials (rows)
                    scatter(x_offset(j), data_points(k), point_size, colors_participants(k, :), 'filled', 'MarkerFaceAlpha', 0.6);
                end
            end
        end
    else
        for i = 1:size(means, 2) % Loop through each group
            % Error bars
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            errorbar(x_offset, means(:,i), stds(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1.5);
        end
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:3); 
    ylim([0 0.92])
    xticklabels(categories);
    ylabel('Normalized Priority');
    legend({'Null/Low', 'Medium', 'High'});
    title('Perceived Walking Speed Priority');
    grid on;
    hold off;

    %% Survey metabolics by condition - Plot 56
    
    survey_metabolics_balance0 = [];
    survey_metabolics_balance1 = [];
    survey_metabolics_balance2 = [];
    survey_metabolics_speed0 = [];
    survey_metabolics_speed1 = [];
    survey_metabolics_speed2 = [];
    survey_metabolics_accuracy0 = [];
    survey_metabolics_accuracy1 = [];
    survey_metabolics_accuracy2 = [];

    for i=1:length(list_labels)
        current_label = list_labels{i};

        if contains(current_label, 's0')
            survey_metabolics_speed0 = [survey_metabolics_speed0, data_labels.(current_label).survey_metabolics(:,1)];
        end
        if contains(current_label, 's1')
            survey_metabolics_speed1 = [survey_metabolics_speed1, data_labels.(current_label).survey_metabolics(:,1)];
        end
        if contains(current_label, 's2')
            survey_metabolics_speed2 = [survey_metabolics_speed2, data_labels.(current_label).survey_metabolics(:,1)];
        end

        if contains(current_label, 'a0')
            survey_metabolics_accuracy0 = [survey_metabolics_accuracy0, data_labels.(current_label).survey_metabolics(:,1)];
        end
        if contains(current_label, 'a1')
            survey_metabolics_accuracy1 = [survey_metabolics_accuracy1, data_labels.(current_label).survey_metabolics(:,1)];
        end
        if contains(current_label, 'a2')
            survey_metabolics_accuracy2 = [survey_metabolics_accuracy2, data_labels.(current_label).survey_metabolics(:,1)];
        end

        if contains(current_label, 'b0')
            survey_metabolics_balance0 = [survey_metabolics_balance0, data_labels.(current_label).survey_metabolics(:,1)];
        end
        if contains(current_label, 'b1')
            survey_metabolics_balance1 = [survey_metabolics_balance1, data_labels.(current_label).survey_metabolics(:,1)];
        end
        if contains(current_label, 'b2')
            survey_metabolics_balance2 = [survey_metabolics_balance2, data_labels.(current_label).survey_metabolics(:,1)];
        end
    end

    survey_metabolics_balance0_std = std(mean(survey_metabolics_balance0,2));
    survey_metabolics_balance1_std = std(mean(survey_metabolics_balance1,2));
    survey_metabolics_balance2_std = std(mean(survey_metabolics_balance2,2));
    survey_metabolics_speed0_std = std(mean(survey_metabolics_speed0,2));
    survey_metabolics_speed1_std = std(mean(survey_metabolics_speed1,2));
    survey_metabolics_speed2_std = std(mean(survey_metabolics_speed2,2));
    survey_metabolics_accuracy0_std = std(mean(survey_metabolics_accuracy0,2));
    survey_metabolics_accuracy1_std = std(mean(survey_metabolics_accuracy1,2));
    survey_metabolics_accuracy2_std = std(mean(survey_metabolics_accuracy2,2));

    survey_metabolics_balance0_avg = mean(survey_metabolics_balance0(:));
    survey_metabolics_balance1_avg = mean(survey_metabolics_balance1(:));
    survey_metabolics_balance2_avg = mean(survey_metabolics_balance2(:));
    survey_metabolics_speed0_avg = mean(survey_metabolics_speed0(:));
    survey_metabolics_speed1_avg = mean(survey_metabolics_speed1(:));
    survey_metabolics_speed2_avg = mean(survey_metabolics_speed2(:));
    survey_metabolics_accuracy0_avg = mean(survey_metabolics_accuracy0(:));
    survey_metabolics_accuracy1_avg = mean(survey_metabolics_accuracy1(:));
    survey_metabolics_accuracy2_avg = mean(survey_metabolics_accuracy2(:));

    means = [survey_metabolics_balance0_avg, survey_metabolics_balance1_avg, survey_metabolics_balance2_avg; 
             survey_metabolics_speed0_avg, survey_metabolics_speed1_avg, survey_metabolics_speed2_avg; 
             survey_metabolics_accuracy0_avg, survey_metabolics_accuracy1_avg, survey_metabolics_accuracy2_avg];
    
    stds = [survey_metabolics_balance0_std, survey_metabolics_balance1_std, survey_metabolics_balance2_std; 
            survey_metabolics_speed0_std, survey_metabolics_speed1_std, survey_metabolics_speed2_std; 
            survey_metabolics_accuracy0_std, survey_metabolics_accuracy1_std, survey_metabolics_accuracy2_std];
        
    figure;
    hold on;
    categories = {'Visual Perturbation', 'Speed Prompt', 'Accuracy Prompt'};
    x = 1:3; % Groups
    bar_handles = bar(x, means, 'grouped'); 

    if plot_eachparticipant == 1
        points_means = {mean(survey_metabolics_balance0,2), mean(survey_metabolics_balance1,2), mean(survey_metabolics_balance2,2); 
                    mean(survey_metabolics_speed0,2), mean(survey_metabolics_speed1,2), mean(survey_metabolics_speed2,2); 
                    mean(survey_metabolics_accuracy0,2), mean(survey_metabolics_accuracy1,2), mean(survey_metabolics_accuracy2,2)};
        
        num_trials = size(points_means{1,1}, 1); % Assuming all groups have the same number of trials
        %colors = lines(num_trials); % Use 'lines' colormap to generate distinct colors
        %colors = min(colors * 1.4, 1); % Increase intensity, limit to 1
        for i = 1:size(means, 2) % Loop through balance, speed, accuracy columns
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars
            
            % Scatter plot for individual points with unique colors per row
            for j = 1:3 % Loop through categories (Balance, Speed, Accuracy)
                data_points = points_means{j, i}; % Extract row-wise mean values
                for k = 1:num_trials % Loop through trials (rows)
                    scatter(x_offset(j), data_points(k), point_size, colors_participants(k, :), 'filled', 'MarkerFaceAlpha', 0.6);
                end
            end
        end
    else
        for i = 1:size(means, 2) % Loop through each group
            % Error bars
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            errorbar(x_offset, means(:,i), stds(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1.5);
        end
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:3); 
    xticklabels(categories);
    ylabel('Normalized Priority');
    legend({'Null/Low', 'Medium', 'High'});
    title('Perceived Energy Expenditure Reduction Priority');
    grid on;
    hold off;

    %% Survey difficulty by condition - Plot 57
    
    survey_difficulty_balance0 = [];
    survey_difficulty_balance1 = [];
    survey_difficulty_balance2 = [];
    survey_difficulty_speed0 = [];
    survey_difficulty_speed1 = [];
    survey_difficulty_speed2 = [];
    survey_difficulty_accuracy0 = [];
    survey_difficulty_accuracy1 = [];
    survey_difficulty_accuracy2 = [];

    for i=1:length(list_labels)
        current_label = list_labels{i};

        if contains(current_label, 's0')
            survey_difficulty_speed0 = [survey_difficulty_speed0, data_labels.(current_label).survey_difficulty(:,1)];
        end
        if contains(current_label, 's1')
            survey_difficulty_speed1 = [survey_difficulty_speed1, data_labels.(current_label).survey_difficulty(:,1)];
        end
        if contains(current_label, 's2')
            survey_difficulty_speed2 = [survey_difficulty_speed2, data_labels.(current_label).survey_difficulty(:,1)];
        end

        if contains(current_label, 'a0')
            survey_difficulty_accuracy0 = [survey_difficulty_accuracy0, data_labels.(current_label).survey_difficulty(:,1)];
        end
        if contains(current_label, 'a1')
            survey_difficulty_accuracy1 = [survey_difficulty_accuracy1, data_labels.(current_label).survey_difficulty(:,1)];
        end
        if contains(current_label, 'a2')
            survey_difficulty_accuracy2 = [survey_difficulty_accuracy2, data_labels.(current_label).survey_difficulty(:,1)];
        end

        if contains(current_label, 'b0')
            survey_difficulty_balance0 = [survey_difficulty_balance0, data_labels.(current_label).survey_difficulty(:,1)];
        end
        if contains(current_label, 'b1')
            survey_difficulty_balance1 = [survey_difficulty_balance1, data_labels.(current_label).survey_difficulty(:,1)];
        end
        if contains(current_label, 'b2')
            survey_difficulty_balance2 = [survey_difficulty_balance2, data_labels.(current_label).survey_difficulty(:,1)];
        end
    end

    survey_difficulty_balance0_std = std(mean(survey_difficulty_balance0,2));
    survey_difficulty_balance1_std = std(mean(survey_difficulty_balance1,2));
    survey_difficulty_balance2_std = std(mean(survey_difficulty_balance2,2));
    survey_difficulty_speed0_std = std(mean(survey_difficulty_speed0,2));
    survey_difficulty_speed1_std = std(mean(survey_difficulty_speed1,2));
    survey_difficulty_speed2_std = std(mean(survey_difficulty_speed2,2));
    survey_difficulty_accuracy0_std = std(mean(survey_difficulty_accuracy0,2));
    survey_difficulty_accuracy1_std = std(mean(survey_difficulty_accuracy1,2));
    survey_difficulty_accuracy2_std = std(mean(survey_difficulty_accuracy2,2));

    survey_difficulty_balance0_avg = mean(survey_difficulty_balance0(:));
    survey_difficulty_balance1_avg = mean(survey_difficulty_balance1(:));
    survey_difficulty_balance2_avg = mean(survey_difficulty_balance2(:));
    survey_difficulty_speed0_avg = mean(survey_difficulty_speed0(:));
    survey_difficulty_speed1_avg = mean(survey_difficulty_speed1(:));
    survey_difficulty_speed2_avg = mean(survey_difficulty_speed2(:));
    survey_difficulty_accuracy0_avg = mean(survey_difficulty_accuracy0(:));
    survey_difficulty_accuracy1_avg = mean(survey_difficulty_accuracy1(:));
    survey_difficulty_accuracy2_avg = mean(survey_difficulty_accuracy2(:));

    means = [survey_difficulty_balance0_avg, survey_difficulty_balance1_avg, survey_difficulty_balance2_avg; 
             survey_difficulty_speed0_avg, survey_difficulty_speed1_avg, survey_difficulty_speed2_avg; 
             survey_difficulty_accuracy0_avg, survey_difficulty_accuracy1_avg, survey_difficulty_accuracy2_avg];
    
    stds = [survey_difficulty_balance0_std, survey_difficulty_balance1_std, survey_difficulty_balance2_std; 
            survey_difficulty_speed0_std, survey_difficulty_speed1_std, survey_difficulty_speed2_std; 
            survey_difficulty_accuracy0_std, survey_difficulty_accuracy1_std, survey_difficulty_accuracy2_std];
        
    figure;
    hold on;
    categories = {'Visual Perturbation', 'Speed Prompt', 'Accuracy Prompt'};
    x = 1:3; % Groups
    bar_handles = bar(x, means, 'grouped'); 

    if plot_eachparticipant == 1
        points_means = {mean(survey_difficulty_balance0,2), mean(survey_difficulty_balance1,2), mean(survey_difficulty_balance2,2); 
                    mean(survey_difficulty_speed0,2), mean(survey_difficulty_speed1,2), mean(survey_difficulty_speed2,2); 
                    mean(survey_difficulty_accuracy0,2), mean(survey_difficulty_accuracy1,2), mean(survey_difficulty_accuracy2,2)};
        
        num_trials = size(points_means{1,1}, 1); % Assuming all groups have the same number of trials
        %colors = lines(num_trials); % Use 'lines' colormap to generate distinct colors
        %colors = min(colors * 1.4, 1); % Increase intensity, limit to 1
        for i = 1:size(means, 2) % Loop through balance, speed, accuracy columns
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars
            
            % Scatter plot for individual points with unique colors per row
            for j = 1:3 % Loop through categories (Balance, Speed, Accuracy)
                data_points = points_means{j, i}; % Extract row-wise mean values
                for k = 1:num_trials % Loop through trials (rows)
                    scatter(x_offset(j), data_points(k), point_size, colors_participants(k, :), 'filled', 'MarkerFaceAlpha', 0.6);
                end
            end
        end
    else
        for i = 1:size(means, 2) % Loop through each group
            % Error bars
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            errorbar(x_offset, means(:,i), stds(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1.5);
        end
    end

    difficulty_levels = {'Extremely Easy', 'Somewhat Easy', 'Neither Easy nor Difficult', 'Somewhat Difficult', 'Extremely Difficult'};
    yticks(1:5);
    yticklabels(difficulty_levels);
    ylim([1 5.1])
    
    set(gca, 'FontSize', 22);
    xticks(1:3); 
    xticklabels(categories);
    %ylabel('Perceived Difficulty');
    legend({'Null/Low', 'Medium', 'High'});
    title('Perceived Difficulty');
    grid on;
    hold off;

    %% All surveys but difficulty in each bar - Plot 58

    values = [survey_balance_balance0_avg, survey_walking_speed_balance0_avg, survey_foot_placement_balance0_avg, survey_metabolics_balance0_avg;
              survey_balance_balance1_avg, survey_walking_speed_balance1_avg, survey_foot_placement_balance1_avg, survey_metabolics_balance1_avg;
              survey_balance_balance2_avg, survey_walking_speed_balance2_avg, survey_foot_placement_balance2_avg, survey_metabolics_balance2_avg;
              NaN, NaN, NaN, NaN;  % Spacer
              survey_balance_speed0_avg, survey_walking_speed_speed0_avg, survey_foot_placement_speed0_avg, survey_metabolics_speed0_avg;
              survey_balance_speed1_avg, survey_walking_speed_speed1_avg, survey_foot_placement_speed1_avg, survey_metabolics_speed1_avg;
              survey_balance_speed2_avg, survey_walking_speed_speed2_avg, survey_foot_placement_speed2_avg, survey_metabolics_speed2_avg;
              NaN, NaN, NaN, NaN;  % Spacer
              survey_balance_accuracy0_avg, survey_walking_speed_accuracy0_avg, survey_foot_placement_accuracy0_avg, survey_metabolics_accuracy0_avg;
              survey_balance_accuracy1_avg, survey_walking_speed_accuracy1_avg, survey_foot_placement_accuracy1_avg, survey_metabolics_accuracy1_avg;
              survey_balance_accuracy2_avg, survey_walking_speed_accuracy2_avg, survey_foot_placement_accuracy2_avg, survey_metabolics_accuracy2_avg];

    colors_2 = [0.2 0.4 0.8;   % Soft Blue  
          0.8 0.3 0.3;   % Soft Red  
          0.3 0.7 0.3;   % Soft Green  
          0.6 0.6 0.6];  % Medium Gray 

    figure;
    b = bar(values, 'stacked');
    for i = 1:length(b)
        b(i).FaceColor = colors_2(i, :);
    end
    
    xticks([1 2 3 5 6 7 9 10 11]); % Ignore the NaN rows
    xticklabels({'No Perturbation', 'Medium Perturbation', 'High Perturbation', ...
                 'Slow Speed', 'Medium Speed', 'High Speed', ...
                 'No Accuracy', 'Medium Accuracy', 'High Accuracy'});
    
    ylim([0 1.25]);       
    yticks(0:0.2:1);    
    set(gca, 'FontSize', 22);
    ylabel('Normalized Priority');
    title('Perceived Prioritization by Category');
    legend({'Balance', 'Walking Speed', 'Foot Placement', 'Energy Expenditure Reduction'}, 'Location', 'northeast', 'FontSize', 14);


    %% Step error by accuracy for all participants - Plot 59

    figure
    hold on
    
    accuracy_levels = {'Null', 'Medium', 'High'};  
    colors_2 = lines(length(bmh_list));  
    %softer_colors = colors * 0.6 + 0.4; 
    softer_colors = colors_participants * 0.3 + 0.7; % Softer colors
    all_participants_errors = []; 

    legend_handles = []; 
    
    for i = 1:length(bmh_list) % One iteration for each BMH
        bmh = bmh_list(i);
        s3a0b3 = [];
        s3a1b3 = [];
        s3a2b3 = [];

        for j = 3:29
            if double(imported_data.(bmh).prompt.accuracy(j,2)) == 0
                s3a0b3 = [s3a0b3, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10]; % From mm to cm
            elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 1
                s3a1b3 = [s3a1b3, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10]; % From mm to cm
            elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 2
                s3a2b3 = [s3a2b3, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10]; % From mm to cm
            end
        end
        
        errors_std = [std(s3a0b3), std(s3a1b3), std(s3a2b3)]; % std calculation before subtracting the baseline(s)
        errors_participant = [mean(s3a0b3), mean(s3a1b3), mean(s3a2b3)];
        all_participants_errors = [all_participants_errors; errors_participant];        
        if shaded_plot == 0
            errorbar(1:3, errors_participant, errors_std, '--o', 'LineWidth', 2.5, 'MarkerSize', 5, 'Color', softer_colors(i, :), 'MarkerFaceColor', softer_colors(i, :)); 
        else
            x = 1:3;
            fill([x, fliplr(x)], [errors_participant + errors_std, fliplr(errors_participant - errors_std)], ...
                 softer_colors(i, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        
            h = plot(1:3, errors_participant, '--o', 'LineWidth', 2.5, 'MarkerSize', 5, ...
                 'Color', softer_colors(i, :), 'MarkerFaceColor', softer_colors(i, :), 'DisplayName', string(bmh)); 
            legend_handles = [legend_handles, h]; 
        end
    end
    
    % Average calculation
    average_errors = mean(all_participants_errors, 1);
    average_std = std(all_participants_errors, 0, 1);
    if shaded_plot == 0
        red_line = errorbar(1:3, average_errors, average_std, '--o', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0]); % Red average line
        bmh_list = [bmh_list; "Average"]; % Add average for the legend
        legend(red_line, 'Average', 'FontSize', 14); % Adjust the font size as needed
        bmh_list(strcmp(bmh_list, "Average")) = []; % Remove average to bmh_list
    else
        h_avg = plot(1:3, average_errors, '--o', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0], 'DisplayName', 'Average');
        legend_handles = [legend_handles, h_avg];
        legend(legend_handles, 'FontSize', 14, 'Location', 'southwest');
    end

    % No baselines because here we are sorting by balance, not speed
    title("Step Error for all Participants")
    
    set(gca, 'FontSize', 22);
    xlabel('Accuracy Prompt');
    ylabel('Step Error (cm)');
    xticks(1:3);
    xticklabels(accuracy_levels);
    grid on;
    hold off

    %% Step error high accuracy across speeds - Plot 60

    figure;
    hold on;

    speed_levels = {'Slow', 'Medium', 'Fast'}; 

    colors_2 = lines(length(bmh_list));  
    %softer_colors = colors * 0.6 + 0.4; 
    softer_colors = colors_participants * 0.3 + 0.7; % Softer colors 
    all_participants_errors = []; 

    for i = 1:length(bmh_list) % For each bmh
        bmh = bmh_list(i);
        s3a2b0 = [];
        s3a2b1 = [];
        s3a2b2 = [];

        for j = 3:29
            if double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 && double(imported_data.(bmh).prompt.speed(j,2)) == 0
                s3a2b0 = [s3a2b0, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10]; %From mm to cm
            elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 && double(imported_data.(bmh).prompt.speed(j,2)) == 1
                s3a2b1 = [s3a2b1, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10]; %From mm to cm
            elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 2 && double(imported_data.(bmh).prompt.speed(j,2)) == 2
                s3a2b2 = [s3a2b2, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10]; %From mm to cm
            end
        end

        errors_std = [std(s3a2b0), std(s3a2b1), std(s3a2b2)];
        errors_participant = [mean(s3a2b0), mean(s3a2b1), mean(s3a2b2)];
        errors_participant = reshape(errors_participant, 1, []); 
        all_participants_errors = [all_participants_errors; errors_participant];
        errorbar(1:3, errors_participant, errors_std, '--o', 'LineWidth', 2.5, 'MarkerSize', 5, 'Color', softer_colors(i, :), 'MarkerFaceColor', softer_colors(i, :)); 
    end

    % Average calculation
    average_errors = mean(all_participants_errors, 1);
    average_std = std(all_participants_errors, 0, 1);
    red_line = errorbar(1:3, average_errors, average_std, '--o', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0]); % Red average line
    bmh_list = [bmh_list; "Average"]; % Add average for the legend

    legend(red_line, 'Average', 'FontSize', 14); % Adjust the font size as needed

    % No baselines because here we are sorting by balance, not speed
    title("Step Error for all Participants for High Accuracy Condition")

    set(gca, 'FontSize', 22);
    xlabel('Speed Prompt');
    ylabel('Step Error (cm)');
    xticks(1:3);
    xticklabels(speed_levels);
    bmh_list(strcmp(bmh_list, "Average")) = []; % Remove average to bmh_list
    grid on;
    hold off;    

    %% Step error medium accuracy across speeds - Plot 61

    figure;
    hold on;

    speed_levels = {'Slow', 'Medium', 'Fast'}; 

    colors_2 = lines(length(bmh_list));  
    %softer_colors = colors * 0.6 + 0.4; 
    softer_colors = colors_participants * 0.3 + 0.7; % Softer colors 
    all_participants_errors = []; 

    for i = 1:length(bmh_list) % For each bmh
        bmh = bmh_list(i);
        s3a2b0 = [];
        s3a2b1 = [];
        s3a2b2 = [];

        for j = 3:29
            if double(imported_data.(bmh).prompt.accuracy(j,2)) == 1 && double(imported_data.(bmh).prompt.speed(j,2)) == 0
                s3a2b0 = [s3a2b0, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10]; %From mm to cm
            elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 1 && double(imported_data.(bmh).prompt.speed(j,2)) == 1
                s3a2b1 = [s3a2b1, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10]; %From mm to cm
            elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 1 && double(imported_data.(bmh).prompt.speed(j,2)) == 2
                s3a2b2 = [s3a2b2, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10]; %From mm to cm
            end
        end

        errors_std = [std(s3a2b0), std(s3a2b1), std(s3a2b2)];
        errors_participant = [mean(s3a2b0), mean(s3a2b1), mean(s3a2b2)];
        errors_participant = reshape(errors_participant, 1, []); 
        all_participants_errors = [all_participants_errors; errors_participant];
        errorbar(1:3, errors_participant, errors_std, '--o', 'LineWidth', 2.5, 'MarkerSize', 5, 'Color', softer_colors(i, :), 'MarkerFaceColor', softer_colors(i, :)); 
    end

    % Average calculation
    average_errors = mean(all_participants_errors, 1);
    average_std = std(all_participants_errors, 0, 1);
    red_line = errorbar(1:3, average_errors, average_std, '--o', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0]); % Red average line
    bmh_list = [bmh_list; "Average"]; % Add average for the legend

    legend(red_line, 'Average', 'FontSize', 14); % Adjust the font size as needed

    % No baselines because here we are sorting by balance, not speed
    title("Step Error for all Participants for Medium Accuracy Condition")

    set(gca, 'FontSize', 22);
    xlabel('Speed Prompt');
    ylabel('Step Error (cm)');
    xticks(1:3);
    xticklabels(speed_levels);
    bmh_list(strcmp(bmh_list, "Average")) = []; % Remove average to bmh_list
    grid on;
    hold off;    

    %% Step error by balance for all participants - Plot 62

    figure
    hold on
    
    balance_levels = {'None', 'Medium', 'High'};  
    colors_2 = lines(length(bmh_list));  
    %softer_colors = colors * 0.6 + 0.4; 
    softer_colors = colors_participants * 0.3 + 0.7; % Softer colors
    all_participants_errors = []; 

    legend_handles = []; 
    
    for i = 1:length(bmh_list) % One iteration for each BMH
        bmh = bmh_list(i);
        s3a0b3 = [];
        s3a1b3 = [];
        s3a2b3 = [];

        for j = 3:29
            if double(imported_data.(bmh).prompt.balance(j,2)) == 0
                s3a0b3 = [s3a0b3, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10]; % From mm to cm
            elseif double(imported_data.(bmh).prompt.balance(j,2)) == 1
                s3a1b3 = [s3a1b3, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10]; % From mm to cm
            elseif double(imported_data.(bmh).prompt.balance(j,2)) == 2
                s3a2b3 = [s3a2b3, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10]; % From mm to cm
            end
        end
        
        errors_std = [std(s3a0b3), std(s3a1b3), std(s3a2b3)]; % std calculation before subtracting the baseline(s)
        errors_participant = [mean(s3a0b3), mean(s3a1b3), mean(s3a2b3)];
        all_participants_errors = [all_participants_errors; errors_participant];        
        if shaded_plot == 0
            errorbar(1:3, errors_participant, errors_std, '--o', 'LineWidth', 2.5, 'MarkerSize', 5, 'Color', softer_colors(i, :), 'MarkerFaceColor', softer_colors(i, :)); 
        else
            x = 1:3;
            fill([x, fliplr(x)], [errors_participant + errors_std, fliplr(errors_participant - errors_std)], ...
                 softer_colors(i, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        
            h = plot(1:3, errors_participant, '--o', 'LineWidth', 2.5, 'MarkerSize', 5, ...
                 'Color', softer_colors(i, :), 'MarkerFaceColor', softer_colors(i, :), 'DisplayName', string(bmh)); 
            legend_handles = [legend_handles, h]; 
        end
    end
    
    % Average calculation
    average_errors = mean(all_participants_errors, 1);
    average_std = std(all_participants_errors, 0, 1);
    if shaded_plot == 0
        red_line = errorbar(1:3, average_errors, average_std, '--o', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0]); % Red average line
        bmh_list = [bmh_list; "Average"]; % Add average for the legend
        %legend(bmh_list, 'FontSize', 14, 'Location', 'southwest'); % Adjust the font size as needed
        legend(red_line, 'Average', 'FontSize', 14); % Adjust the font size as needed
        bmh_list(strcmp(bmh_list, "Average")) = []; % Remove average to bmh_list
    else
        h_avg = plot(1:3, average_errors, '--o', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0], 'DisplayName', 'Average');
        legend_handles = [legend_handles, h_avg];
        legend(legend_handles, 'FontSize', 14, 'Location', 'southwest');
    end

    % No baselines because here we are sorting by balance, not speed
    title("Step Error for all Participants")
    
    set(gca, 'FontSize', 22);
    xlabel('Visual Perturbation');
    ylabel('Step Error (cm)');
    xticks(1:3);
    xticklabels(balance_levels);
    grid on;
    hold off

    %% Step error by speed for all participants - Plot 63

    figure
    hold on
    
    speed_levels = {'Slow', 'Medium', 'Fast'};  
    colors_2 = lines(length(bmh_list));  
    %softer_colors = colors * 0.6 + 0.4; 
    softer_colors = colors_participants * 0.3 + 0.7; % Softer colors
    all_participants_errors = []; 

    legend_handles = []; 
    
    for i = 1:length(bmh_list) % One iteration for each BMH
        bmh = bmh_list(i);
        s3a0b3 = [];
        s3a1b3 = [];
        s3a2b3 = [];

        for j = 3:29
            if double(imported_data.(bmh).prompt.speed(j,2)) == 0
                s3a0b3 = [s3a0b3, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10]; % From mm to cm
            elseif double(imported_data.(bmh).prompt.speed(j,2)) == 1
                s3a1b3 = [s3a1b3, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10]; % From mm to cm
            elseif double(imported_data.(bmh).prompt.speed(j,2)) == 2
                s3a2b3 = [s3a2b3, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10]; % From mm to cm
            end
        end
        
        errors_std = [std(s3a0b3), std(s3a1b3), std(s3a2b3)]; % std calculation before subtracting the baseline(s)
        errors_participant = [mean(s3a0b3), mean(s3a1b3), mean(s3a2b3)];
        all_participants_errors = [all_participants_errors; errors_participant];        
        if shaded_plot == 0
            errorbar(1:3, errors_participant, errors_std, '--o', 'LineWidth', 2.5, 'MarkerSize', 5, 'Color', softer_colors(i, :), 'MarkerFaceColor', softer_colors(i, :)); 
        else
            x = 1:3;
            fill([x, fliplr(x)], [errors_participant + errors_std, fliplr(errors_participant - errors_std)], ...
                 softer_colors(i, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        
            h = plot(1:3, errors_participant, '--o', 'LineWidth', 2.5, 'MarkerSize', 5, ...
                 'Color', softer_colors(i, :), 'MarkerFaceColor', softer_colors(i, :), 'DisplayName', string(bmh)); 
            legend_handles = [legend_handles, h]; 
        end
    end
    
    % Average calculation
    average_errors = mean(all_participants_errors, 1);
    average_std = std(all_participants_errors, 0, 1);
    if shaded_plot == 0
        red_line = errorbar(1:3, average_errors, average_std, '--o', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0]); % Red average line
        bmh_list = [bmh_list; "Average"]; % Add average for the legend
        %legend(bmh_list, 'FontSize', 14, 'Location', 'southwest'); % Adjust the font size as needed
        legend(red_line, 'Average', 'FontSize', 14); % Adjust the font size as needed
        bmh_list(strcmp(bmh_list, "Average")) = []; % Remove average to bmh_list
    else
        h_avg = plot(1:3, average_errors, '--o', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0], 'DisplayName', 'Average');
        legend_handles = [legend_handles, h_avg];
        legend(legend_handles, 'FontSize', 14, 'Location', 'southwest');
    end

    % No baselines because here we are sorting by balance, not speed
    title("Step Error for all Participants")
    
    set(gca, 'FontSize', 22);
    xlabel('Speed Prompt');
    ylabel('Step Error (cm)');
    xticks(1:3);
    xticklabels(speed_levels);
    grid on;
    hold off

    %% Step error medium accuracy across balance - Plot 64

    figure;
    hold on;

    colors_2 = lines(length(bmh_list));  
    %softer_colors = colors * 0.6 + 0.4; 
    softer_colors = colors_participants * 0.3 + 0.7; % Softer colors 
    all_participants_errors = []; 

    for i = 1:length(bmh_list) % For each bmh
        bmh = bmh_list(i);
        s3a2b0 = [];
        s3a2b1 = [];
        s3a2b2 = [];

        for j = 3:29
            if double(imported_data.(bmh).prompt.accuracy(j,2)) == 1 && double(imported_data.(bmh).prompt.balance(j,2)) == 0
                s3a2b0 = [s3a2b0, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10]; %From mm to cm
            elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 1 && double(imported_data.(bmh).prompt.balance(j,2)) == 1
                s3a2b1 = [s3a2b1, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10]; %From mm to cm
            elseif double(imported_data.(bmh).prompt.accuracy(j,2)) == 1 && double(imported_data.(bmh).prompt.balance(j,2)) == 2
                s3a2b2 = [s3a2b2, imported_data.(bmh).overall_errors.straight_mean_absolute(j)/10]; %From mm to cm
            end
        end

        errors_std = [std(s3a2b0), std(s3a2b1), std(s3a2b2)];
        errors_participant = [mean(s3a2b0), mean(s3a2b1), mean(s3a2b2)];
        errors_participant = reshape(errors_participant, 1, []); 
        all_participants_errors = [all_participants_errors; errors_participant];
        errorbar(1:3, errors_participant, errors_std, '--o', 'LineWidth', 2.5, 'MarkerSize', 5, 'Color', softer_colors(i, :), 'MarkerFaceColor', softer_colors(i, :)); 
    end

    % Average calculation
    average_errors = mean(all_participants_errors, 1);
    average_std = std(all_participants_errors, 0, 1);
    red_line = errorbar(1:3, average_errors, average_std, '--o', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [1, 0, 0], 'MarkerFaceColor', [1, 0, 0]); % Red average line
    bmh_list = [bmh_list; "Average"]; % Add average for the legend

    legend(red_line, 'Average', 'FontSize', 14); % Adjust the font size as needed

    % No baselines because here we are sorting by balance, not speed
    title("Step Error for all Participants for Medium Accuracy Condition")

    set(gca, 'FontSize', 22);
    xlabel('Visual Perturbation');
    ylabel('Step Error (cm)');
    xticks(1:3);
    xticklabels(balance_levels);
    bmh_list(strcmp(bmh_list, "Average")) = []; % Remove average to bmh_list
    grid on;
    hold off;    

    %% Signed error by condition - Plot 65

    step_error_signed_balance0 = [];
    step_error_signed_balance1 = [];
    step_error_signed_balance2 = [];
    step_error_signed_speed0 = [];
    step_error_signed_speed1 = [];
    step_error_signed_speed2 = [];
    step_error_signed_accuracy0 = [];
    step_error_signed_accuracy1 = [];
    step_error_signed_accuracy2 = [];

    for i=1:length(list_labels)
        current_label = list_labels{i};

        if contains(current_label, 's0')
            step_error_signed_speed0 = [step_error_signed_speed0, data_labels.(current_label).step_error_signed(:,1)];
        end
        if contains(current_label, 's1')
            step_error_signed_speed1 = [step_error_signed_speed1, data_labels.(current_label).step_error_signed(:,1)];
        end
        if contains(current_label, 's2')
            step_error_signed_speed2 = [step_error_signed_speed2, data_labels.(current_label).step_error_signed(:,1)];
        end

        if contains(current_label, 'a0')
            step_error_signed_accuracy0 = [step_error_signed_accuracy0, data_labels.(current_label).step_error_signed(:,1)];
        end
        if contains(current_label, 'a1')
            step_error_signed_accuracy1 = [step_error_signed_accuracy1, data_labels.(current_label).step_error_signed(:,1)];
        end
        if contains(current_label, 'a2')
            step_error_signed_accuracy2 = [step_error_signed_accuracy2, data_labels.(current_label).step_error_signed(:,1)];
        end

        if contains(current_label, 'b0')
            step_error_signed_balance0 = [step_error_signed_balance0, data_labels.(current_label).step_error_signed(:,1)];
        end
        if contains(current_label, 'b1')
            step_error_signed_balance1 = [step_error_signed_balance1, data_labels.(current_label).step_error_signed(:,1)];
        end
        if contains(current_label, 'b2')
            step_error_signed_balance2 = [step_error_signed_balance2, data_labels.(current_label).step_error_signed(:,1)];
        end
    end

    step_error_signed_balance0_std = std(mean(step_error_signed_balance0,2));
    step_error_signed_balance1_std = std(mean(step_error_signed_balance1,2));
    step_error_signed_balance2_std = std(mean(step_error_signed_balance2,2));
    step_error_signed_speed0_std = std(mean(step_error_signed_speed0,2));
    step_error_signed_speed1_std = std(mean(step_error_signed_speed1,2));
    step_error_signed_speed2_std = std(mean(step_error_signed_speed2,2));
    step_error_signed_accuracy0_std = std(mean(step_error_signed_accuracy0,2));
    step_error_signed_accuracy1_std = std(mean(step_error_signed_accuracy1,2));
    step_error_signed_accuracy2_std = std(mean(step_error_signed_accuracy2,2));

    step_error_signed_balance0_avg = mean(step_error_signed_balance0(:));
    step_error_signed_balance1_avg = mean(step_error_signed_balance1(:));
    step_error_signed_balance2_avg = mean(step_error_signed_balance2(:));
    step_error_signed_speed0_avg = mean(step_error_signed_speed0(:));
    step_error_signed_speed1_avg = mean(step_error_signed_speed1(:));
    step_error_signed_speed2_avg = mean(step_error_signed_speed2(:));
    step_error_signed_accuracy0_avg = mean(step_error_signed_accuracy0(:));
    step_error_signed_accuracy1_avg = mean(step_error_signed_accuracy1(:));
    step_error_signed_accuracy2_avg = mean(step_error_signed_accuracy2(:));

    means = [step_error_signed_balance0_avg, step_error_signed_balance1_avg, step_error_signed_balance2_avg; 
             step_error_signed_speed0_avg, step_error_signed_speed1_avg, step_error_signed_speed2_avg; 
             step_error_signed_accuracy0_avg, step_error_signed_accuracy1_avg, step_error_signed_accuracy2_avg];

    stds = [step_error_signed_balance0_std, step_error_signed_balance1_std, step_error_signed_balance2_std; 
            step_error_signed_speed0_std, step_error_signed_speed1_std, step_error_signed_speed2_std; 
            step_error_signed_accuracy0_std, step_error_signed_accuracy1_std, step_error_signed_accuracy2_std];

    figure;
    hold on;
    categories = {'Visual Perturbation', 'Speed Prompt', 'Accuracy Prompt'};
    x = 1:3; % Groups
    bar_handles = bar(x, means, 'grouped'); 

    if plot_eachparticipant == 1
        step_error_signed_means = {mean(step_error_signed_balance0,2), mean(step_error_signed_balance1,2), mean(step_error_signed_balance2,2); 
                    mean(step_error_signed_speed0,2), mean(step_error_signed_speed1,2), mean(step_error_signed_speed2,2); 
                    mean(step_error_signed_accuracy0,2), mean(step_error_signed_accuracy1,2), mean(step_error_signed_accuracy2,2)};
        
        num_trials = size(step_error_signed_means{1,1}, 1); % Assuming all groups have the same number of trials
        %colors = lines(num_trials); % Use 'lines' colormap to generate distinct colors
        %colors = min(colors * 1.4, 1); % Increase intensity, limit to 1
        for i = 1:size(means, 2) % Loop through balance, speed, accuracy columns
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars
            
            % Scatter plot for individual points with unique colors per row
            for j = 1:3 % Loop through categories (Balance, Speed, Accuracy)
                data_points = step_error_signed_means{j, i}; % Extract row-wise mean values
                for k = 1:num_trials % Loop through trials (rows)
                    scatter(x_offset(j), data_points(k), point_size, colors_participants(k, :), 'filled', 'MarkerFaceAlpha', 0.6);
                end
            end
        end
    else
        for i = 1:size(means, 2)
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars with separate lower and upper limits
        end
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:3); 
    ylim([-17 20])
    xticklabels(categories);
    %ylabel('Step Error - Signed (cm)');
    ylabel('Mean Value (cm)')
    legend({'Null/Low', 'Medium', 'High'});
    title('Signed Step Error');
    grid on;
    hold off;

    %% Signed error var by condition - Plot 66

    step_error_signed_var_balance0 = [];
    step_error_signed_var_balance1 = [];
    step_error_signed_var_balance2 = [];
    step_error_signed_var_speed0 = [];
    step_error_signed_var_speed1 = [];
    step_error_signed_var_speed2 = [];
    step_error_signed_var_accuracy0 = [];
    step_error_signed_var_accuracy1 = [];
    step_error_signed_var_accuracy2 = [];

    for i=1:length(list_labels)
        current_label = list_labels{i};

        if contains(current_label, 's0')
            step_error_signed_var_speed0 = [step_error_signed_var_speed0, data_labels.(current_label).step_error_signed_std(:,1)];
        end
        if contains(current_label, 's1')
            step_error_signed_var_speed1 = [step_error_signed_var_speed1, data_labels.(current_label).step_error_signed_std(:,1)];
        end
        if contains(current_label, 's2')
            step_error_signed_var_speed2 = [step_error_signed_var_speed2, data_labels.(current_label).step_error_signed_std(:,1)];
        end

        if contains(current_label, 'a0')
            step_error_signed_var_accuracy0 = [step_error_signed_var_accuracy0, data_labels.(current_label).step_error_signed_std(:,1)];
        end
        if contains(current_label, 'a1')
            step_error_signed_var_accuracy1 = [step_error_signed_var_accuracy1, data_labels.(current_label).step_error_signed_std(:,1)];
        end
        if contains(current_label, 'a2')
            step_error_signed_var_accuracy2 = [step_error_signed_var_accuracy2, data_labels.(current_label).step_error_signed_std(:,1)];
        end

        if contains(current_label, 'b0')
            step_error_signed_var_balance0 = [step_error_signed_var_balance0, data_labels.(current_label).step_error_signed_std(:,1)];
        end
        if contains(current_label, 'b1')
            step_error_signed_var_balance1 = [step_error_signed_var_balance1, data_labels.(current_label).step_error_signed_std(:,1)];
        end
        if contains(current_label, 'b2')
            step_error_signed_var_balance2 = [step_error_signed_var_balance2, data_labels.(current_label).step_error_signed_std(:,1)];
        end
    end

    step_error_signed_var_balance0_std = std(mean(step_error_signed_var_balance0,2));
    step_error_signed_var_balance1_std = std(mean(step_error_signed_var_balance1,2));
    step_error_signed_var_balance2_std = std(mean(step_error_signed_var_balance2,2));
    step_error_signed_var_speed0_std = std(mean(step_error_signed_var_speed0,2));
    step_error_signed_var_speed1_std = std(mean(step_error_signed_var_speed1,2));
    step_error_signed_var_speed2_std = std(mean(step_error_signed_var_speed2,2));
    step_error_signed_var_accuracy0_std = std(mean(step_error_signed_var_accuracy0,2));
    step_error_signed_var_accuracy1_std = std(mean(step_error_signed_var_accuracy1,2));
    step_error_signed_var_accuracy2_std = std(mean(step_error_signed_var_accuracy2,2));

    step_error_signed_var_balance0_avg = mean(step_error_signed_var_balance0(:));
    step_error_signed_var_balance1_avg = mean(step_error_signed_var_balance1(:));
    step_error_signed_var_balance2_avg = mean(step_error_signed_var_balance2(:));
    step_error_signed_var_speed0_avg = mean(step_error_signed_var_speed0(:));
    step_error_signed_var_speed1_avg = mean(step_error_signed_var_speed1(:));
    step_error_signed_var_speed2_avg = mean(step_error_signed_var_speed2(:));
    step_error_signed_var_accuracy0_avg = mean(step_error_signed_var_accuracy0(:));
    step_error_signed_var_accuracy1_avg = mean(step_error_signed_var_accuracy1(:));
    step_error_signed_var_accuracy2_avg = mean(step_error_signed_var_accuracy2(:));

    means = [step_error_signed_var_balance0_avg, step_error_signed_var_balance1_avg, step_error_signed_var_balance2_avg; 
             step_error_signed_var_speed0_avg, step_error_signed_var_speed1_avg, step_error_signed_var_speed2_avg; 
             step_error_signed_var_accuracy0_avg, step_error_signed_var_accuracy1_avg, step_error_signed_var_accuracy2_avg];

    stds = [step_error_signed_var_balance0_std, step_error_signed_var_balance1_std, step_error_signed_var_balance2_std; 
            step_error_signed_var_speed0_std, step_error_signed_var_speed1_std, step_error_signed_var_speed2_std; 
            step_error_signed_var_accuracy0_std, step_error_signed_var_accuracy1_std, step_error_signed_var_accuracy2_std];

    figure;
    hold on;
    categories = {'Visual Perturbation', 'Speed Prompt', 'Accuracy Prompt'};
    x = 1:3; % Groups
    bar_handles = bar(x, means, 'grouped'); 

    if plot_eachparticipant == 1
        step_error_signed_var_means = {mean(step_error_signed_var_balance0,2), mean(step_error_signed_var_balance1,2), mean(step_error_signed_var_balance2,2); 
                    mean(step_error_signed_var_speed0,2), mean(step_error_signed_var_speed1,2), mean(step_error_signed_var_speed2,2); 
                    mean(step_error_signed_var_accuracy0,2), mean(step_error_signed_var_accuracy1,2), mean(step_error_signed_var_accuracy2,2)};
        
        num_trials = size(step_error_signed_var_means{1,1}, 1); % Assuming all groups have the same number of trials
        %colors = lines(num_trials); % Use 'lines' colormap to generate distinct colors
        %colors = min(colors * 1.4, 1); % Increase intensity, limit to 1
        for i = 1:size(means, 2) % Loop through balance, speed, accuracy columns
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars
            
            % Scatter plot for individual points with unique colors per row
            for j = 1:3 % Loop through categories (Balance, Speed, Accuracy)
                data_points = step_error_signed_var_means{j, i}; % Extract row-wise mean values
                for k = 1:num_trials % Loop through trials (rows)
                    scatter(x_offset(j), data_points(k), point_size, colors_participants(k, :), 'filled', 'MarkerFaceAlpha', 0.6);
                end
            end
        end
    else
        for i = 1:size(means, 2)
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars with separate lower and upper limits
        end
    end
    
    set(gca, 'FontSize', 22);
    ylim([0 35])
    xticks(1:3); 
    xticklabels(categories);
    %ylabel('Signed Step Error Variability (cm)');
    ylabel('Mean Value (cm)')
    legend({'Null/Low', 'Medium', 'High'});
    title('Signed Step Error Variability');
    grid on;
    hold off;

    %% Metabolics normalized - Plot 67

    metabolics_normalized_balance0 = [];
    metabolics_normalized_balance1 = [];
    metabolics_normalized_balance2 = [];
    metabolics_normalized_speed0 = [];
    metabolics_normalized_speed1 = [];
    metabolics_normalized_speed2 = [];
    metabolics_normalized_accuracy0 = [];
    metabolics_normalized_accuracy1 = [];
    metabolics_normalized_accuracy2 = [];

    for i=1:length(list_labels)
        current_label = list_labels{i};

        if contains(current_label, 's0')
            metabolics_normalized_speed0 = [metabolics_normalized_speed0, data_labels.(current_label).metabolics.normalized(:,1)];
        end
        if contains(current_label, 's1')
            metabolics_normalized_speed1 = [metabolics_normalized_speed1, data_labels.(current_label).metabolics.normalized(:,1)];
        end
        if contains(current_label, 's2')
            metabolics_normalized_speed2 = [metabolics_normalized_speed2, data_labels.(current_label).metabolics.normalized(:,1)];
        end

        if contains(current_label, 'a0')
            metabolics_normalized_accuracy0 = [metabolics_normalized_accuracy0, data_labels.(current_label).metabolics.normalized(:,1)];
        end
        if contains(current_label, 'a1')
            metabolics_normalized_accuracy1 = [metabolics_normalized_accuracy1, data_labels.(current_label).metabolics.normalized(:,1)];
        end
        if contains(current_label, 'a2')
            metabolics_normalized_accuracy2 = [metabolics_normalized_accuracy2, data_labels.(current_label).metabolics.normalized(:,1)];
        end

        if contains(current_label, 'b0')
            metabolics_normalized_balance0 = [metabolics_normalized_balance0, data_labels.(current_label).metabolics.normalized(:,1)];
        end
        if contains(current_label, 'b1')
            metabolics_normalized_balance1 = [metabolics_normalized_balance1, data_labels.(current_label).metabolics.normalized(:,1)];
        end
        if contains(current_label, 'b2')
            metabolics_normalized_balance2 = [metabolics_normalized_balance2, data_labels.(current_label).metabolics.normalized(:,1)];
        end
    end

    metabolics_normalized_balance0_std = std(mean(metabolics_normalized_balance0,2));
    metabolics_normalized_balance1_std = std(mean(metabolics_normalized_balance1,2));
    metabolics_normalized_balance2_std = std(mean(metabolics_normalized_balance2,2));
    metabolics_normalized_speed0_std = std(mean(metabolics_normalized_speed0,2));
    metabolics_normalized_speed1_std = std(mean(metabolics_normalized_speed1,2));
    metabolics_normalized_speed2_std = std(mean(metabolics_normalized_speed2,2));
    metabolics_normalized_accuracy0_std = std(mean(metabolics_normalized_accuracy0,2));
    metabolics_normalized_accuracy1_std = std(mean(metabolics_normalized_accuracy1,2));
    metabolics_normalized_accuracy2_std = std(mean(metabolics_normalized_accuracy2,2));

    metabolics_normalized_balance0_avg = mean(metabolics_normalized_balance0(:));
    metabolics_normalized_balance1_avg = mean(metabolics_normalized_balance1(:));
    metabolics_normalized_balance2_avg = mean(metabolics_normalized_balance2(:));
    metabolics_normalized_speed0_avg = mean(metabolics_normalized_speed0(:));
    metabolics_normalized_speed1_avg = mean(metabolics_normalized_speed1(:));
    metabolics_normalized_speed2_avg = mean(metabolics_normalized_speed2(:));
    metabolics_normalized_accuracy0_avg = mean(metabolics_normalized_accuracy0(:));
    metabolics_normalized_accuracy1_avg = mean(metabolics_normalized_accuracy1(:));
    metabolics_normalized_accuracy2_avg = mean(metabolics_normalized_accuracy2(:));

    means = [metabolics_normalized_balance0_avg, metabolics_normalized_balance1_avg, metabolics_normalized_balance2_avg; 
             metabolics_normalized_speed0_avg, metabolics_normalized_speed1_avg, metabolics_normalized_speed2_avg; 
             metabolics_normalized_accuracy0_avg, metabolics_normalized_accuracy1_avg, metabolics_normalized_accuracy2_avg];

    stds = [metabolics_normalized_balance0_std, metabolics_normalized_balance1_std, metabolics_normalized_balance2_std; 
            metabolics_normalized_speed0_std, metabolics_normalized_speed1_std, metabolics_normalized_speed2_std; 
            metabolics_normalized_accuracy0_std, metabolics_normalized_accuracy1_std, metabolics_normalized_accuracy2_std];

    figure;
    hold on;
    categories = {'Visual Perturbation', 'Speed Prompt', 'Accuracy Prompt'};
    x = 1:3; % Groups
    bar_handles = bar(x, means, 'grouped'); 

    if plot_eachparticipant == 1
        metabolics_normalized_means = {mean(metabolics_normalized_balance0,2), mean(metabolics_normalized_balance1,2), mean(metabolics_normalized_balance2,2); 
                    mean(metabolics_normalized_speed0,2), mean(metabolics_normalized_speed1,2), mean(metabolics_normalized_speed2,2); 
                    mean(metabolics_normalized_accuracy0,2), mean(metabolics_normalized_accuracy1,2), mean(metabolics_normalized_accuracy2,2)};
        
        num_trials = size(metabolics_normalized_means{1,1}, 1); % Assuming all groups have the same number of trials
        %colors = lines(num_trials); % Use 'lines' colormap to generate distinct colors
        %colors = min(colors * 1.4, 1); % Increase intensity, limit to 1
        for i = 1:size(means, 2) % Loop through balance, speed, accuracy columns
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars
            
            % Scatter plot for individual points with unique colors per row
            for j = 1:3 % Loop through categories (Balance, Speed, Accuracy)
                data_points = metabolics_normalized_means{j, i}; % Extract row-wise mean values
                for k = 1:num_trials % Loop through trials (rows)
                    scatter(x_offset(j), data_points(k), point_size, colors_participants(k, :), 'filled', 'MarkerFaceAlpha', 0.6);
                end
            end
        end
    else
        for i = 1:size(means, 2)
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars with separate lower and upper limits
        end
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:3); 
    xticklabels(categories);
    ylabel('Mean Value (W/Kg)')
    legend({'Null/Low', 'Medium', 'High'});
    title('Normalized Energy Expenditure');
    grid on;
    hold off;

    %% Metabolics not normalized - Plot 68

    metabolics_not_normalized_balance0 = [];
    metabolics_not_normalized_balance1 = [];
    metabolics_not_normalized_balance2 = [];
    metabolics_not_normalized_speed0 = [];
    metabolics_not_normalized_speed1 = [];
    metabolics_not_normalized_speed2 = [];
    metabolics_not_normalized_accuracy0 = [];
    metabolics_not_normalized_accuracy1 = [];
    metabolics_not_normalized_accuracy2 = [];

    for i=1:length(list_labels)
        current_label = list_labels{i};

        if contains(current_label, 's0')
            metabolics_not_normalized_speed0 = [metabolics_not_normalized_speed0, data_labels.(current_label).metabolics.not_normalized(:,1)];
        end
        if contains(current_label, 's1')
            metabolics_not_normalized_speed1 = [metabolics_not_normalized_speed1, data_labels.(current_label).metabolics.not_normalized(:,1)];
        end
        if contains(current_label, 's2')
            metabolics_not_normalized_speed2 = [metabolics_not_normalized_speed2, data_labels.(current_label).metabolics.not_normalized(:,1)];
        end

        if contains(current_label, 'a0')
            metabolics_not_normalized_accuracy0 = [metabolics_not_normalized_accuracy0, data_labels.(current_label).metabolics.not_normalized(:,1)];
        end
        if contains(current_label, 'a1')
            metabolics_not_normalized_accuracy1 = [metabolics_not_normalized_accuracy1, data_labels.(current_label).metabolics.not_normalized(:,1)];
        end
        if contains(current_label, 'a2')
            metabolics_not_normalized_accuracy2 = [metabolics_not_normalized_accuracy2, data_labels.(current_label).metabolics.not_normalized(:,1)];
        end

        if contains(current_label, 'b0')
            metabolics_not_normalized_balance0 = [metabolics_not_normalized_balance0, data_labels.(current_label).metabolics.not_normalized(:,1)];
        end
        if contains(current_label, 'b1')
            metabolics_not_normalized_balance1 = [metabolics_not_normalized_balance1, data_labels.(current_label).metabolics.not_normalized(:,1)];
        end
        if contains(current_label, 'b2')
            metabolics_not_normalized_balance2 = [metabolics_not_normalized_balance2, data_labels.(current_label).metabolics.not_normalized(:,1)];
        end
    end

    metabolics_not_normalized_balance0_std = std(mean(metabolics_not_normalized_balance0,2));
    metabolics_not_normalized_balance1_std = std(mean(metabolics_not_normalized_balance1,2));
    metabolics_not_normalized_balance2_std = std(mean(metabolics_not_normalized_balance2,2));
    metabolics_not_normalized_speed0_std = std(mean(metabolics_not_normalized_speed0,2));
    metabolics_not_normalized_speed1_std = std(mean(metabolics_not_normalized_speed1,2));
    metabolics_not_normalized_speed2_std = std(mean(metabolics_not_normalized_speed2,2));
    metabolics_not_normalized_accuracy0_std = std(mean(metabolics_not_normalized_accuracy0,2));
    metabolics_not_normalized_accuracy1_std = std(mean(metabolics_not_normalized_accuracy1,2));
    metabolics_not_normalized_accuracy2_std = std(mean(metabolics_not_normalized_accuracy2,2));

    metabolics_not_normalized_balance0_avg = mean(metabolics_not_normalized_balance0(:));
    metabolics_not_normalized_balance1_avg = mean(metabolics_not_normalized_balance1(:));
    metabolics_not_normalized_balance2_avg = mean(metabolics_not_normalized_balance2(:));
    metabolics_not_normalized_speed0_avg = mean(metabolics_not_normalized_speed0(:));
    metabolics_not_normalized_speed1_avg = mean(metabolics_not_normalized_speed1(:));
    metabolics_not_normalized_speed2_avg = mean(metabolics_not_normalized_speed2(:));
    metabolics_not_normalized_accuracy0_avg = mean(metabolics_not_normalized_accuracy0(:));
    metabolics_not_normalized_accuracy1_avg = mean(metabolics_not_normalized_accuracy1(:));
    metabolics_not_normalized_accuracy2_avg = mean(metabolics_not_normalized_accuracy2(:));

    means = [metabolics_not_normalized_balance0_avg, metabolics_not_normalized_balance1_avg, metabolics_not_normalized_balance2_avg; 
             metabolics_not_normalized_speed0_avg, metabolics_not_normalized_speed1_avg, metabolics_not_normalized_speed2_avg; 
             metabolics_not_normalized_accuracy0_avg, metabolics_not_normalized_accuracy1_avg, metabolics_not_normalized_accuracy2_avg];

    stds = [metabolics_not_normalized_balance0_std, metabolics_not_normalized_balance1_std, metabolics_not_normalized_balance2_std; 
            metabolics_not_normalized_speed0_std, metabolics_not_normalized_speed1_std, metabolics_not_normalized_speed2_std; 
            metabolics_not_normalized_accuracy0_std, metabolics_not_normalized_accuracy1_std, metabolics_not_normalized_accuracy2_std];

    figure;
    hold on;
    categories = {'Visual Perturbation', 'Speed Prompt', 'Accuracy Prompt'};
    x = 1:3; % Groups
    bar_handles = bar(x, means, 'grouped'); 

    if plot_eachparticipant == 1
        metabolics_not_normalized_means = {mean(metabolics_not_normalized_balance0,2), mean(metabolics_not_normalized_balance1,2), mean(metabolics_not_normalized_balance2,2); 
                    mean(metabolics_not_normalized_speed0,2), mean(metabolics_not_normalized_speed1,2), mean(metabolics_not_normalized_speed2,2); 
                    mean(metabolics_not_normalized_accuracy0,2), mean(metabolics_not_normalized_accuracy1,2), mean(metabolics_not_normalized_accuracy2,2)};
        
        num_trials = size(metabolics_not_normalized_means{1,1}, 1); % Assuming all groups have the same number of trials
        %colors = lines(num_trials); % Use 'lines' colormap to generate distinct colors
        %colors = min(colors * 1.4, 1); % Increase intensity, limit to 1
        for i = 1:size(means, 2) % Loop through balance, speed, accuracy columns
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars
            
            % Scatter plot for individual points with unique colors per row
            for j = 1:3 % Loop through categories (Balance, Speed, Accuracy)
                data_points = metabolics_not_normalized_means{j, i}; % Extract row-wise mean values
                for k = 1:num_trials % Loop through trials (rows)
                    scatter(x_offset(j), data_points(k), point_size, colors_participants(k, :), 'filled', 'MarkerFaceAlpha', 0.6);
                end
            end
        end
    else
        for i = 1:size(means, 2)
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars with separate lower and upper limits
        end
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:3); 
    xticklabels(categories);
    ylabel('Mean Value (W)')
    legend({'Null/Low', 'Medium', 'High'});
    title('Energy Expenditure');
    grid on;
    hold off;

    %% Survey metabolics by condition NOT NORMALIZED - Plot 69
    
    survey_metabolics_balance0 = [];
    survey_metabolics_balance1 = [];
    survey_metabolics_balance2 = [];
    survey_metabolics_speed0 = [];
    survey_metabolics_speed1 = [];
    survey_metabolics_speed2 = [];
    survey_metabolics_accuracy0 = [];
    survey_metabolics_accuracy1 = [];
    survey_metabolics_accuracy2 = [];

    for i=1:length(list_labels)
        current_label = list_labels{i};

        if contains(current_label, 's0')
            survey_metabolics_speed0 = [survey_metabolics_speed0, data_labels.(current_label).survey_metabolics(:,1)];
        end
        if contains(current_label, 's1')
            survey_metabolics_speed1 = [survey_metabolics_speed1, data_labels.(current_label).survey_metabolics(:,1)];
        end
        if contains(current_label, 's2')
            survey_metabolics_speed2 = [survey_metabolics_speed2, data_labels.(current_label).survey_metabolics(:,1)];
        end

        if contains(current_label, 'a0')
            survey_metabolics_accuracy0 = [survey_metabolics_accuracy0, data_labels.(current_label).survey_metabolics(:,1)];
        end
        if contains(current_label, 'a1')
            survey_metabolics_accuracy1 = [survey_metabolics_accuracy1, data_labels.(current_label).survey_metabolics(:,1)];
        end
        if contains(current_label, 'a2')
            survey_metabolics_accuracy2 = [survey_metabolics_accuracy2, data_labels.(current_label).survey_metabolics(:,1)];
        end

        if contains(current_label, 'b0')
            survey_metabolics_balance0 = [survey_metabolics_balance0, data_labels.(current_label).survey_metabolics(:,1)];
        end
        if contains(current_label, 'b1')
            survey_metabolics_balance1 = [survey_metabolics_balance1, data_labels.(current_label).survey_metabolics(:,1)];
        end
        if contains(current_label, 'b2')
            survey_metabolics_balance2 = [survey_metabolics_balance2, data_labels.(current_label).survey_metabolics(:,1)];
        end
    end

    survey_metabolics_balance0_std = std(mean(survey_metabolics_balance0,2));
    survey_metabolics_balance1_std = std(mean(survey_metabolics_balance1,2));
    survey_metabolics_balance2_std = std(mean(survey_metabolics_balance2,2));
    survey_metabolics_speed0_std = std(mean(survey_metabolics_speed0,2));
    survey_metabolics_speed1_std = std(mean(survey_metabolics_speed1,2));
    survey_metabolics_speed2_std = std(mean(survey_metabolics_speed2,2));
    survey_metabolics_accuracy0_std = std(mean(survey_metabolics_accuracy0,2));
    survey_metabolics_accuracy1_std = std(mean(survey_metabolics_accuracy1,2));
    survey_metabolics_accuracy2_std = std(mean(survey_metabolics_accuracy2,2));

    survey_metabolics_balance0_avg = mean(survey_metabolics_balance0(:));
    survey_metabolics_balance1_avg = mean(survey_metabolics_balance1(:));
    survey_metabolics_balance2_avg = mean(survey_metabolics_balance2(:));
    survey_metabolics_speed0_avg = mean(survey_metabolics_speed0(:));
    survey_metabolics_speed1_avg = mean(survey_metabolics_speed1(:));
    survey_metabolics_speed2_avg = mean(survey_metabolics_speed2(:));
    survey_metabolics_accuracy0_avg = mean(survey_metabolics_accuracy0(:));
    survey_metabolics_accuracy1_avg = mean(survey_metabolics_accuracy1(:));
    survey_metabolics_accuracy2_avg = mean(survey_metabolics_accuracy2(:));

    means = [survey_metabolics_balance0_avg, survey_metabolics_balance1_avg, survey_metabolics_balance2_avg; 
             survey_metabolics_speed0_avg, survey_metabolics_speed1_avg, survey_metabolics_speed2_avg; 
             survey_metabolics_accuracy0_avg, survey_metabolics_accuracy1_avg, survey_metabolics_accuracy2_avg];
    
    stds = [survey_metabolics_balance0_std, survey_metabolics_balance1_std, survey_metabolics_balance2_std; 
            survey_metabolics_speed0_std, survey_metabolics_speed1_std, survey_metabolics_speed2_std; 
            survey_metabolics_accuracy0_std, survey_metabolics_accuracy1_std, survey_metabolics_accuracy2_std];
        
    figure;
    hold on;
    categories = {'Visual Perturbation', 'Speed Prompt', 'Accuracy Prompt'};
    x = 1:3; % Groups
    bar_handles = bar(x, means, 'grouped'); 

    if plot_eachparticipant == 1
        points_means = {mean(survey_metabolics_balance0,2), mean(survey_metabolics_balance1,2), mean(survey_metabolics_balance2,2); 
                    mean(survey_metabolics_speed0,2), mean(survey_metabolics_speed1,2), mean(survey_metabolics_speed2,2); 
                    mean(survey_metabolics_accuracy0,2), mean(survey_metabolics_accuracy1,2), mean(survey_metabolics_accuracy2,2)};
        
        num_trials = size(points_means{1,1}, 1); % Assuming all groups have the same number of trials
        %colors = lines(num_trials); % Use 'lines' colormap to generate distinct colors
        %colors = min(colors * 1.4, 1); % Increase intensity, limit to 1
        for i = 1:size(means, 2) % Loop through balance, speed, accuracy columns
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            lower_error = min(means(:,i), stds(:,i)); % Ensure no negative error
            upper_error = stds(:,i); % Upper limit stays as std
            errorbar(x_offset, means(:,i), lower_error, upper_error, 'k', 'linestyle', 'none', 'LineWidth', 1.5); % Plot error bars
            
            % Scatter plot for individual points with unique colors per row
            for j = 1:3 % Loop through categories (Balance, Speed, Accuracy)
                data_points = points_means{j, i}; % Extract row-wise mean values
                for k = 1:num_trials % Loop through trials (rows)
                    scatter(x_offset(j), data_points(k), point_size, colors_participants(k, :), 'filled', 'MarkerFaceAlpha', 0.6);
                end
            end
        end
    else
        for i = 1:size(means, 2) % Loop through each group
            % Error bars
            x_offset = bar_handles(i).XEndPoints; % Get bar positions
            errorbar(x_offset, means(:,i), stds(:,i), 'k', 'linestyle', 'none', 'LineWidth', 1.5);
        end
    end
    
    set(gca, 'FontSize', 22);
    xticks(1:3); 
    xticklabels(categories);
    ylabel('Normalized Priority');
    legend({'Null/Low', 'Medium', 'High'});
    title('Perceived Energy Expenditure Reduction Priority');
    grid on;
    hold off;

end