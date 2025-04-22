 % Group the metabolic data by accuracy, speed, and balance. Metabolic value for each trial (doesn't include 4 short ones)
 function[] = Metabolics(imported_data,error,colors,metabolic_data,baseline_errors,normalize,bmh)

    for i=1:size(fieldnames(imported_data),1) %through bmhs
        names = string(fieldnames(imported_data));
        bmh = names(i);
        weights.(bmh) = imported_data.(bmh).participant_data.weight;
    end

    if nargin < 7
        bmh_text = "Across all Participants";
        participants = "multiple";
    else
        participants = "single";
        bmh_text = string(bmh);
    end

    if normalize == 1
        normalize_text = "Normalized to Weight";
    else
        normalize_text = "Not Normalized";
    end

    num_bars = 3;
    
    accuracy_levels = ["Low", "Medium", "High", "Unknown"];
    speed_levels = ["Slow", "Medium", "Fast", "Unknown"]; % Speed categories
    balance_levels = ["None", "Medium", "High", "Unknown"]; % Balance categories
    
    %% Classify trials according to their labels

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

    %metabolics = [];
    walking_speed = [];
    step_width = [];
    step_length = [];
    step_error = [];

    if participants == "single"
        metabolics = metabolic_data.(bmh);
        walking_speed = imported_data.(bmh).walkingspeed;
        step_width = imported_data.(bmh).meanwidthstraights;
        step_width_var = imported_data.(bmh).straightsvariability;  
        step_length = imported_data.(bmh).meanlength_straights;
        step_error = imported_data.(bmh).overall_errors.straight_mean_absolute;
    elseif participants == "multiple"
        for i=1:size(fieldnames(imported_data),1)
            names = string(fieldnames(imported_data));
            bmh = names(i);
            metabolics(i,:) = table2array(metabolic_data.(bmh));
            walking_speed(i,:) = imported_data.(bmh).walkingspeed;
            step_width(i,:) = imported_data.(bmh).meanwidthstraights;
            step_width_var(i,:) = imported_data.(bmh).straightsvariability;
            step_length(i,:) = imported_data.(bmh).meanlength_straights;
            step_error(i,:) = imported_data.(bmh).overall_errors.straight_mean_absolute;
        end
   end

    list_labels = fieldnames(data_labels); %labels stored in the struct

    bmh_list = string(fieldnames(baseline_errors));
    if size(bmh_list) > size(fieldnames(weights))
        error("The weight of a participant is missing. Introduce it in Metabolics.m")
    end

    for i=1:size(fieldnames(imported_data),1) %through bmhs
        names = string(fieldnames(imported_data));
        bmh = names(i);
        numtrials = size(imported_data.(bmh).labels,2);
        labels = imported_data.(bmh).labels;  % Format: 'sxaybz'

        for j = 3:29 %through trials
            current_label = labels{j}; % from the ordered label array
            trialname = sprintf('Trial%04d', j); 
            for k=1:size(list_labels,1)
                label_struct = list_labels{k};
                if strcmp(current_label, label_struct)
                    data_labels.(current_label).walking_speed(i,:) = walking_speed(i,j);
                    data_labels.(current_label).metabolics(i,:) = metabolics(i,j);
                    data_labels.(current_label).step_width(i,:) = step_width(i,j);
                    data_labels.(current_label).step_width_var(i,:) = step_width_var(i,j);
                    data_labels.(current_label).step_length(i,:) = step_length(i,j);
                    data_labels.(current_label).step_error(i,:) = step_error(i,j);
                else
                    continue
                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ACCURACY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure();
    hold on;
    boxnum = 1;
    labels = {};
    bar_positions = 1:num_bars; % Only plotting the average, no need for separate conditions
    
    for i = 1:num_bars
        accuracy = accuracy_levels(i); 
        fns.accuracy.(accuracy) = fieldnames(error.byaccuracy.(accuracy));
        accuracy_data.(accuracy) = [];
    
        for j = 1:length(fns.accuracy.(accuracy)) % Loop over all trials from every accuracy level
            trialname = fns.accuracy.(accuracy){j};
            trial_num = str2double(extractBetween(trialname, 7, 9));
            if participants == "multiple"
                bmh = string(extractBetween(trialname, 11, 15));
            end
            if trial_num > 31 || trial_num < 3 % Skip trials that don't have metabolic data
                continue
            end
            if participants == "single"
                data = metabolic_data.(bmh)(1, trial_num);
                title({'Energy Expenditure Grouped by Accuracy - ' + normalize_text + ' - ' + bmh_text}, 'FontSize', 20);
            elseif participants == "multiple"
                data = table2array(metabolic_data.(bmh)(1, trial_num));
                title({'Energy Expenditure Grouped by Accuracy ' + normalize_text + ' - ' + bmh_text}, 'FontSize', 20);
            end

            if normalize == 1
                data = data / weights.(bmh);
            end
            accuracy_data.(accuracy) = [accuracy_data.(accuracy), data];
        end
    
        % Calculate the mean and std for each accuracy level and plot a single bar
        mean_data = mean(accuracy_data.(accuracy));
        std_data = std(accuracy_data.(accuracy));
        bar(bar_positions(i), mean_data, 'FaceColor', colors(i, :), 'EdgeColor', 'black');
        errorbar(bar_positions(i), mean_data, std_data, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
        labels{end+1} = accuracy;
    end
    
    set(gca, 'XTick', bar_positions, 'XTickLabel', labels, 'FontSize', 14);

    if normalize == 1
        ylabel('Normalized Energy Expenditure (W/Kg)', 'FontSize', 14);
    else
        ylabel('Energy Expenditure (W)', 'FontSize', 14);
    end
    
    xlabel('Accuracy Levels', 'FontSize', 14);
    hold off;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SPEED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure();
    hold on;
    boxnum = 1;
    labels = {};
    bar_positions = 1:num_bars;
    
    for i = 1:num_bars
        speed = speed_levels(i); 
        fns.speed.(speed) = fieldnames(error.byspeed.(speed));
        speed_data.(speed) = [];
    
        for j = 1:length(fns.speed.(speed)) 
            trialname = fns.speed.(speed){j};
            trial_num = str2double(extractBetween(trialname, 7, 9));
            if participants == "multiple"
                bmh = string(extractBetween(trialname, 11, 15));
            end
            if trial_num > 31 || trial_num < 3  % Skip trials that don't have metabolic data
                continue
            end
    
            if participants == "single"
                data = metabolic_data.(bmh)(1, trial_num);
                title({'Energy Expenditure Grouped by Speed - ' + normalize_text + ' - ' + bmh_text}, 'FontSize', 14);
            elseif participants == "multiple"
                data = table2array(metabolic_data.(bmh)(1, trial_num));
                title({'Energy Expenditure Grouped by Speed - ' + normalize_text + ' - ' + bmh_text}, 'FontSize', 14);
            end

            if normalize == 1
                data = data / weights.(bmh);
            end

            speed_data.(speed) = [speed_data.(speed), data];
        end
    
        % Calculate the mean and std for each speed level and plot a single bar
        mean_data = mean(speed_data.(speed));
        std_data = std(speed_data.(speed));
        bar(bar_positions(i), mean_data, 'FaceColor', colors(i, :), 'EdgeColor', 'black');
        errorbar(bar_positions(i), mean_data, std_data, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
        labels{end+1} = speed;
    end
    
    set(gca, 'XTick', 1:num_bars, 'XTickLabel', labels, 'FontSize', 15);

    if normalize == 1
        ylabel('Normalized Energy Expenditure (W/Kg)', 'FontSize', 15);
    else
        ylabel('Energy Expenditure (W)', 'FontSize', 15);
    end

    xlabel('Speed Levels', 'FontSize', 15);

    hold off;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BALANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    figure();
    hold on;
    boxnum = 1;
    labels = {};
    bar_positions = 1:num_bars;
    
    for i = 1:num_bars
        balance = balance_levels(i); 
        fns.balance.(balance) = fieldnames(error.bybalance.(balance));
        balance_data.(balance) = [];
    
        for j = 1:length(fns.balance.(balance)) 
            trialname = fns.balance.(balance){j};
            trial_num = str2double(extractBetween(trialname, 7, 9));
            if participants == "multiple"
                bmh = string(extractBetween(trialname, 11, 15));
            end
            if trial_num > 31 || trial_num < 3  % Skip trials that don't have metabolic data
                continue
            end
    
            if participants == "single"
                data = metabolic_data.(bmh)(1, trial_num);
                title({'Energy Expenditure Grouped by Balance - ' + normalize_text + ' - ' + bmh_text}, 'FontSize', 14);
            elseif participants == "multiple"
                data = table2array(metabolic_data.(bmh)(1, trial_num));
                title({'Energy Expenditure Grouped by Balance - ' + normalize_text + ' - ' + bmh_text}, 'FontSize', 14);
            end

            if normalize == 1
                data = data / weights.(bmh);
            end

            balance_data.(balance) = [balance_data.(balance), data];
        end
    
        % Calculate the mean and std for each balance level and plot a single bar
        mean_data = mean(balance_data.(balance));
        std_data = std(balance_data.(balance));
        bar(bar_positions(i), mean_data, 'FaceColor', colors(i, :), 'EdgeColor', 'black');
        errorbar(bar_positions(i), mean_data, std_data, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
        labels{end+1} = balance;
    end
    
    set(gca, 'XTick', 1:num_bars, 'XTickLabel', labels, 'FontSize', 15);
    
    if normalize == 1
        ylabel('Normalized Energy Expenditure (W/Kg)', 'FontSize', 15);
    else
        ylabel('Energy Expenditure (W)', 'FontSize', 15);
    end

    xlabel('Balance Levels', 'FontSize', 15);
    hold off;


    %% Metabolics for each trial (label); Bar plot for each label
    
    means = zeros(size(list_labels, 1), 1); 
    stds = zeros(size(list_labels, 1), 1); 

    for i = 1:size(list_labels, 1)
        label = string(list_labels{i}); 
        if normalize == 0
            means(i) = mean(data_labels.(label).metabolics); 
            stds(i) = std(data_labels.(label).metabolics);
        elseif normalize == 1
            for j=1:size(data_labels.(label).metabolics)
                norm(j) = mean(data_labels.(label).metabolics)/weights.(bmh_list(j)); 
            end
            means(i) = mean(norm);
            stds(i) = std(norm);
        end
    end
    
    figure;
    bar(means); 
    hold on;
    errorbar(1:size(list_labels, 1), means, stds, 'k.', 'LineWidth', 1.5); 
    xlabel('Labels'); 
    if normalize == 0
        ylabel('Energy Expenditure (W)');
    elseif normalize == 1
        ylabel('Normalized Energy Expenditure (W/Kg)');
    end
    
    xticks(1:size(list_labels, 1)); 
    xticklabels(list_labels); 
    title('Mean Metabolics for Each Label - ' + normalize_text + ' - ' + bmh_text);
    set(gca, 'FontSize', 16);
    
    hold off;

    %% By condition (PLOT 5)

    balance_avg = zeros(1, 3);
    balance_std = zeros(1, 3);
    speed_avg = zeros(1, 3);
    speed_std = zeros(1, 3);
    accuracy_avg = zeros(1, 3);
    accuracy_std = zeros(1, 3);
    
    balance_avg(1) = mean(balance_data.None);
    balance_std(1) = std(balance_data.None);
    
    balance_avg(2) = mean(balance_data.Medium);
    balance_std(2) = std(balance_data.Medium);
    
    balance_avg(3) = mean(balance_data.High);
    balance_std(3) = std(balance_data.High);
    
    speed_avg(1) = mean(speed_data.Slow);
    speed_std(1) = std(speed_data.Slow);
    
    speed_avg(2) = mean(speed_data.Medium);
    speed_std(2) = std(speed_data.Medium);
    
    speed_avg(3) = mean(speed_data.Fast);
    speed_std(3) = std(speed_data.Fast);
    
    accuracy_avg(1) = mean(accuracy_data.Low);
    accuracy_std(1) = std(accuracy_data.Low);
    
    accuracy_avg(2) = mean(accuracy_data.Medium);
    accuracy_std(2) = std(accuracy_data.Medium);
    
    accuracy_avg(3) = mean(accuracy_data.High);
    accuracy_std(3) = std(accuracy_data.High);
    
    avg_data = [balance_avg; speed_avg; accuracy_avg]; % Rows: Conditions, Columns: Levels
    std_data = [balance_std; speed_std; accuracy_std];
    
    figure;
    bar_handle = bar(avg_data); % Create grouped bars
    hold on;
    
    % Add error bars
    num_groups = size(avg_data, 1); % Number of conditions (balance, speed, accuracy)
    num_bars = size(avg_data, 2); % Number of levels (none/low, medium, high)
    group_width = min(0.8, num_bars/(num_bars + 1.5)); % Adjust group width
    
    for i = 1:num_bars
        % X positions of bars within each group
        x = (1:num_groups) - group_width/2 + (2*i-1) * group_width / (2*num_bars);
        errorbar(x, avg_data(:,i), std_data(:,i), 'k.', 'LineWidth', 1.2); % Add error bars
    end
    
    xticks(1:num_groups);
    xticklabels({'Balance Perturbation', 'Speed', 'Accuracy'}); % X-axis labels
    %xlabel('Condition');

    if normalize == 0
        ylabel('Energy Expenditure (W)');
        title('Energy Expenditure by Condition');
    elseif normalize == 1
        ylabel('Mean Value (W/kg)');
        title('Normalized Energy Expenditure');
    end
    
    set(gca, 'FontSize', 16);
    legend({'None/Low', 'Medium', 'High/Fast'});
    grid on;
    hold off;

    %% Metabolics comparison (PLOT 6)

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

    speed_balance0_std = std(speed_balance0(:));
    speed_balance1_std = std(speed_balance1(:));
    speed_balance2_std = std(speed_balance2(:));
    speed_speed0_std = std(speed_speed0(:));
    speed_speed1_std = std(speed_speed1(:));
    speed_speed2_std = std(speed_speed2(:));
    speed_accuracy0_std = std(speed_accuracy0(:));
    speed_accuracy1_std = std(speed_accuracy1(:));
    speed_accuracy2_std = std(speed_accuracy2(:));

    speed_balance0 = mean(speed_balance0(:));
    speed_balance1 = mean(speed_balance1(:));
    speed_balance2 = mean(speed_balance2(:));
    speed_speed0 = mean(speed_speed0(:));
    speed_speed1 = mean(speed_speed1(:));
    speed_speed2 = mean(speed_speed2(:));
    speed_accuracy0 = mean(speed_accuracy0(:));
    speed_accuracy1 = mean(speed_accuracy1(:));
    speed_accuracy2 = mean(speed_accuracy2(:));

    % Define the speed range (0.6 to 1.4 m/s)
    S = linspace(0.6, 1.4, 100);
    
    % Calculate the energy expenditure (EE) for the given range
    EE = 1.44 * 1.94 * S.^0.43 + 0.24 * S.^4;
    
    % Example means and standard deviations for speed condition (replace with your actual data)
    speed_means = [speed_speed0, speed_speed1, speed_speed2]; % Speeds for the speed condition
    metabolics_byspeed_avg = speed_avg; % Metabolics for the speed condition
    metabolics_byspeed_std = speed_std; % Standard deviation of metabolics
    
    figure;
    plot(S, EE, 'LineWidth', 2, 'DisplayName', 'Expected Values (Literature)');
    hold on;
    
    % Plot the dashed line for speed
    plot(speed_means, metabolics_byspeed_avg, '--', 'LineWidth', 1.5, 'DisplayName', 'Experimental Values', 'Color', 'r'); % Red line for speed
    
    % Plot the individual points
    scatter(speed_means, metabolics_byspeed_avg, 50, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off'); % Exclude points from legend
    
    % Add shaded area for standard deviation
    fill([speed_means, fliplr(speed_means)], [metabolics_byspeed_avg + metabolics_byspeed_std, fliplr(speed_avg - metabolics_byspeed_std)], ...
         'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off'); % Exclude shaded area from legend
    
    set(gca, 'FontSize', 16);
    xlabel('Walking Speed (m/s)');
    ylabel('Energy Expenditure (W/kg)');
    title('Energy Expenditure with Speed Condition');
    legend('show', 'Location', 'Best'); % Show only the desired legend entries
    xlim([0.82 1.25])
    grid on;
    hold off;
    
    s1a0b0_speed_avg = mean(data_labels.s1a0b0.walking_speed);
    s1a0b0_speed_std = std(data_labels.s1a0b0.walking_speed);

    s1a0b0_metabolics_avg = mean(data_labels.s1a0b0.metabolics);
    
    if normalize == 1
        s1a0b0_metabolics_avg = s1a0b0_metabolics_avg / weights.(bmh);
    end

    s1a0b0_metabolics_std = std(data_labels.s1a0b0.metabolics / weights.(bmh));

    EE_s1a0b0 = 1.44 * 1.94 * s1a0b0_speed_avg.^0.43 + 0.24 * s1a0b0_speed_avg.^4;

    disp(['Walking speed at s1a0b0: ', num2str(s1a0b0_speed_avg), ' +- ', num2str(s1a0b0_speed_std), ...
      '. Recorded metabolic value: ', num2str(s1a0b0_metabolics_avg), ' +- ', num2str(s1a0b0_metabolics_std)]);

    disp(['Metabolic value for that speed according to the literature: ', num2str(EE_s1a0b0)]);

 end % End of the function


