function[] = BoxPlot_error(imported_data,error,n,section,foot,colors,subtract_speeds,subtract_all,bmh)

    num_bars = 3; %3 to NOT plot 3s (unknown bar; short trials and calibration), 4 to plot EVERYTHING
    raincloud = 0;
    
    % Check if bmh was provided; if not, set it to 0
    if nargin < 9
        bmh_text = "Across all Participants";
        participants = "multiple";
    else
        participants = "single";
        bmh_text = string(bmh);
    end

    if subtract_speeds == 1 && subtract_all == 1
        error("subtract_speeds and subtract_all cannot be both 1 at the same time")
    end

    if n==3
        n_text = "Longitudinal Error (Absolute Value)";
    elseif n==4
        n_text = "Longitudinal Error";
    elseif n==5
        n_text = "Total Error";
    elseif n==6
        n_text = "Lateral Error";
    end

    if section==1
        section_text = "All the Circuit";
        m_path = "overall_mean_signed"; %Only used when plotting the signed error
    elseif section==2
        section_text = "Straight Sections";
        m_path = "straight_mean_signed"; %Only used when plotting the signed error
    elseif section==3
        section_text = "Curved Sections";
        m_path = "curve_mean_signed"; %Only used when plotting the signed error
    end

    if foot==1
        foot_text = "Both Feet";
    elseif foot==2
        foot_text = "Right Foot";
    elseif foot==3
        foot_text = "Left Foot";
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ACCURACY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure();
    hold on;
    boxnum = 1;
    labels = {}; % Initialize a cell array to store labels
    accuracy_levels = ["Low", "Medium", "High", "Unknown"]; % Accuracy categories

    for i = 1:num_bars % loop through accuracy
        accuracy = accuracy_levels(i); % Use pre-defined accuracy levels
        fns.accuracy.(accuracy) = fieldnames(error.byaccuracy.(accuracy));
        
        for j = 1:length(fns.accuracy.(accuracy)) % loop through trials for each accuracy
            trialname = fns.accuracy.(accuracy){j};
            data = error.byaccuracy.(accuracy).(trialname); % get data
            boxplot(data / 10, 'Positions', boxnum, 'Widths', 0.2, 'Colors', colors(i, :), 'Symbol', '_k'); % it plots here
            labels{boxnum} = strrep(trialname, 'Trial00', ''); % Store the label for the current box plot
            boxnum = boxnum + 1;
        end
    end
    
    % Create legend manually
    p = arrayfun(@(x) plot(NaN, NaN, 'Color', colors(x, :)), 1:num_bars); 
    legend(p, accuracy_levels(1:num_bars));
    set(gca, 'XTick', 1:boxnum-1, 'XTickLabel', labels, 'FontSize', 15);
    ylabel('Step Error (cm)', 'FontSize', 15);
    xlabel('Trial Number', 'FontSize', 15);
    
    title({'Step Error Grouped by Accuracy Prompt - ' + bmh_text + ' - ' + n_text; section_text + ' - ' + foot_text}, 'FontSize', 20);
    hold off;
    
    if n == 6
        figure;
        hold on;
        
        avg_all = []; % Initialize for overall averages
        labels = {}; % To store trial labels
        boxnum = 1; % To track x-axis position
        
        % Loop through accuracies
        for i = 1:num_bars
            accuracy = accuracy_levels(i);
            
            for j = 1:length(fns.accuracy.(accuracy))
                trialname = fns.accuracy.(accuracy){j};
                data = error.byaccuracy.(accuracy).(trialname); % Get data for the current trial
                avg_all(end+1) = mean(data / 10); % Convert from section to cm
                
                labels{boxnum} = strrep(trialname, 'Trial00', ''); % Store label
                boxnum = boxnum + 1;
            end
        end
        
        % Plot the bar graph for the overall averages
        bar_handle = bar(avg_all, 'FaceColor', 'flat');  % Set 'FaceColor' to 'flat'
        
        % Assign colors based on accuracy
        for i = 1:num_bars
            accuracy = accuracy_levels(i);
            for j = 1:length(fns.accuracy.(accuracy))
                trial_idx = find(strcmp(labels, strrep(fns.accuracy.(accuracy){j}, 'Trial00', '')));
                bar_handle.CData(trial_idx, :) = colors(i, :); % Assign color
            end
        end
        
        % Customize plot
        set(gca, 'XTick', 1:length(labels), 'XTickLabel', labels, 'FontSize', 15);
        ylabel('Average Step Error (cm)', 'FontSize', 15);
        xlabel('Trial Number', 'FontSize', 15);
        
        % Create legend manually
        p = arrayfun(@(x) plot(NaN, NaN, 'Color', colors(x, :)), 1:num_bars);
        legend(p, accuracy_levels(1:num_bars));
        title('Average Lateral Step Error Grouped by Accuracy Condition - ' + bmh_text, 'FontSize', 20);
        hold off;
    end

    % BAR PLOT
    
    % Pre-allocate arrays for means and standard deviations
    means = zeros(1, 3); % Now only 3 categories (Low, Medium, High)
    stds = zeros(1, 3);
    all_data = cell(1, 3); % Store data for each accuracy level
    
    % Loop through the accuracies and compute mean, std, and collect data for each (without Unknown)
    for i = 1:num_bars
        accuracy = accuracy_levels(i);
        accuracy_data.(accuracy) = []; % Initialize array to collect all data for the current accuracy level
    
        for j = 1:length(fns.accuracy.(accuracy)) % Loop through trials for each accuracy
            trialname = fns.accuracy.(accuracy){j};
            if participants == "multiple"
                trial_num = str2double(extractBetween(trialname, 7, 9));  % Extract trial number
            else
                trial_num = str2double(strrep(trialname, 'Trial00', '')); % Extract trial number
            end
    
            if trial_num >= 3 && trial_num <= 29
                data = error.byaccuracy.(accuracy).(trialname) / 10; % Convert from mm to cm
                accuracy_data.(accuracy) = [accuracy_data.(accuracy); data]; % Append current trial data
            end
        end
        
        % Store data for each accuracy level for plotting
        all_data{i} = accuracy_data.(accuracy);
        means(i) = mean(accuracy_data.(accuracy));
        stds(i) = std(accuracy_data.(accuracy));
    end
    
    % Create the bar graph with error bars
    figure;
    b = bar(means); % Bar graph of means
    b.FaceColor = 'flat'; % Enable color customization
    for k = 1:num_bars
        b.CData(k, :) = colors(k, :); % Assign a different color to each bar
    end
    hold on;
    errorbar(1:num_bars, means, stds, 'k', 'linestyle', 'none'); % Add error bars
    
    % Customize plot
    set(gca, 'XTick', 1:num_bars, 'XTickLabel', accuracy_levels(1:num_bars), 'FontSize', 15);
    ylabel('Step Error (cm)', 'FontSize', 15);
    xlabel('Accuracy Level', 'FontSize', 15);
    title({'Step Error Grouped by Accuracy - ' + bmh_text + ' - ' + n_text; section_text + ' - ' + foot_text}, 'FontSize', 20);
    
    % RAINCLOUD PLOT (Refined with reduction factor and closer bars)
    if raincloud == 1
        reduction_factor = 4; % Change this value to reduce points by 1/x (e.g., 2 for half, 3 for a third, etc.)
        figure;
        hold on;
        
        % Define offset to bring bars closer
        x_offset = 0.2; % Adjust to make bars closer, smaller values bring them closer
        
        for i = 1:num_bars
            % Generate a narrower "Gaussian" effect with kernel density estimation
            [f, xi] = ksdensity(all_data{i}, 'Bandwidth', 0.2); % Set bandwidth for narrower distribution
            f = f / max(f) * 0.1; % Scale to a smaller size (height and width)
        
            % Plot Gaussian "cloud" on one side of the bar
            fill(f + i * x_offset, xi, colors(i, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        
            % Overlay jittered individual data points, reduced by the defined factor
            num_points = round(length(all_data{i}) / reduction_factor); % Calculate reduced number of points
            sample_data = all_data{i}(randperm(length(all_data{i}), num_points)); % Select subset of points
            scatter(i * x_offset * ones(size(sample_data)) + 0.05 * (rand(size(sample_data)) - 0.5), ...
                sample_data, 15, colors(i, :), 'filled', 'MarkerFaceAlpha', 0.6);
        
            % Plot boxplot for central tendency and spread
            boxplot(all_data{i}, 'positions', i * x_offset, 'Widths', 0.05, 'Colors', 'k', 'Symbol', '');
        end
        
        % Adjust plot spacing for closer bars
        set(gca, 'XTick', x_offset * (1:num_bars), 'XTickLabel', accuracy_levels(1:num_bars), 'FontSize', 15);
        ylabel('Step Error (cm)', 'FontSize', 15);
        xlabel('Accuracy Level', 'FontSize', 15);
        title({'Step Error Grouped by Accuracy - ' + bmh_text + ' - ' + n_text; section_text + ' - ' + foot_text}, 'FontSize', 20);
        
        xlim([0.5 * x_offset, (num_bars + 0.5) * x_offset]); % Bring bars even closer by narrowing x-axis limits
        hold off;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SPEEDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure();
    hold on;
    boxnum = 1;
    labels = {}; % Initialize a cell array to store labels
    speed_levels = ["Slow", "Medium", "Fast", "Unknown"]; % Speed categories
    subtract_text = "";
    
    for i = 1:num_bars % Loop through speeds
        speed = speed_levels(i); % Use pre-defined speed levels
        fns.speed.(speed) = fieldnames(error.byspeed.(speed));
    
        for j = 1:length(fns.speed.(speed)) % Loop through trials for each speed
            trialname = fns.speed.(speed){j};
            data = error.byspeed.(speed).(trialname); % Get data
            if participants == "multiple"
                bmh = regexp(trialname, 'BMH\d{2}', 'match', 'once');
            end

            [speed_s0a2b0, speed_s1a2b0, speed_s2a2b0] = BaselineErrors(imported_data,bmh,n,section); % Extract the baseline for this participant

            if subtract_all == 1 && (n == 3 || n == 4)
                warning('Both feet used; not implemented yet for just right or left foot');
                val = (speed_s0a2b0 + speed_s1a2b0 + speed_s2a2b0) / 3;
                boxplot((data + val) / 10, 'Positions', boxnum, 'Widths', 0.2, 'Colors', colors(i, :), 'Symbol', '_k'); % it plots here
                subtract_text = " - Mean Subtracted";              
            elseif subtract_speeds == 1 && (n == 3 || n == 4)
                warning('Both feet used; not implemented yet for just right or left foot');
                if speed == "Slow"
                    val = speed_s0a2b0;
                    means(1) = means(1) + val; 
                elseif speed == "Medium"
                    val = speed_s1a2b0;
                    means(2) = means(2) + val; 
                elseif speed == "Fast"
                    val = speed_s2a2b0;
                    means(3) = means(3) + val;   
                end
                boxplot((data + val) / 10, 'Positions', boxnum, 'Widths', 0.2, 'Colors', colors(i, :), 'Symbol', '_k'); % it plots here
                subtract_text = " - Mean Subtracted According to Speed";   
            else %both variables are 0
                boxplot(data / 10, 'Positions', boxnum, 'Widths', 0.2, 'Colors', colors(i, :), 'Symbol', '_k'); % it plots here
            end

            labels{boxnum} = strrep(trialname, 'Trial00', ''); % Store the label for the current box plot
            boxnum = boxnum + 1;
        end
    end
    
    % Create legend manually
    p = arrayfun(@(x) plot(NaN, NaN, 'Color', colors(x, :)), 1:num_bars);
    legend(p, speed_levels);
    set(gca, 'XTick', 1:boxnum-1, 'XTickLabel', labels, 'FontSize', 15);
    ylabel('Step Error (cm)', 'FontSize', 15);
    xlabel('Trial Number', 'FontSize', 15);
    
    title({'Step Error Grouped by Speed Prompt - ' + bmh_text + ' - ' + n_text; section_text + ' - ' + foot_text + subtract_text}, 'FontSize', 20);
    hold off;
    
    if n == 6
        figure;
        hold on;
    
        avg_all = []; % Initialize for overall averages
        labels = {}; % To store trial labels
        boxnum = 1; % To track x-axis position
    
        % Loop through speeds
        for i = 1:num_bars
            speed = speed_levels(i); % Use pre-defined speed levels
    
            for j = 1:length(fns.speed.(speed))
                trialname = fns.speed.(speed){j};
                data = error.byspeed.(speed).(trialname); % Get data for the current trial
                avg_all(end+1) = mean(data / 10); % Convert from section to cm
    
                labels{boxnum} = strrep(trialname, 'Trial00', ''); % Store label
                boxnum = boxnum + 1;
            end
        end
    
        % Plot the bar graph for the overall averages
        bar_handle = bar(avg_all, 'FaceColor', 'flat');  % Set 'FaceColor' to 'flat'
    
        % Assign colors based on speed
        for i = 1:num_bars
            speed = speed_levels(i); % Use pre-defined speed levels
            for j = 1:length(fns.speed.(speed))
                trial_idx = find(strcmp(labels, strrep(fns.speed.(speed){j}, 'Trial00', '')));
                bar_handle.CData(trial_idx, :) = colors(i, :); % Assign color
            end
        end
    
        % Customize plot
        set(gca, 'XTick', 1:length(labels), 'XTickLabel', labels, 'FontSize', 15);
        ylabel('Average Step Error (cm)', 'FontSize', 15);
        xlabel('Trial Number', 'FontSize', 15);
    
        % Create legend manually
        p = arrayfun(@(x) plot(NaN, NaN, 'Color', colors(x, :)), 1:num_bars);
        legend(p, speed_levels);
        title('Average Lateral Step Error Grouped by Speed Condition - ' + bmh_text, 'FontSize', 20);
        hold off;
    end
    
    % BAR PLOT
 
    means = zeros(1, 3); % Now only 3 categories (Slow, Medium, Fast)
    stds = zeros(1, 3);
    
    % Loop through the speeds and compute mean and std for each (without Unknown)
    for i = 1:num_bars % Loop through each speed
        speed = speed_levels(i); % Use pre-defined speed levels
        speed_data.(speed) = []; % Initialize array to collect all data for the current speed
    
        for j = 1:length(fns.speed.(speed)) % Loop through trials for each speed
            trialname = fns.speed.(speed){j};
            if participants == "multiple"
                trial_num = str2double(extractBetween(trialname, 7, 9));  % Extract the characters in positions 7 to 9 (trial number)
                bmh = regexp(trialname, 'BMH\d{2}', 'match', 'once');
            elseif participants == "single"
                trial_num = str2double(strrep(trialname, 'Trial00', '')); % Extract trial number
            end

            if trial_num >= 3 && trial_num <= 29 %As BMHs alternate, call the baseline function in each repetition and subtract what corresponds to that Trial and bmh
                [speed_s0a2b0, speed_s1a2b0, speed_s2a2b0] = BaselineErrors(imported_data,bmh,n,section); % Extract the baseline for this participant
                % When susbtracting the error, substract to each individual step and then do the average of each trial
                if subtract_all == 1 && (n == 3 || n == 4)
                    val = (speed_s0a2b0 + speed_s1a2b0 +speed_s2a2b0) / 3; % In mm
                    if n == 3
                        data = abs((error.byspeed.(speed).(trialname)) - val) / 10; % Convert from mm to cm. For all steps
                    elseif n == 4
                        data = (error.byspeed.(speed).(trialname) + val) / 10; % Convert from mm to cm. For all steps
                    end                
                elseif subtract_speeds == 1 && (n == 3 || n == 4)
                    if i == 1 % Slow speed
                        val = speed_s0a2b0;
                    elseif i == 2 % Medium speed
                        val = speed_s1a2b0;
                    elseif i == 3 % Fast speed
                        val = speed_s2a2b0;
                    end
                    if n == 3
                        data = abs((error.byspeed.(speed).(trialname)) - val) / 10; % Convert from mm to cm
                    elseif n == 4
                        data = (error.byspeed.(speed).(trialname) + val) / 10; % Convert from mm to cm
                    end
                else % Nothing is substracted
                    data = error.byspeed.(speed).(trialname) / 10; % Convert from mm to cm
                end 
                speed_data.(speed) = [speed_data.(speed); data]; % Append current trial data
            end
        end
        % Compute mean and std for the current speed
        means(i) = mean(speed_data.(speed)); % Mean of each trial
        stds(i) = std(speed_data.(speed));
    end
    
    % Substraction after averaging each trial 
    % if subtract_all == 1 && (n == 3 || n == 4)
    %     warning('Both feet used; not implemented yet for just right or left foot');
    %     %val = imported_data.(bmh).overall_errors.balance0.(m_path)/10; % convert to cm
    %     val = (speed_s0a2b0 + speed_s1a2b0 +speed_s2a2b0) / 3;
    %     if n == 3 %absolute
    %         means = means - val/10; % convert to cm
    %     elseif n == 4 %signed
    %         means = means + val/10; % convert to cm
    %     end
    %     stds = std(means);
    % 
    % elseif subtract_speeds == 1 && (n== 3 || n == 4)
    %     warning('Both feet used; not implemented yet for just right or left foot');
    %     % Low
    %     val = speed_s0a2b0;
    %     if n == 3
    %         means(1) = means(1) - val/10; 
    %     elseif n == 4
    %         means(1) = means(1) + val/10; 
    %     end      
    % 
    %     % Medium
    %     val = speed_s1a2b0;
    %     if n == 3
    %         means(2) = means(2) - val/10; 
    %     elseif n == 4
    %         means(2) = means(2) + val/10; 
    %     end
    % 
    %     % High
    %     val = speed_s2a2b0;
    %     if n == 3
    %         means(3) = means(3) - val/10;   
    %     elseif n == 4
    %         means(3) = means(3) + val/10;   
    %     end   
    % end

    % Create the bar graph with error bars
    figure;
    b = bar(means); % Bar graph of means
    b.FaceColor = 'flat'; % Enable color customization
    for k = 1:num_bars
        b.CData(k, :) = colors(k, :); % Assign a different color to each bar
    end
    hold on;

    % Ensure error bars do not go below zero
    lower_error = max(0, means - stds); % Clip lower error to 0
    upper_error = means + stds;         % Calculate upper error
    errorbar(1:num_bars, means, means - lower_error, upper_error - means, 'k', 'linestyle', 'none'); % Add asymmetric error bars
    
    % Customize plot
    set(gca, 'XTick', 1:num_bars, 'XTickLabel', {'Slow', 'Medium', 'Fast'}, 'FontSize', 15);
    ylabel('Step Error (cm)', 'FontSize', 15);
    xlabel('Speed Category', 'FontSize', 15);
    
    title({'Step Error Grouped by Speed Category - ' + bmh_text + ' - ' + n_text; section_text + ' - ' + foot_text + subtract_text}, 'FontSize', 20);

    % RAINCLOUD PLOT
    if raincloud == 1
        % Define the reduction factor (e.g., 2 to reduce points by half)
        reduction_factor = 8;
        
        % Pre-allocate a cell array to hold all speed data for the raincloud plot
        all_speed_data = cell(1, num_bars); % Store data for each speed level
        
        % Loop through the speeds and collect all data for the raincloud plot
        for i = 1:num_bars
            speed = speed_levels(i); % Use pre-defined speed levels
            speed_data = []; % Initialize array to collect all data for the current speed
        
            for j = 1:length(fns.speed.(speed)) % Loop through trials for each speed
                trialname = fns.speed.(speed){j};
                if participants == "multiple"
                    trial_num = str2double(extractBetween(trialname, 7, 9));  % Extract the characters in positions 7 to 9 (trial number)
                else
                    trial_num = str2double(strrep(trialname, 'Trial00', '')); % Extract trial number
                end
        
                if trial_num >= 3 && trial_num <= 29
                    data = error.byspeed.(speed).(trialname) / 10; % Convert from mm to cm
                    speed_data = [speed_data; data]; % Append current trial data
                end
            end
            
            % Apply the reduction factor to the speed data
            if ~isempty(speed_data)
                all_speed_data{i} = speed_data / reduction_factor; % Reduce data by the specified factor
            else
                all_speed_data{i} = []; % Handle case where no data is available
            end
        end
        
        % Create the raincloud plot
        figure;
        hold on;
        
        % Define offset to separate rainclouds for each speed category
        x_offset = 0.2; % Adjust to separate the plots horizontally
        
        for i = 1:num_bars
            if ~isempty(all_speed_data{i}) % Check if there is data to plot
                % Generate a narrower "Gaussian" effect with kernel density estimation
                [f, xi] = ksdensity(all_speed_data{i}, 'Bandwidth', 0.2); % Set bandwidth for distribution
                f = f / max(f) * 0.1; % Scale to a smaller size (height)
        
                % Plot Gaussian "cloud" on one side of the bar
                fill(f + i * x_offset, xi, colors(i, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        
                % Overlay jittered individual data points
                scatter(i * x_offset * ones(size(all_speed_data{i})) + 0.05 * (rand(size(all_speed_data{i})) - 0.5), ...
                    all_speed_data{i}, 15, colors(i, :), 'filled', 'MarkerFaceAlpha', 0.6);
        
                % Plot boxplot for central tendency and spread
                boxplot(all_speed_data{i}, 'positions', i * x_offset, 'Widths', 0.05, 'Colors', 'k', 'Symbol', '');
            end
        end
        
        % Adjust plot spacing for clearer visualization
        set(gca, 'XTick', x_offset * (1:num_bars), 'XTickLabel', {'Slow', 'Medium', 'Fast'}, 'FontSize', 15);
        ylabel('Step Error (cm)', 'FontSize', 15);
        xlabel('Speed Category', 'FontSize', 15);
        title({'Step Error Grouped by Speed Category - ' + bmh_text + ' - ' + n_text; section_text + ' - ' + foot_text + subtract_text}, 'FontSize', 20);
        
        xlim([0.5 * x_offset, (num_bars + 0.5) * x_offset]); % Adjust limits for better spacing
        hold off;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BALANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure();
    hold on;
    boxnum = 1;
    labels = {}; % Initialize a cell array to store labels
    balance_levels = ["None", "Medium", "High", "Unknown"]; % Balance categories
    
    for i = 1:num_bars % Loop through balance levels
        balance = balance_levels(i); % Use pre-defined balance levels
        fns.balance.(balance) = fieldnames(error.bybalance.(balance));
    
        for j = 1:length(fns.balance.(balance)) % Loop through trials for each balance level
            trialname = fns.balance.(balance){j};
            data = error.bybalance.(balance).(trialname); % Get data
            boxplot(data / 10, 'Positions', boxnum, 'Widths', 0.2, 'Colors', colors(i, :), 'Symbol', '_k'); % it plots here
            labels{boxnum} = strrep(trialname, 'Trial00', ''); % Store the label for the current box plot
            boxnum = boxnum + 1;
        end
    end
    
    % Create legend manually
    p = arrayfun(@(x) plot(NaN, NaN, 'Color', colors(x, :)), 1:4);
    legend(p, balance_levels(1:num_bars));
    set(gca, 'XTick', 1:boxnum-1, 'XTickLabel', labels, 'FontSize', 15);
    ylabel('Step Error (cm)', 'FontSize', 15);
    xlabel('Trial Number', 'FontSize', 15);
    
    title({'Step Error Grouped by Balance Prompt - ' + bmh_text + ' - ' + n_text; section_text + ' - ' + foot_text}, 'FontSize', 20);
    hold off;
    
    if n == 6
        figure;
        hold on;
    
        avg_all = []; % Initialize for overall averages
        labels = {}; % To store trial labels
        boxnum = 1; % To track x-axis position
    
        % Loop through balance levels
        for i = 1:num_bars
            balance = balance_levels(i); % Use pre-defined balance levels
    
            for j = 1:length(fns.balance.(balance))
                trialname = fns.balance.(balance){j};
                data = error.bybalance.(balance).(trialname); % Get data for the current trial
                avg_all(end+1) = mean(data / 10); % Convert from section to cm
    
                labels{boxnum} = strrep(trialname, 'Trial00', ''); % Store label
                boxnum = boxnum + 1;
            end
        end
    
        % Plot the bar graph for the overall averages
        bar_handle = bar(avg_all, 'FaceColor', 'flat');  % Set 'FaceColor' to 'flat'
    
        % Assign colors based on balance
        for i = 1:num_bars
            balance = balance_levels(i); % Use pre-defined balance levels
            for j = 1:length(fns.balance.(balance))
                trial_idx = find(strcmp(labels, strrep(fns.balance.(balance){j}, 'Trial00', '')));
                bar_handle.CData(trial_idx, :) = colors(i, :); % Assign color
            end
        end
    
        % Customize plot
        set(gca, 'XTick', 1:length(labels), 'XTickLabel', labels, 'FontSize', 15);
        ylabel('Average Step Error (cm)', 'FontSize', 15);
        xlabel('Trial Number', 'FontSize', 15);
    
        % Create legend manually
        p = arrayfun(@(x) plot(NaN, NaN, 'Color', colors(x, :)), 1:4);
        legend(p, balance_levels(1:num_bars));
        title('Average Lateral Step Error Grouped by Balance Condition - ' + bmh_text, 'FontSize', 20);
        hold off;
    end
    
    % BAR PLOT
    
    % Pre-allocate arrays for means and standard deviations
    means = zeros(1, 3); % Now only 3 categories (Low, Medium, High)
    stds = zeros(1, 3);
    
    % Loop through the balance levels and compute mean and std for each (without Unknown)
    for i = 1:3
        balance = balance_levels(i); % Use pre-defined balance levels
        balance_data.(balance) = []; % Initialize array to collect all data for the current balance level
    
        for j = 1:length(fns.balance.(balance)) % Loop through trials for each balance level
            trialname = fns.balance.(balance){j};
            if participants == "multiple"
                trial_num = str2double(extractBetween(trialname, 7, 9));  % Extract the characters in positions 7 to 9 (trial number)
            else
                trial_num = str2double(strrep(trialname, 'Trial00', '')); % Extract trial number
            end
    
            if trial_num >= 3 && trial_num <= 29
                data = error.bybalance.(balance).(trialname) / 10; % Convert from mm to section
                balance_data.(balance) = [balance_data.(balance); data]; % Append current trial data
            end
        end
    
        % Compute mean and std for the current balance level
        means(i) = mean(balance_data.(balance));
        stds(i) = std(balance_data.(balance));
    end
    
    % It makes no sense here subtracting the accuracy based on the speed, since here we are plotting the data regarding different accuracies

    % Create the bar graph with error bars
    figure;
    b = bar(means); % Bar graph of means
    b.FaceColor = 'flat'; % Enable color customization
    for k = 1:3
        b.CData(k, :) = colors(k, :); % Assign a different color to each bar
    end
    hold on;
    errorbar(1:3, means, stds, 'k', 'linestyle', 'none'); % Add error bars
    
    % Customize plot
    set(gca, 'XTick', 1:3, 'XTickLabel', {'Low', 'Medium', 'High'}, 'FontSize', 15);
    ylabel('Step Error (cm)', 'FontSize', 15);
    xlabel('Balance Category', 'FontSize', 15);
    
    title({'Step Error Grouped by Balance Category - ' + bmh_text + ' - ' + n_text; section_text + ' - ' + foot_text}, 'FontSize', 20);



    % RAINCLOUD PLOT
    if raincloud == 1
        % Define the reduction factor (e.g., 2 to reduce points by half)
        reduction_factor = 2;
        
        % Pre-allocate a cell array to hold all balance data for the raincloud plot
        all_balance_data = cell(1, 3); % Store data for each balance level
        
        % Loop through the balance levels and collect all data for the raincloud plot
        for i = 1:3
            balance = balance_levels(i); % Use pre-defined balance levels
            balance_data.(balance) = []; % Initialize array to collect all data for the current balance level
        
            for j = 1:length(fns.balance.(balance)) % Loop through trials for each balance level
                trialname = fns.balance.(balance){j};
                if participants == "multiple"
                    trial_num = str2double(extractBetween(trialname, 7, 9));  % Extract the characters in positions 7 to 9 (trial number)
                else
                    trial_num = str2double(strrep(trialname, 'Trial00', '')); % Extract trial number
                end
        
                if trial_num >= 3 && trial_num <= 29
                    data = error.bybalance.(balance).(trialname) / 100; % Convert from mm to cm
                    balance_data.(balance) = [balance_data; data]; % Append current trial data
                end
            end
            
            % Apply the reduction factor to the balance data
            if ~isempty(balance_data.(balance))
                all_balance_data{i} = balance_data.(balance)(1:reduction_factor:end); % Reduce points by taking every nth point
            else
                all_balance_data{i} = []; % Handle case where no data is available
            end
        end
        
        % Create the raincloud plot
        figure;
        hold on;
        
        % Define offset to separate rainclouds for each balance category
        x_offset = 0.2; % Adjust to separate the plots horizontally
        
        for i = 1:3
            if ~isempty(all_balance_data{i}) % Check if there is data to plot
                % Create the histogram "cloud" distribution
                [counts, edges] = histcounts(all_balance_data{i}, 'BinWidth', 0.2); % Bin width for each balance level
                bin_centers = edges(1:end-1) + diff(edges) / 2; % Calculate bin centers for plotting
                counts = counts / max(counts) * 0.1; % Normalize to fit height
        
                % Plot "cloud" as filled area on one side only (right side of the axis)
                fill([i * x_offset + counts, i * x_offset * ones(size(counts))], ...
                     [bin_centers, fliplr(bin_centers)], colors(i, :), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        
                % Overlay jittered individual data points
                scatter(i * x_offset * ones(size(all_balance_data{i})) + 0.05 * (rand(size(all_balance_data{i})) - 0.5), ...
                    all_balance_data{i}, 15, colors(i, :), 'filled', 'MarkerFaceAlpha', 0.6);
        
                % Plot boxplot for central tendency and spread
                boxplot(all_balance_data{i}, 'positions', i * x_offset, 'Widths', 0.05, 'Colors', 'k', 'Symbol', '');
            end
        end
        
        % Adjust plot spacing for clearer visualization
        set(gca, 'XTick', x_offset * (1:3), 'XTickLabel', {'Low', 'Medium', 'High'}, 'FontSize', 15);
        ylabel('Step Error (cm)', 'FontSize', 15);
        xlabel('Balance Category', 'FontSize', 15);
        title({'Step Error Grouped by Balance Category - ' + bmh_text + ' - ' + n_text; section_text + ' - ' + foot_text}, 'FontSize', 20);
        
        xlim([0.5 * x_offset, (3 + 0.5) * x_offset]); % Adjust limits for better spacing
        hold off;
    end

    %% All in one plot - Step Error - Plot 7

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
    
    % Compute averages and standard deviations for Accuracy
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
    
    set(gca, 'FontSize', 16);
    xticks(1:num_groups);
    xticklabels({'Balance Perturbation', 'Speed', 'Accuracy'}); % X-axis labels
    %xlabel('Condition');
    ylabel('Step Error (cm)');
    legend({'None/Low', 'Medium', 'High'});
    title('Step Error by Condition');
    grid on;
    hold off;

end %End of the function