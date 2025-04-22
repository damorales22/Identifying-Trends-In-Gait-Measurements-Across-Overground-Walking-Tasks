function[] = Multiple_customized(imported_data,error,n,section,foot,speed_input,accuracy_input,balance_input,subtract_speeds,subtract_all,bmh)

    if nargin < 11
        participants = "multiple";
    else
        participants = "single";
    end

    if n == 3
        n_text = "Longitudinal Error (Absolute Value)";
    elseif n == 4
        n_text = "Longitudinal Error";
    elseif n == 5
        n_text = "Total Error";
    elseif n == 6
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
    
    if foot == 1
        foot_text = "Both Feet";
    elseif foot == 2
        foot_text = "Right Foot";
    elseif foot == 3
        foot_text = "Left Foot";
    end

%% See what trials fulfill the user conditions

    common_labels = {};
    labels_accuracy = {};
    labels_balance = {};
    labels_speed = {};
    
    % Check conditions for accuracy, balance, and speed, set labels if matched
    if accuracy_input == 0
        labels_accuracy = fieldnames(error.byaccuracy.Low);
    elseif accuracy_input == 1
        labels_accuracy = fieldnames(error.byaccuracy.Medium);
    elseif accuracy_input == 2
        labels_accuracy = fieldnames(error.byaccuracy.High);       
    elseif accuracy_input == 3
        labels_accuracy = union(union(fieldnames(error.byaccuracy.Low), fieldnames(error.byaccuracy.Medium), 'stable'), fieldnames(error.byaccuracy.High), 'stable'); % All labels if accuracy_input is 3
    end
    
    if balance_input == 0
        labels_balance = fieldnames(error.bybalance.None);
    elseif balance_input == 1
        labels_balance = fieldnames(error.bybalance.Medium);
    elseif balance_input == 2
        labels_balance = fieldnames(error.bybalance.High);       
    elseif balance_input == 3
        labels_balance = union(union(fieldnames(error.bybalance.None), fieldnames(error.bybalance.Medium), 'stable'), fieldnames(error.bybalance.High), 'stable'); % All labels if balance_input is 3
    end
    
    if speed_input == 0
        labels_speed = fieldnames(error.byspeed.Slow);
    elseif speed_input == 1
        labels_speed = fieldnames(error.byspeed.Medium);
    elseif speed_input == 2
        labels_speed = fieldnames(error.byspeed.Fast);       
    elseif speed_input == 3
        labels_speed = union(union(fieldnames(error.byspeed.Slow), fieldnames(error.byspeed.Medium), 'stable'), fieldnames(error.byspeed.Fast), 'stable');
    end
    
    common_labels = intersect(labels_accuracy, intersect(labels_speed, labels_balance, 'stable'), 'stable'); %Find the common elements in the three arrays

    % Filter common_labels to retain only trials within the range 3-29
    for k = length(common_labels):-1:1 % Iterate in reverse to avoid indexing issues
        trial_name = common_labels{k};
        if participants == "single"
            trial_num = str2double(trial_name(end-1:end)); % Extract trial number from 'Trial00aa'
        elseif participants == "multiple"
            trial_num = str2double(extractBetween(trial_name, 7, 9));
        end
        if trial_num < 3 || trial_num > 29 % Remove trials outside of 3-29 range
            common_labels(k) = []; % Remove the trial from common_labels
        end
    end

%% Processing and plotting
    max_combined_data = 0;
    min_combined_data = 0;  % Track the minimum value

    conditions = struct();

    for i = 1:length(common_labels)

        trialname = common_labels{i};  

        if participants == "multiple"
            bmh = string(extractBetween(trialname, 11, 15));
            trialname = string(extractBetween(trialname, 1, 9));    
        end

        trial_num = str2double(extractBetween(trialname, 7, 9));

        if speed_input == 0
            speed = "Slow";
            condition_speed = "Speed_Slow";
        elseif speed_input == 1
            speed = "Medium";
            condition_speed = "Speed_Medium";
        elseif speed_input == 2
            speed = "Fast";
            condition_speed = "Speed_Fast";
        elseif speed_input == 3
            speed = imported_data.(bmh).prompt.speed(trial_num,1);
            condition_speed = "speed_" + speed;
        end

        if accuracy_input == 0
            accuracy = "Low";
            condition_accuracy = "Accuracy_Low";
        elseif accuracy_input == 1
            accuracy = "Medium";
            condition_accuracy = "Accuracy_Medium";
        elseif accuracy_input == 2
            accuracy = "High";
            condition_accuracy = "Accuracy_High";
        elseif accuracy_input == 3
            accuracy = imported_data.(bmh).prompt.accuracy(trial_num,1);
        end

        if balance_input == 0
            balance = "None";
            condition_balance = "Balance_None";
        elseif balance_input == 1
            balance = "Medium";
            condition_balance = "Balance_Medium";
        elseif balance_input == 2
            balance = "High";
            condition_balance = "Balance_High";
        elseif balance_input == 3
            balance = imported_data.(bmh).prompt.balance(trial_num,1);
        end

        right_data = imported_data.(bmh).matches.(trialname).right.match_straight;
        left_data = imported_data.(bmh).matches.(trialname).left.match_straight;

        % Use right foot, left foot, or both for each trial in the loop
        if foot == 1
            combined_data = [right_data(:, n); left_data(:, n)];
        elseif foot == 2
            combined_data = right_data(:, n);
        elseif foot == 3
            combined_data = left_data(:, n);
        end

        [speed_s0a2b0, speed_s1a2b0, speed_s2a2b0] = BaselineErrors(imported_data,bmh,n,section); % Extract the baseline for this participant
        
        %Cambiar
        subtract_text = "";
        if subtract_all == 1 && n == 4 && subtract_speeds == 0 %substract the mean of all speeds when 320
            val = imported_data.(bmh).overall_errors.balance0.(m_path);
            combined_data = combined_data - val;
            subtract_text = " - Mean Subtracted";
        elseif subtract_speeds == 1 && n == 4 %substract depending on the speed (120, 220, or 320)
            if speed == "Slow"
                val = speed_s0a2b0;
                combined_data = abs(combined_data - val);
                subtract_text = " - Mean Subtracted According to Speed";
            elseif speed == "Medium"
                val = speed_s1a2b0;
                combined_data = abs(combined_data - val);
                subtract_text = " - Mean Subtracted According to Speed";
            elseif speed == "Fast"
                val = speed_s2a2b0;
                combined_data = abs(combined_data - val);
                subtract_text = " - Mean Subtracted According to Speed";
            end
        end

        % Store data to be used in the plot across participants
        conditions.(condition_speed).(condition_accuracy).(condition_balance).(bmh) = combined_data;

        % Update global maximum and minimum values for x-axis scaling
        max_combined_data = max(max_combined_data, max(combined_data));  % Update the global maximum
        min_combined_data = min(min_combined_data, min(combined_data));  % Update the global minimum

        %Calculate means for control
        common_labels_parameters.(trialname).mean_error = mean(combined_data);

        % Define bar width and bin edges with a scale of 0 ± 20
        width = 20;  % Step size of 20
        if bmh == "BMH06" && trialname == "Trial0016" && n == 4
            bin_edges = (-360:20:360);
        else
            bin_edges = floor(min_combined_data/20)*20:width:ceil(max_combined_data/20)*20;
        end

        % Create a new figure for the histogram
        figure;
        hold on;
        histogram(combined_data, 'BinEdges', bin_edges, 'FaceColor', 'b', 'EdgeColor', 'k', ...
            'DisplayName', 'Combined Data', 'FaceAlpha', 0.7);

        if n == 4 && condition_accuracy ~= "Accuracy_Low" % fit and plot a Gaussian curve
            % Fit a Gaussian distribution
            mu = mean(combined_data);  % Mean of the data
            sigma = std(combined_data);  % Standard deviation of the data

            common_labels_parameters.(trialname).gaussian.mu = mu;
            common_labels_parameters.(trialname).gaussian.sigma = sigma;

            % Create a range of x values to plot the Gaussian curve
            
            x_values = linspace(min(bin_edges), max(bin_edges), 100);
            gaussian_curve = (1/(sigma*sqrt(2*pi))) * exp(-(x_values-mu).^2 / (2*sigma^2));

            % Scale the Gaussian curve to match the histogram
            gaussian_curve = gaussian_curve * (length(combined_data) * width);  % Scale based on total data count and bin width

            % Plot the Gaussian curve
            plot(x_values, gaussian_curve, 'r', 'LineWidth', 2, 'DisplayName', 'Gaussian Fit');
            plot(mu, max(gaussian_curve), 'p', 'MarkerSize', 18, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'DisplayName', 'Gaussian Max'); %plot a star at the maximum
        end

        % Set the title and labels for each plot
        xlabel('Step Error (cm)');
        ylabel('Number of Steps');

        % Set x-axis limits and labels with 0 ± 20 scale
        xlim([floor(min_combined_data/20)*20 ceil(max_combined_data/20)*20]);  % Ensure x-limits match the bin edges

        if n == 4 %positive and negative values
            ylim([0 20]);
            xline(0, 'k--', 'LineWidth', 2, 'DisplayName', 'x = 0'); % Draw a vertical line at x = 0
        else
            ylim([0 25]);
        end

        if n == 4 && condition_accuracy == "Accuracy_Low"
            xlim([-460 460])
            ylim([0 10])
        end

        xticks(bin_edges);  % Set the x-ticks to match the bin edges
        xticklabels(string(bin_edges/10));  % Convert bin edges to string for labeling; convert from mm to cm
        set(gca, 'FontSize', 18);

        % For a nicer plot
        % if n==4
        %     title("Step Error Distribution (Signed Error)")
        %     xlim([-400 400])
        % elseif n==3
        %     title("Step Error Distribution (Absolute Error)")
        %     xlim([0 400])
        %     ylim([0 25])
        % end

        if bmh == "BMH10" && trialname == "Trial0018"
            if n == 4
                title("Signed Step Error with High Accuracy and No Visual Disturbance")
                xlim([-400 400])
                ylim([0 20])
            elseif n == 3
                title("Absolute Step Error with High Accuracy and No Visual Disturbance")
            end
        elseif bmh == "BMH10" && trialname == "Trial0024"
            if n == 4
                title("Signed Step Error with High Accuracy and No Visual Disturbance")
                xlim([-400 400])
                ylim([0 18])
            elseif n == 3
                ylim([0 18])
                title("Absolute Step Error with High Accuracy and No Visual Disturbance")
            end    
        elseif bmh == "BMH02" && trialname == "Trial0020"
            if n == 4
                title("Signed Step Error with High Accuracy and No Visual Disturbance")
                xlim([-380 380])
                ylim([0 15])
            elseif n == 3
                title("Absolute Step Error with High Accuracy and No Visual Disturbance")
                xlim([0 380])
                ylim([0 21])
            end 
        elseif bmh == "BMH06" && trialname == "Trial0016"
            if n == 4
                title("Signed Step Error with High Accuracy and No Visual Disturbance")
                xlim([-360 360])
                ylim([0 21])
            elseif n == 3
                title("Absolute Step Error with High Accuracy and No Visual Disturbance")
            end 
        else
            title({bmh + ' - ' + string(trialname) + ' - ' + n_text + ' - ' + section_text + ' - ' + foot_text + subtract_text; "Speed: " + speed + " - Accuracy: " + accuracy + " - Visual Disturbance: " +  balance}, 'FontSize', 20);
        end
        hold off;
    end

%% Global plot for each unique condition across all trials
    if participants == "multiple"
        close all
        
        fieldnames_conditions = fieldnames(conditions); % Retrieve field names for each condition

        % Loop through each speed, accuracy, and balance combination
        for spd = 1:length(fieldnames_conditions)
            speed_field = fieldnames_conditions{spd};
            accuracy_fields = fieldnames(conditions.(speed_field));

            for acc = 1:length(accuracy_fields)
                accuracy_field = accuracy_fields{acc};
                balance_fields = fieldnames(conditions.(speed_field).(accuracy_field));

                for bal = 1:length(balance_fields)
                    balance_field = balance_fields{bal};
                    combined_data = []; % Aggregate data for the current speed-accuracy-balance combination
                    bmh_fields = fieldnames(conditions.(speed_field).(accuracy_field).(balance_field));

                    for bmh_idx = 1:length(bmh_fields)
                        bmh_field = bmh_fields{bmh_idx};
                        combined_data = [combined_data; conditions.(speed_field).(accuracy_field).(balance_field).(bmh_field)];
                    end

                    % Plot only if there is data for this condition
                    if ~isempty(combined_data)
                        % Define bin edges based on the combined min and max data for this condition
                        bin_edges = floor(min_combined_data/20)*20 : 20 : ceil(max_combined_data/20)*20;

                        figure;
                        hold on;
                        histogram_data = histogram(combined_data, 'BinEdges', bin_edges, 'FaceColor', 'b', 'EdgeColor', 'k', 'DisplayName', 'Condition Combined Data', 'FaceAlpha', 0.7);

                        max_y_value = max(histogram_data.Values);  % Calculate dynamic y-limit based on the highest bin count
                        y_limit = ceil(max_y_value / 5) * 5;  % Round up to the nearest 5 for a cleaner y-axis limit

                        if n == 4 % Gaussian
                            mu_condition = mean(combined_data); 
                            sigma_condition = std(combined_data);  

                            x_values = linspace(min(bin_edges), max(bin_edges), 100);
                            gaussian_curve = (1/(sigma_condition*sqrt(2*pi))) * exp(-(x_values - mu_condition).^2 / (2*sigma_condition^2));
                            gaussian_curve = gaussian_curve * (length(combined_data) * 20);  % Scale to histogram

                            plot(x_values, gaussian_curve, 'r', 'LineWidth', 2, 'DisplayName', 'Gaussian Fit');
                            plot(mu_condition, max(gaussian_curve), 'p', 'MarkerSize', 18, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'DisplayName', 'Gaussian Max');
                            
                            xline(0, 'k--', 'LineWidth', 2, 'DisplayName', 'x = 0'); %  x = 0
                        end

                        xlabel('Step Error (cm)');
                        ylabel('Number of Steps');

                        speed_label = extractAfter(speed_field, "_");
                        accuracy_label = extractAfter(accuracy_field, "_");
                        balance_label = extractAfter(balance_field, "_");
                        
                        if accuracy_input == 0 && balance_input == 0 && n==3
                            %title("Accumulated Absolute Error for Medium Speed Prompt, No Visual Disturbance, and Null Accuracy Prompt")
                            title('Accumulated Step Error - Medium Speed, No Visual Disturbance, Null Accuracy')
                            xlim([0 460])
                        else
                            title({["Across all Participants - " + n_text + ' - ' + section_text + ' - ' + foot_text + subtract_text]; "Speed: " + speed_label + " - Accuracy: " + accuracy_label + " - Visual Disturbance: " + balance_label}, 'FontSize', 20);
                            xlim([floor(min_combined_data/20)*20, ceil(max_combined_data/20)*20]);
                        end

                        ylim([0, y_limit]);
                        xticks(bin_edges);
                        xticklabels(string(bin_edges/10));
                        %legend('show');
                        set(gca, 'FontSize', 18);
                        hold off;
                    end
                end
            end
        end
    end
end

