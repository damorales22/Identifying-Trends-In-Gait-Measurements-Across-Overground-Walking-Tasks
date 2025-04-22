function[] =  Error_single_trial(imported_data,n,section,foot,trial_analyze,subtract_speeds,subtract_all,bmh,speed_s0a2b0,speed_s1a2b0,speed_s2a2b0)

    if section == 1
        right_data = imported_data.(bmh).matches.(trial_analyze).right.match;
        left_data = imported_data.(bmh).matches.(trial_analyze).left.matcht;
    elseif section == 2
        right_data = imported_data.(bmh).matches.(trial_analyze).right.match_straight;
        left_data = imported_data.(bmh).matches.(trial_analyze).left.match_straight;
    elseif section == 3 
        right_data = imported_data.(bmh).matches.(trial_analyze).right.match_curves;
        left_data = imported_data.(bmh).matches.(trial_analyze).left.match_curves;
    end
    
    trial_num = str2double(strrep(trial_analyze, 'Trial00', '')); % Extract trial number
    speed = imported_data.(bmh).prompt.speed(trial_num,1);
    balance = imported_data.(bmh).prompt.balance(trial_num,1);
    accuracy = imported_data.(bmh).prompt.accuracy(trial_num,1);

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
    
    %% Right and left foot plot
    % hold on
    % scatter(1:length(right_data(:,1)), right_data(:,n), 'filled');
    % scatter(1:length(left_data(:,1)), left_data(:,n), 'filled');
    % hold off
    % 
    % legend('Right foot', 'Left foot');
    % xlabel('Step Number');
    % ylabel('Step Error (section)');
    % title("Step Error per Step - " + trial_analyze + " - Straight Sections");
    % title({"Step Error per Step - " + trial_analyze + ' - ' + n_text + ' - ' + section_text}, 'FontSize', 20);

    %% Bar plot
    
    max_combined_data = 0;  
    min_combined_data = 0;  % Track the minimum value

    if foot == 1
        combined_data = [right_data(:, n); left_data(:, n)];
    elseif foot == 2
        combined_data = right_data(:, n);
    elseif foot == 3
        combined_data = left_data(:, n);
    end

    subtract_text = "";

    if subtract_all == 1 && n == 4 && subtract_speeds == 0 %substract the mean of all speeds when 320
        val = imported_data.(bmh).overall_errors.balance0.(m_path);
        combined_data = combined_data - val;
        subtract_text = " - Mean subtracted";
    elseif subtract_speeds == 1 && n == 4 %substract depending on the speed (120, 220, or 320)
        if speed == "Slow"
            val = speed_s0a2b0;
            combined_data = combined_data - val;
            subtract_text = " - Mean Subtracted According to Speed";
        elseif speed == "Medium"
            val = speed_s1a2b0;
            combined_data = combined_data - val;
            subtract_text = " - Mean Subtracted According to Speed";
        elseif speed == "Fast"
            val = speed_s2a2b0;
            combined_data = combined_data - val;
            subtract_text = " - Mean Subtracted According to Speed";
        end
    end

    % Update global maximum and minimum values for x-axis scaling
    max_combined_data = max(max_combined_data, max(combined_data));  % Update the global maximum
    min_combined_data = min(min_combined_data, min(combined_data));  % Update the global minimum
    
    % Define bar width and bin edges using min/max scaling similar to the previous code
    width = 20;
    bin_edges = floor(min_combined_data/20)*20:width:ceil(max_combined_data/20)*20;
    
    % Create figure and plot histogram
    figure;
    hold on;
    histogram(combined_data, 'BinEdges', bin_edges, 'FaceColor', 'b', 'EdgeColor', 'k', 'DisplayName', 'Combined Data', 'FaceAlpha', 0.7);
    
    % If n == 4, fit and plot a Gaussian curve
    if n == 4
        % Fit a Gaussian distribution
        mu = mean(combined_data);  % Mean of the data
        sigma = std(combined_data);  % Standard deviation of the data
    
        % Create a range of x values to plot the Gaussian curve
        x_values = linspace(min(bin_edges), max(bin_edges), 100);
        gaussian_curve = (1/(sigma*sqrt(2*pi))) * exp(-(x_values-mu).^2 / (2*sigma^2));
    
        % Scale the Gaussian curve to match the histogram
        gaussian_curve = gaussian_curve * (length(combined_data) * width);  % Scale based on total data count and bin width
    
        % Plot the Gaussian curve
        plot(x_values, gaussian_curve, 'r', 'LineWidth', 2, 'DisplayName', 'Gaussian Fit');
        plot(mu, max(gaussian_curve), 'p', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'y', 'DisplayName', 'Gaussian Max'); %plot a star at the maximum
    end
    
    hold off;
    
    % Set title, labels, and x-axis limits
    xlabel('Step Error (mm)');
    ylabel('Number of Counts');   
    title({bmh + ' - ' + trial_analyze + ' - ' + n_text + ' - ' + section_text + ' - ' + foot_text + subtract_text; "Speed: " + speed + " - Accuracy: " + accuracy + " - Balance disturbance: " +  balance}, 'FontSize', 20);

    % Set x-axis limits and labels with 0 ± 20 scale
    xlim([floor(min_combined_data/20)*20 ceil(max_combined_data/20)*20]);  % Ensure x-limits match the bin edges

    if n == 4 %positive and negative values
        ylim([0 20]);
        xline(0, 'k--', 'LineWidth', 2, 'DisplayName', 'x = 0'); % Draw a vertical line at x = 0
    else
        ylim([0 25]);
    end
    
    xticks(bin_edges);  % Set the x-ticks to match the bin edges
    xticklabels(string(bin_edges));  % Convert bin edges to string for labeling

end