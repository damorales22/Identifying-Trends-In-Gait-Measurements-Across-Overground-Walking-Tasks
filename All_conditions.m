function [] = All_conditions(data_labels,plot_number,subtract_all,subtract_speeds,baseline_errors,show_cluster_analysis,participants,show_line,show_formula,plot_option,normalize_metabolics,bmh)

    if participants == "multiple"
        bmh_text = "Across all Participants";
    else
        bmh_text = string(bmh);
    end

    list_labels = fieldnames(data_labels); %labels stored in the struct

    %% Colors and shapes
    balance_codes = zeros(1, 27);
    accuracy_codes = zeros(1, 27);
    speed_codes = zeros(1, 27);

    for i = 1:length(list_labels)
        trial = list_labels{i};
        
        % Extract the speed, accuracy, balance from the trial string (sxaybz)
        s = str2double(trial(2));  % Speed
        a = str2double(trial(4));  % Accuracy
        b = str2double(trial(6));  % Balance
        
        % According to the order stored in the struct
        accuracy_codes(i) = a; % Map accuracy to color (0=green, 1=red, 2=blue)
        balance_codes(i) = b; % Map balance to shape (0=Circle, 1=Square, 2=Triangle)  
        speed_codes(i) = s;
    end

    if plot_option == 1 % Color: accuracy. Shape: balance
        color_codes = accuracy_codes;
        shape_codes = balance_codes;
        shapes = {'o', 's', '^'};  % Circle, Square, Triangle (for balance)
        colors = {'g', 'r', 'b'};  % Green, Red, Blue (for accuracy)
        legend_labels = {'Accuracy: None, Visual Perturbation: None', 'Accuracy: None, Visual Perturbation: Medium', ...
                        'Accuracy: None, Visual Perturbation: High', 'Accuracy: Medium, Visual Perturbation: None', ...
                        'Accuracy: Medium, Visual Perturbation: Medium', 'Accuracy: Medium, Visual Perturbation: High', ...
                        'Accuracy: High, Visual Perturbation: None', 'Accuracy: High, Visual Perturbation: Medium', ...
                        'Accuracy: High, Visual Perturbation: High'};
        vars_order = ["a0b0","a0b1","a0b2","a1b0","a1b1","a1b2","a2b0","a2b1","a2b2"];
        clusters_color_silhouette = {'No Accuracy', 'Medium Accuracy', 'High Accuracy'};
    elseif plot_option == 2 % Color: balance. Shape: speed
        color_codes = balance_codes;
        shape_codes = speed_codes;
        shapes = {'p', 'd', '<'};  
        colors = {'#00B5B8', '#FF00FF', '#FFA500'};  % Cyan (#00B5B8), Magenta (#FF00FF), Orange (#FFA500)
        legend_labels = {'Visual Perturbation: None, Speed: Slow', 'Visual Perturbation: None, Speed: Medium', ...
                        'Visual Perturbation: None, Speed: Fast', 'Visual Perturbation: Medium, Speed: Slow', ...
                        'Visual Perturbation: Medium, Speed: Medium', 'Visual Perturbation: Medium, Speed: Fast', ...
                        'Visual Perturbation: High, Speed: Slow', 'Visual Perturbation: High, Speed: Medium', ...
                        'Visual Perturbation: High, Speed: Fast'};
        vars_order = ["b0s0","b0s1","b0s2","b1s0","b1s1","b1s2","b2s0","b2s1","b2s2"];
        clusters_color_silhouette = {'No Visual Perturbation', 'Medium Visual Perturbation', 'High Visual Perturbation'};
    elseif plot_option == 3 % Color: speed. Shape: accuracy
        color_codes = speed_codes;
        shape_codes = accuracy_codes;
        shapes = {'h', 'o', 'v'};  
        colors = {'#DAA520', '#FF1493', '#32CD32'};  % Gold (#DAA520), DeepPink (#FF1493), LimeGreen (#32CD32)
        legend_labels = {'Speed: Slow, Accuracy: None', 'Speed: Slow, Accuracy: Medium', ...
                        'Speed: Slow, Accuracy: High', 'Speed: Medium, Accuracy: None', ...
                        'Speed: Medium, Accuracy: Medium', 'Speed: Medium, Accuracy: High', ...
                        'Speed: Fast, Accuracy: None', 'Speed: Fast, Accuracy: Medium', ...
                        'Speed: Fast, Accuracy: High'};
        vars_order = ["s0a0","s0a1","s0a2","s1a0","s1a1","s1a2","s2a0","s2a1","s2a2"];
        clusters_color_silhouette = {'Slow Speed', 'Medium Speed', 'Fast Speed'};
    else
        error("Not valid plot option input")
    end

%% Plot
     % plot_number = ...
     % 1: speed vs metabolics.      2: speed vs step width.        3: speed vs step length.        4: speed vs step error.
     % 5: metabolics vs step width. 6: metabolics vs step length.  7: metabolics vs step error.    8: step width vs step length. 
     % 9: step width vs step error. 10: step length vs step error. 11: speed vs width variability. 12: metabolics vs width variability. 
     % 13: step error vs width variability. 14: step error vs head angle. 15: step error vs head angle (just colors)
     % 16: step error vs head angle (shape and colors). 17: step error vs trunk angle. 18: step error vs trunk angle (shape and colors). 
     % 19: speed vs stride time. 20: trunk angle vs step width. 21: trunk angle vs step width var. 22: stride time vs step length var. 
     % 23: stride time var vs step length var. 24: stride time var vs step error 25: step width var vs step length. 
     % 26: step width var vs step length var. 27: step width var vs step width

    walking_speed_avg = struct();
    walking_speed_var_avg = struct();
    metabolics_avg = struct();
    step_width_avg = struct();
    step_width_var_avg = struct();
    step_length_avg = struct();
    step_length_var_avg = struct();
    stride_time_avg = struct();
    stride_time_var_avg = struct();
    step_error_avg = struct();
    step_error_var_avg = struct();
    head_angle_avg = struct();
    trunk_angle_avg = struct();

    figure;
    hold on;

    points_a0b0 = []; % Array to store the data points of each plot
    points_a0b1 = [];
    points_a0b2 = [];
    points_a1b0 = [];
    points_a1b1 = [];
    points_a1b2 = [];
    points_a2b0 = [];
    points_a2b1 = [];
    points_a2b2 = [];
    count1 = 1;
    count2 = 1;
    count3 = 1;
    count4 = 1;
    count5 = 1;
    count6 = 1;
    count7 = 1;
    count8 = 1;
    count9 = 1;

    var1_all = [];
    var2_all = [];

    if plot_number >= 1 && plot_number <= 4 || plot_number == 11 || plot_number == 19
        for i = 1:length(list_labels)
    
            current_label = list_labels{i}; % Use curly braces {} for cell array indexing
            
            shape = shapes{shape_codes(i) + 1};  % Get shape 
            color = colors{color_codes(i) + 1};  % Get color 

            if normalize_metabolics == 0
                metabolics_avg = mean(data_labels.(current_label).metabolics.not_normalized, 1);       % Mean across all participants for that label
            elseif normalize_metabolics == 1
                metabolics_avg = mean(data_labels.(current_label).metabolics.normalized, 1);   % Mean across all participants for that label
            end

            walking_speed_avg = mean(data_labels.(current_label).walking_speed, 1); % Mean across all participants for that label
            step_width_avg = mean(data_labels.(current_label).step_width, 1);       % Mean across all participants for that label
            step_length_avg = mean(data_labels.(current_label).step_length, 1);     % Mean across all participants for that label
            step_width_var_avg = mean(data_labels.(current_label).step_length, 1);   % Mean across all participants for that label
            head_angle_avg = mean(data_labels.(current_label).head_angle, 1);  % Mean across all participants for that label
            trunk_angle_avg = mean(data_labels.(current_label).trunk_angle, 1);  % Mean across all participants for that label
            walking_speed_var_avg = nanmean(data_labels.(current_label).walking_speed_std, 1); 
            step_length_var_avg = mean(data_labels.(current_label).step_length_var, 1); 
            stride_time_avg = mean(data_labels.(current_label).stride_time, 1);
            stride_time_var_avg = mean(data_labels.(current_label).stride_time_var, 1);
            step_error_avg = mean(data_labels.(current_label).step_error, 1); 
            step_error_var_avg = mean(data_labels.(current_label).step_error_std, 1); 
            
            if subtract_all == 1
                % The bmh folders have to be sorted by their number
                for j=1:size(data_labels.(current_label).step_error,1) %To subtract the baseline of each bmh
                    bmh = bmh_list(j);
                    step_error_(j) = abs(data_labels.(current_label).step_error(j) - ((baseline_errors.(bmh).speed_s0a2b0 + baseline_errors.(bmh).speed_s1a2b0 +baseline_errors.(bmh).speed_s2a2b0)/3));
                end
                step_error_avg = mean(step_error_); % Mean across all participants - mean error
            elseif subtract_speeds == 1
                if startsWith(list_labels(i), 's0')
                    for j=1:size(data_labels.(current_label).step_error,1)
                        bmh = bmh_list(j);
                        step_error_(j) = abs(data_labels.(current_label).step_error(j) - baseline_errors.(bmh).speed_s0a2b0);
                    end
                step_error_avg = mean(step_error_); % Mean across all participants - mean error
                elseif startsWith(list_labels(i), 's1')
                    for j=1:size(data_labels.(current_label).step_error,1)
                        bmh = bmh_list(j);
                        step_error_(j) = abs(data_labels.(current_label).step_error(j) - baseline_errors.(bmh).speed_s1a2b0);
                    end
                    step_error_avg = mean(step_error_); % Mean across all participants - mean error
                elseif startsWith(list_labels(i), 's2')
                    for j=1:size(data_labels.(current_label).step_error,1)
                        bmh = bmh_list(j);
                        step_error_(j) = abs(data_labels.(current_label).step_error(j) - baseline_errors.(bmh).speed_s2a2b0);
                    end
                    step_error_avg = mean(step_error_); % Mean across all participants - mean error
                end
            else % No subtraction of the baseline
                step_error_avg = mean(data_labels.(current_label).step_error, 1);   % Mean across all participants for that label
            end

            if plot_number == 1

                %scatter(walking_speed_avg, metabolics_avg, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color); % plots one by one
                scatter(walking_speed_avg, metabolics_avg, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, ...
                'MarkerFaceAlpha', 0.12, ... % Reduce opacity
                'MarkerEdgeAlpha', 0.1);   % Reduce edge opacity

                xlabel('Walking Speed (m/s)');
                if normalize_metabolics == 0
                    ylabel('Energy Expenditure (W)');
                elseif normalize_metabolics == 1
                    ylabel('Energy Expenditure (W/kg)');
                end
                title('Energy Expenditure by Walking Speed');

                var1 = walking_speed_avg;
                var2 = metabolics_avg;
                var1_all = [var1_all; var1]; % Accumulate data for regression line
                var2_all = [var2_all; var2];


            elseif plot_number == 2
                scatter(walking_speed_avg, step_width_avg, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, ...
                'MarkerFaceAlpha', 0.12, ... % Reduce opacity
                'MarkerEdgeAlpha', 0.1);   % Reduce edge opacity

                xlabel('Walking Speed (m/s)');
                ylabel('Step Width (cm)');
                title('Step Width by Walking Speed');

                var1 = walking_speed_avg;
                var2 = step_width_avg;
                var1_all = [var1_all; var1]; % Accumulate data for regression line
                var2_all = [var2_all; var2];
                
            elseif plot_number == 3
                scatter(walking_speed_avg, step_length_avg, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, ...
                'MarkerFaceAlpha', 0.12, ... % Reduce opacity
                'MarkerEdgeAlpha', 0.1);   % Reduce edge opacity

                xlabel('Walking Speed (m/s)');
                ylabel('Step Length (cm)');
                title('Step Length by Walking Speed');

                var1 = walking_speed_avg;
                var2 = step_length_avg;
                var1_all = [var1_all; var1]; % Accumulate data for regression line
                var2_all = [var2_all; var2];
    
            elseif plot_number == 4
                scatter(walking_speed_avg, step_error_avg, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, ...
                'MarkerFaceAlpha', 0.12, ... % Reduce opacity
                'MarkerEdgeAlpha', 0.1);   % Reduce edge opacity

                xlabel('Walking Speed (m/s)');
                ylabel('Step Error (cm)');
                if subtract_all == 1
                    title('Step Error by Walking Speed - Mean Error Subtracted for All');
                elseif subtract_speeds == 1
                    title('Step Error by Walking Speed - Mean Error Subtracted for Each Speed');
                else
                    title('Step Error by Walking Speed');
                end

                var1 = walking_speed_avg;
                var2 = step_error_avg;
                var1_all = [var1_all; var1]; % Accumulate data for regression line
                var2_all = [var2_all; var2];

            elseif plot_number == 11
                scatter(walking_speed_avg, step_width_var_avg, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, ...
                'MarkerFaceAlpha', 0.12, ... % Reduce opacity
                'MarkerEdgeAlpha', 0.1);   % Reduce edge opacity

                xlabel('Walking Speed (m/s)');
                ylabel('Step Width Variability (cm)');
                title('Step Width Variability by Walking Speed');

                var1 = walking_speed_avg;
                var2 = step_width_var_avg;
                var1_all = [var1_all; var1]; % Accumulate data for regression line
                var2_all = [var2_all; var2];

                elseif plot_number == 19
                scatter(walking_speed_avg, stride_time_avg, 200, 'Marker', shape, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, ...
                'MarkerFaceAlpha', 0.12, ... % Reduce opacity
                'MarkerEdgeAlpha', 0.1);   % Reduce edge opacity

                xlabel('Walking Speed (m/s)');
                ylabel('Stride Duration (s)');
                title('Stride Duration by Walking Speed');

                var1 = walking_speed_avg;
                var2 = stride_time_avg;
                var1_all = [var1_all; var1]; % Accumulate data for regression line
                var2_all = [var2_all; var2];
            end

            % Store all points
            if shape_codes(i) == 0 && color_codes(i) == 0
                points_a0b0(count1,1) = var1;
                points_a0b0(count1,2) = var2;
                count1 = count1+1;
            elseif shape_codes(i) == 1 && color_codes(i) == 0
                points_a0b1(count2,1) = var1;
                points_a0b1(count2,2) = var2;
                count2 = count2+1;
            elseif shape_codes(i) == 2 && color_codes(i) == 0
                points_a0b2(count3,1) = var1;
                points_a0b2(count3,2) = var2;
                count3 = count3+1;
            elseif shape_codes(i) == 0 && color_codes(i) == 1
                points_a1b0(count4,1) = var1;
                points_a1b0(count4,2) = var2;
                count4 = count4+1;
            elseif shape_codes(i) == 1 && color_codes(i) == 1
                points_a1b1(count5,1) = var1;
                points_a1b1(count5,2) = var2;
                count5 = count5+1;
            elseif shape_codes(i) == 2 && color_codes(i) == 1
                points_a1b2(count6,1) = var1;
                points_a1b2(count6,2) = var2;
                count6 = count6+1;
            elseif shape_codes(i) == 0 && color_codes(i) == 2
                points_a2b0(count7,1) = var1;
                points_a2b0(count7,2) = var2;
                count7 = count7+1;
            elseif shape_codes(i) == 1 && color_codes(i) == 2
                points_a2b1(count8,1) = var1;
                points_a2b1(count8,2) = var2;
                count8 = count8+1;
            elseif shape_codes(i) == 2 && color_codes(i) == 2
                points_a2b2(count9,1) = var1;
                points_a2b2(count9,2) = var2;
                count9 = count9+1;
            end
        end

        %% Plot averages
        avg_a0b0 = mean(points_a0b0,1);
        scatter(avg_a0b0(1), avg_a0b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});
    
        avg_a0b1 = mean(points_a0b1,1);
        scatter(avg_a0b1(1), avg_a0b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});
    
        avg_a0b2 = mean(points_a0b2,1);
        scatter(avg_a0b2(1), avg_a0b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});
    
        avg_a1b0 = mean(points_a1b0,1);
        scatter(avg_a1b0(1), avg_a1b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});
    
        avg_a1b1 = mean(points_a1b1,1);
        scatter(avg_a1b1(1), avg_a1b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});
    
        avg_a1b2 = mean(points_a1b2,1);
        scatter(avg_a1b2(1), avg_a1b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2});
    
        avg_a2b0 = mean(points_a2b0,1);
        scatter(avg_a2b0(1), avg_a2b0(2), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});
    
        avg_a2b1 = mean(points_a2b1,1);
        scatter(avg_a2b1(1), avg_a2b1(2), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3});
    
        avg_a2b2 = mean(points_a2b2,1);
        scatter(avg_a2b2(1), avg_a2b2(2), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); 

        legend_handles = gobjects(9, 1);  % Preallocate graphics objects for legend entries
        
        grid on;
        ax = gca;               % Get current axis
        ax.FontSize = 20;       % Set axis number font size (e.g., 12)
        
        xlim([min(var1_all) - range(var1_all) * 0.04, max(var1_all) + range(var1_all) * 0.04]);  
        ylim([min(var2_all) - range(var2_all) * 0.04, max(var2_all) + range(var2_all) * 0.04]);  

        % Generate each combination for the legend
        if plot_number ~=15
            index = 1;
            for color = colors
                for shape = shapes
                    legend_handles(index) = scatter(NaN, NaN, 200, 'Marker', shape{1}, 'MarkerEdgeColor', color{1}, 'MarkerFaceColor', color{1});
                    index = index + 1;
                end
            end
        else
            return % dont do cluster analysis 
        end

        legend_handle = legend(legend_handles, legend_labels, 'Location', 'best');
        
        legend_handle.FontSize = 14; % size of the legend

        if show_line == 1
            % Perform linear regression
            coefficients = polyfit(var1_all, var2_all, 1); % Fit line: y = mx + c
            var2_fit = polyval(coefficients, var1_all);   % Predicted y values
            m = coefficients(1); % Slope
            c = coefficients(2); % Intercept
            
            % Compute R^2
            ss_total = sum((var2_all - mean(var2_all)).^2); % Total sum of squares
            ss_residual = sum((var2_all - var2_fit).^2);   % Residual sum of squares
            r_squared = 1 - (ss_residual / ss_total);
            
            % Plot the regression line
            %var1_range = linspace(min(var1_all)-0.01, max(var1_all)+0.01, 100); % Create range for the line
            var1_range = linspace(min(var1_all), max(var1_all), 100);
            %var1_range = linspace(min(var1_all) - range(var1_all) * 0.04, max(var1_all) + range(var1_all) * 0.04, 100); % Create range for the line
            var2_range = polyval(coefficients, var1_range);
            plot(var1_range, var2_range, 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
            
            % Display the formula and R^2 in a location away from the line
            if show_formula == 1
                formula_text = sprintf('y = %.2fx + %.2f\nR^2 = %.2f', m, c, r_squared);
                if plot_number == 1
                    text(max(var1_all)-0.08, min(var2_all)+0.12, formula_text, 'FontSize', 16);
                elseif plot_number == 2
                    text(min(var1_all), min(var2_all)+0.1, formula_text, 'FontSize', 16);
                elseif plot_number == 3
                    text(min(var1_all), max(var2_all)-0.2, formula_text, 'FontSize', 16);
                elseif plot_number == 4
                    text(min(var1_all)-0.08, min(var2_all)+0.2, formula_text, 'FontSize', 16);
                    xlim([0.65 1.4])
                elseif plot_number == 11
                    text(min(var1_all) + 0.005, max(var2_all)-0.4, formula_text, 'FontSize', 16);
                elseif plot_number == 19
                    text(min(var1_all)-0.005, min(var2_all) + 0.03, formula_text, 'FontSize', 16);
                end
            end
        end

        %% Silhouette cluster analysis
        if show_cluster_analysis == 1
            % Silhouette - By accuracy across all conditions (color)
            
            % Old way
            % clusters = [[points_a0b0; points_a0b1; points_a0b2]; [points_a1b0; points_a1b1; points_a1b2]; [points_a2b0; points_a2b1; points_a2b2]];
            % labels_silhouette = [0*ones(9,1); 1*ones(9,1); 2*ones(9,1)];
            % 
            % figure;
            % silhouette(clusters, labels_silhouette);
            % 
            % xlabel('Silhouette Value');
            % %ylabel('Cluster');
            % title('Silhouette Analysis');
            % yticklabels(clusters_color_silhouette);
            % silhouetteValues = silhouette(clusters, labels_silhouette);
            % avgSilhouette = mean(silhouetteValues);
            % disp(['Average Silhouette Score by accuracy (color): ', num2str(avgSilhouette)]);
            % ax = gca;             
            % ax.FontSize = 20;      
            % ylim([1 35])
            % xlim([-1 1])

            clusters = [[points_a0b0; points_a0b1; points_a0b2]; 
                        [points_a1b0; points_a1b1; points_a1b2]; 
                        [points_a2b0; points_a2b1; points_a2b2]];
            
            labels_silhouette = [0*ones(9,1); 1*ones(9,1); 2*ones(9,1)];
            silhouetteValues = silhouette(clusters, labels_silhouette);

            uniqueClusters = flip(unique(labels_silhouette)); 
    
            sortedSilhouette = [];
            sortedLabels = [];
            barPositions = [];
            separation = 1.25; % Distance between clusters
            currentY = 1; 
            
            for i = uniqueClusters'
                clusterIdx = find(labels_silhouette == i); % Índices de los puntos del cluster
                silValues = silhouetteValues(clusterIdx); % Valores silhouette del cluster
                
                [sortedSilValues, sortIdx] = sort(silValues, 'descend'); 
                
                sortedSilhouette = [sortedSilhouette; sortedSilValues];
                sortedLabels = [sortedLabels; i * ones(length(sortedSilValues), 1)];
                
                clusterYPositions = (currentY + length(sortedSilValues) - 1):-1:currentY;
                barPositions = [barPositions, clusterYPositions]; 
                
                currentY = currentY + length(sortedSilValues) + separation;
            end
            
            clusterSizes = histcounts(sortedLabels, 'BinMethod', 'integers');
            clusterPositions = cumsum([1, clusterSizes(1:end-1) + separation]) + clusterSizes / 2 - 0.5;
            
            % Horizontal bars
            figure;
            hold on;
            for i = 1:length(sortedSilhouette)
                clusterID = sortedLabels(i);
                color = colors{clusterID + 1}; % Obtener color del cluster
                barh(barPositions(i), sortedSilhouette(i), 'FaceColor', color, 'EdgeColor', 'none');
            end
            hold off;
            
            xlabel('Silhouette Value');
            %ylabel('Cluster');
            title('Silhouette Analysis');
            yticks(clusterPositions);
            yticklabels({clusters_color_silhouette{3}, clusters_color_silhouette{2}, clusters_color_silhouette{1}});
            
            ax = gca;
            ax.FontSize = 20;
            ylim([0 31.3]);
            xlim([-1 1]);
            
            avgSilhouette = mean(silhouetteValues);
            disp(['Average Silhouette Score by accuracy (color): ', num2str(avgSilhouette)]);
            text(0.99, 0.98, sprintf('Avg Silhouette: %.2f', avgSilhouette),'Units', 'normalized', 'FontSize', 16, 'HorizontalAlignment', 'right');


            % Silhouette - By balance across all conditions (shape)
            clusters = [[points_a0b0; points_a1b0; points_a2b0]; [points_a0b1; points_a1b1; points_a2b1]; [points_a0b2; points_a1b2; points_a2b2]];
            labels_silhouette = [0*ones(9,1); 1*ones(9,1); 2*ones(9,1)];
            
            figure;
            silhouette(clusters, labels_silhouette);
            
            xlabel('Silhouette Value');
            %ylabel('Cluster');
            title('Silhouette Analysis');
    
            silhouetteValues = silhouette(clusters, labels_silhouette);
            avgSilhouette = mean(silhouetteValues);
            disp(['Average Silhouette Score by balance (shape): ', num2str(avgSilhouette)]);
            ax = gca;             
            ax.FontSize = 20; 
    
            % Silhouette - By speed for that condition (shape + color)
            clusters = [points_a0b0; points_a1b0; points_a2b0; points_a0b1; points_a1b1; points_a2b1; points_a0b2; points_a1b2; points_a2b2];
            labels_silhouette = [1*ones(3,1); 2*ones(3,1);3*ones(3,1); 4*ones(3,1); 5*ones(3,1);6*ones(3,1); 7*ones(3,1); 8*ones(3,1); 9*ones(3,1)];
            
            figure;
            silhouette(clusters, labels_silhouette);
            
            xlabel('Silhouette Value');
            %ylabel('Cluster');
            title('Silhouette Plot by Speed for that Condition (Color + Shape)');
    
            silhouetteValues = silhouette(clusters, labels_silhouette);
            avgSilhouette = mean(silhouetteValues);
            disp(['Average Silhouette Score by speed for that condition (color + shape): ', num2str(avgSilhouette)]);
            ax = gca;             
            ax.FontSize = 20; 
    
            % Mean by accuracy (color)
            clusters = [[mean(points_a0b0); mean(points_a0b1); mean(points_a0b2)]; [mean(points_a1b0); mean(points_a1b1); mean(points_a1b2)]; [mean(points_a2b0); mean(points_a2b1); mean(points_a2b2)]];
            labels_silhouette = [0*ones(3,1); 1*ones(3,1); 2*ones(3,1)];
            
            figure;
            silhouette(clusters, labels_silhouette);
            
            xlabel('Silhouette Value');
            %ylabel('Cluster');
            title('Silhouette Plot of the Means by Accuracy (Color)');
    
            silhouetteValues = silhouette(clusters, labels_silhouette);
            avgSilhouette = mean(silhouetteValues);
            disp(['Average Silhouette Score of the means by accuracy (color): ', num2str(avgSilhouette)]);
            ax = gca;             
            ax.FontSize = 20; 
    
            % Mean by shape (balance)
            clusters = [[mean(points_a0b0); mean(points_a1b0); mean(points_a2b0)]; [mean(points_a0b1); mean(points_a1b1); mean(points_a2b1)]; [mean(points_a0b2); mean(points_a1b2); mean(points_a2b2)]];
            labels_silhouette = [0*ones(3,1); 1*ones(3,1); 2*ones(3,1)];
            
            figure;
            silhouette(clusters, labels_silhouette);
            
            xlabel('Silhouette Value');
            %ylabel('Cluster');
            title('Silhouette Plot of the Means by Balance (Shape)');
    
            silhouetteValues = silhouette(clusters, labels_silhouette);
            avgSilhouette = mean(silhouetteValues);
            disp(['Average Silhouette Score of the means by balance (shape): ', num2str(avgSilhouette)]);
            ax = gca;             
            ax.FontSize = 20; 
    
            %% Inter-Cluster vs. Intra-Cluster Density; useless here as it's not normalized, but works
            % Inter and intra for each cluster, and finally the division
    
        %     clusters = [[points_a0b0; points_a0b1; points_a0b2]; [points_a1b0; points_a1b1; points_a1b2]; [points_a2b0; points_a2b1; points_a2b2]];
        %     clusterSize = 3; % Number of points per cluster
        %     numClusters = size(clusters, 1) / clusterSize;
        % 
        %     clusterLabels = ["a0b0", "a1b0", "a2b0", "a0b1", "a1b1", "a2b1", "a0b2", "a1b2", "a2b2"];
        % 
        %     % Initialize arrays for intra- and inter-cluster densities
        %     intraClusterDensities = zeros(numClusters, 1);
        %     interClusterDensities = zeros(numClusters, 1);
        % 
        %     % Loop through each cluster
        %     for i = 1:numClusters
        %         % Extract points for the current cluster
        %         startIdx = (i - 1) * clusterSize + 1;
        %         endIdx = i * clusterSize;
        %         clusterPoints = clusters(startIdx:endIdx, :);
        % 
        %         % 1. Intra-cluster density: average pairwise distance within the cluster
        %         pairwiseDistances = pdist(clusterPoints, 'euclidean');
        %         intraClusterDensities(i) = mean(pairwiseDistances);
        % 
        %         % 2. Inter-cluster separation: minimum average distance to other clusters
        %         minInterClusterDist = inf;
        %         for j = 1:numClusters
        %             if j ~= i
        %                 % Extract points for the other cluster
        %                 otherStartIdx = (j - 1) * clusterSize + 1;
        %                 otherEndIdx = j * clusterSize;
        %                 otherClusterPoints = clusters(otherStartIdx:otherEndIdx, :);
        % 
        %                 % Compute pairwise distances between the two clusters
        %                 interDistances = pdist2(clusterPoints, otherClusterPoints, 'euclidean');
        %                 avgInterDistance = mean(interDistances(:));
        % 
        %                 % Update minimum inter-cluster distance
        %                 minInterClusterDist = min(minInterClusterDist, avgInterDistance);
        %             end
        %         end
        %         interClusterDensities(i) = minInterClusterDist;
        %     end
        % 
        %     % Compute overall density ratio (inter/intra)
        %     densityRatios = interClusterDensities ./ intraClusterDensities;
        %     overallDensityRatio = mean(densityRatios);
        % 
        %     disp('Intra-Cluster and Inter-Cluster Densities:');
        %     resultsTable = table(clusterLabels', intraClusterDensities, interClusterDensities, densityRatios, ...
        % 'VariableNames', {'Cluster', 'IntraDensity', 'InterDensity', 'DensityRatio'});
        %     disp(resultsTable);
        % 
        %     disp(['Overall Density Ratio (Inter/Intra): ', num2str(overallDensityRatio)]);
    
            % Intra-Cluster Density: Average pairwise distance (or othercompactness measure) among the points in that cluster. 
            % Lower values indicate tightly clustered points.
    
            % Inter-Cluster Separation: Minimum average distance to points in other clusters. Higher values indicate better 
            % separation from other clusters.
    
            % Density Ratio: Ratio of inter-cluster to intra-cluster density.  Higher values indicate clusters that are 
            % well-separated relative to their compactness.
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% We don't have speed in the x axis => new processing to group by speeds
    else
        for i = 1:length(list_labels)

            current_label = list_labels{i}; % Use curly braces {} for cell array indexing
            if plot_option == 1
                label_storing = extractBetween(string(current_label),3,6); % Exclude the speed dimension (for speed the previous part of the code)
            elseif plot_option == 2
                label_storing = [current_label(end-1:end), current_label(1:2)];
            elseif plot_option == 3
                label_storing = extractBetween(string(current_label),1,4);
            end

            if ~isfield(walking_speed_avg, label_storing)
                walking_speed_avg.(label_storing) = [];
            end
            if ~isfield(walking_speed_var_avg, label_storing)
                walking_speed_var_avg.(label_storing) = [];
            end
            if ~isfield(metabolics_avg, label_storing)
                metabolics_avg.(label_storing) = [];
            end
            if ~isfield(step_width_avg, label_storing)
                step_width_avg.(label_storing) = [];
            end
            if ~isfield(step_width_var_avg, label_storing)
                step_width_var_avg.(label_storing) = [];
            end
            if ~isfield(step_length_avg, label_storing)
                step_length_avg.(label_storing) = [];
            end
            if ~isfield(step_length_var_avg, label_storing)
                step_length_var_avg.(label_storing) = [];
            end
            if ~isfield(step_error_avg, label_storing)
                step_error_avg.(string(label_storing)) = [];
            end
            if ~isfield(step_error_var_avg, label_storing)
                step_error_var_avg.(string(label_storing)) = [];
            end
            if ~isfield(stride_time_avg, label_storing)
                stride_time_avg.(string(label_storing)) = [];
            end
            if ~isfield(stride_time_var_avg, label_storing)
                stride_time_var_avg.(string(label_storing)) = [];
            end
            if ~isfield(head_angle_avg, label_storing)
                head_angle_avg.(string(label_storing)) = [];
            end
            if ~isfield(trunk_angle_avg, label_storing)
                trunk_angle_avg.(string(label_storing)) = [];
            end
    
            if normalize_metabolics == 0
                metabolics_avg.(label_storing)(end+1) = mean(data_labels.(current_label).metabolics.not_normalized, 1);       % Mean across all participants for that label
            elseif normalize_metabolics == 1
                metabolics_avg.(label_storing)(end+1) = mean(data_labels.(current_label).metabolics.normalized, 1);       % Mean across all participants for that label
            end

            walking_speed_avg.(label_storing)(end+1) = mean(data_labels.(current_label).walking_speed, 1); % Mean across all participants for that label
            step_width_avg.(label_storing)(end+1) = mean(data_labels.(current_label).step_width, 1);       % Mean across all participants for that label
            step_width_var_avg.(label_storing)(end+1) = mean(data_labels.(current_label).step_width_var, 1); % Mean across all participants for that label
            step_length_avg.(label_storing)(end+1) = mean(data_labels.(current_label).step_length, 1);    % Mean across all participants for that label
            head_angle_avg.(label_storing)(end+1) = mean(data_labels.(current_label).head_angle, 1);     % Mean across all participants for that label
            trunk_angle_avg.(label_storing)(end+1) = mean(data_labels.(current_label).trunk_angle, 1);     % Mean across all participants for that label
            walking_speed_var_avg.(label_storing)(end+1) = nanmean(data_labels.(current_label).walking_speed_std, 1);
            step_length_var_avg.(label_storing)(end+1) = mean(data_labels.(current_label).step_length_var, 1);
            stride_time_avg.(label_storing)(end+1) = mean(data_labels.(current_label).stride_time, 1);
            stride_time_var_avg.(label_storing)(end+1) = mean(data_labels.(current_label).stride_time_var, 1);
            %step_error_avg.(label_storing)(end+1) = mean(data_labels.(current_label).step_error, 1);
            step_error_var_avg.(label_storing)(end+1) = mean(data_labels.(current_label).step_error_std, 1);

            clear step_error_

            if subtract_all == 1
                for j=1:size(data_labels.(current_label).step_error,1) %To subtract the baseline of each bmh
                    bmh = bmh_list(j);
                    step_error_(j) = abs(data_labels.(current_label).step_error(j) - ((baseline_errors.(bmh).speed_s0a2b0 + baseline_errors.(bmh).speed_s1a2b0 +baseline_errors.(bmh).speed_s2a2b0)/3));
                end
                step_error_avg.(label_storing)(end+1) = mean(step_error_); % Mean across all participants - mean error

            elseif subtract_speeds == 1
                if startsWith(list_labels(i), 's0')
                    for j=1:size(data_labels.(current_label).step_error,1)
                        bmh = bmh_list(j);
                        step_error_(j) = abs(data_labels.(current_label).step_error(j) - baseline_errors.(bmh).speed_s0a2b0);
                    end
                step_error_avg.(label_storing)(end+1) = mean(step_error_); % Mean across all participants - mean error

                elseif startsWith(list_labels(i), 's1')
                    for j=1:size(data_labels.(current_label).step_error,1)
                        bmh = bmh_list(j);
                        step_error_(j) = abs(data_labels.(current_label).step_error(j) - baseline_errors.(bmh).speed_s1a2b0);
                    end
                    step_error_avg.(label_storing)(end+1) = mean(step_error_); % Mean across all participants - mean error

                elseif startsWith(list_labels(i), 's2')
                    for j=1:size(data_labels.(current_label).step_error,1)
                        bmh = bmh_list(j);
                        step_error_(j) = abs(data_labels.(current_label).step_error(j) - baseline_errors.(bmh).speed_s2a2b0);
                    end
                    step_error_avg.(label_storing)(end+1) = mean(step_error_); % Mean across all participants - mean error
                end

            else % No subtractions
                step_error_avg.(label_storing)(end+1) = mean(data_labels.(current_label).step_error, 1);   % Mean across all participants for that label    
            end
        end

        % All values have the same order, as they are the same labels in the same order

        avg_walking_speed_avg.(vars_order(1)) = mean(walking_speed_avg.(vars_order(1)));
        avg_walking_speed_avg.(vars_order(2)) = mean(walking_speed_avg.(vars_order(2)));
        avg_walking_speed_avg.(vars_order(3)) = mean(walking_speed_avg.(vars_order(3)));
        avg_walking_speed_avg.(vars_order(4)) = mean(walking_speed_avg.(vars_order(4)));
        avg_walking_speed_avg.(vars_order(5)) = mean(walking_speed_avg.(vars_order(5)));
        avg_walking_speed_avg.(vars_order(6)) = mean(walking_speed_avg.(vars_order(6)));
        avg_walking_speed_avg.(vars_order(7)) = mean(walking_speed_avg.(vars_order(7)));
        avg_walking_speed_avg.(vars_order(8)) = mean(walking_speed_avg.(vars_order(8)));
        avg_walking_speed_avg.(vars_order(9)) = mean(walking_speed_avg.(vars_order(9)));

        avg_walking_speed_var_avg.(vars_order(1)) = mean(walking_speed_var_avg.(vars_order(1)));
        avg_walking_speed_var_avg.(vars_order(2)) = mean(walking_speed_var_avg.(vars_order(2)));
        avg_walking_speed_var_avg.(vars_order(3)) = mean(walking_speed_var_avg.(vars_order(3)));
        avg_walking_speed_var_avg.(vars_order(4)) = mean(walking_speed_var_avg.(vars_order(4)));
        avg_walking_speed_var_avg.(vars_order(5)) = mean(walking_speed_var_avg.(vars_order(5)));
        avg_walking_speed_var_avg.(vars_order(6)) = mean(walking_speed_var_avg.(vars_order(6)));
        avg_walking_speed_var_avg.(vars_order(7)) = mean(walking_speed_var_avg.(vars_order(7)));
        avg_walking_speed_var_avg.(vars_order(8)) = mean(walking_speed_var_avg.(vars_order(8)));
        avg_walking_speed_var_avg.(vars_order(9)) = mean(walking_speed_var_avg.(vars_order(9)));

        avg_metabolics_avg.(vars_order(1)) = mean(metabolics_avg.(vars_order(1)));
        avg_metabolics_avg.(vars_order(2)) = mean(metabolics_avg.(vars_order(2)));
        avg_metabolics_avg.(vars_order(3)) = mean(metabolics_avg.(vars_order(3)));
        avg_metabolics_avg.(vars_order(4)) = mean(metabolics_avg.(vars_order(4)));
        avg_metabolics_avg.(vars_order(5)) = mean(metabolics_avg.(vars_order(5)));
        avg_metabolics_avg.(vars_order(6)) = mean(metabolics_avg.(vars_order(6)));
        avg_metabolics_avg.(vars_order(7)) = mean(metabolics_avg.(vars_order(7)));
        avg_metabolics_avg.(vars_order(8)) = mean(metabolics_avg.(vars_order(8)));
        avg_metabolics_avg.(vars_order(9)) = mean(metabolics_avg.(vars_order(9)));
        
        avg_step_width_avg.(vars_order(1)) = mean(step_width_avg.(vars_order(1)));
        avg_step_width_avg.(vars_order(2)) = mean(step_width_avg.(vars_order(2)));
        avg_step_width_avg.(vars_order(3)) = mean(step_width_avg.(vars_order(3)));
        avg_step_width_avg.(vars_order(4)) = mean(step_width_avg.(vars_order(4)));
        avg_step_width_avg.(vars_order(5)) = mean(step_width_avg.(vars_order(5)));
        avg_step_width_avg.(vars_order(6)) = mean(step_width_avg.(vars_order(6)));
        avg_step_width_avg.(vars_order(7)) = mean(step_width_avg.(vars_order(7)));
        avg_step_width_avg.(vars_order(8)) = mean(step_width_avg.(vars_order(8)));
        avg_step_width_avg.(vars_order(9)) = mean(step_width_avg.(vars_order(9)));

        avg_step_width_var_avg.(vars_order(1)) = mean(step_width_var_avg.(vars_order(1)));
        avg_step_width_var_avg.(vars_order(2)) = mean(step_width_var_avg.(vars_order(2)));
        avg_step_width_var_avg.(vars_order(3)) = mean(step_width_var_avg.(vars_order(3)));
        avg_step_width_var_avg.(vars_order(4)) = mean(step_width_var_avg.(vars_order(4)));
        avg_step_width_var_avg.(vars_order(5)) = mean(step_width_var_avg.(vars_order(5)));
        avg_step_width_var_avg.(vars_order(6)) = mean(step_width_var_avg.(vars_order(6)));
        avg_step_width_var_avg.(vars_order(7)) = mean(step_width_var_avg.(vars_order(7)));
        avg_step_width_var_avg.(vars_order(8)) = mean(step_width_var_avg.(vars_order(8)));
        avg_step_width_var_avg.(vars_order(9)) = mean(step_width_var_avg.(vars_order(9)));
        
        avg_step_length_avg.(vars_order(1)) = mean(step_length_avg.(vars_order(1)));
        avg_step_length_avg.(vars_order(2)) = mean(step_length_avg.(vars_order(2)));
        avg_step_length_avg.(vars_order(3)) = mean(step_length_avg.(vars_order(3)));
        avg_step_length_avg.(vars_order(4)) = mean(step_length_avg.(vars_order(4)));
        avg_step_length_avg.(vars_order(5)) = mean(step_length_avg.(vars_order(5)));
        avg_step_length_avg.(vars_order(6)) = mean(step_length_avg.(vars_order(6)));
        avg_step_length_avg.(vars_order(7)) = mean(step_length_avg.(vars_order(7)));
        avg_step_length_avg.(vars_order(8)) = mean(step_length_avg.(vars_order(8)));
        avg_step_length_avg.(vars_order(9)) = mean(step_length_avg.(vars_order(9)));

        avg_step_length_var_avg.(vars_order(1)) = mean(step_length_var_avg.(vars_order(1)));
        avg_step_length_var_avg.(vars_order(2)) = mean(step_length_var_avg.(vars_order(2)));
        avg_step_length_var_avg.(vars_order(3)) = mean(step_length_var_avg.(vars_order(3)));
        avg_step_length_var_avg.(vars_order(4)) = mean(step_length_var_avg.(vars_order(4)));
        avg_step_length_var_avg.(vars_order(5)) = mean(step_length_var_avg.(vars_order(5)));
        avg_step_length_var_avg.(vars_order(6)) = mean(step_length_var_avg.(vars_order(6)));
        avg_step_length_var_avg.(vars_order(7)) = mean(step_length_var_avg.(vars_order(7)));
        avg_step_length_var_avg.(vars_order(8)) = mean(step_length_var_avg.(vars_order(8)));
        avg_step_length_var_avg.(vars_order(9)) = mean(step_length_var_avg.(vars_order(9)));
        
        avg_step_error_avg.(vars_order(1)) = mean(step_error_avg.(vars_order(1)));
        avg_step_error_avg.(vars_order(2)) = mean(step_error_avg.(vars_order(2)));
        avg_step_error_avg.(vars_order(3)) = mean(step_error_avg.(vars_order(3)));
        avg_step_error_avg.(vars_order(4)) = mean(step_error_avg.(vars_order(4)));
        avg_step_error_avg.(vars_order(5)) = mean(step_error_avg.(vars_order(5)));
        avg_step_error_avg.(vars_order(6)) = mean(step_error_avg.(vars_order(6)));
        avg_step_error_avg.(vars_order(7)) = mean(step_error_avg.(vars_order(7)));
        avg_step_error_avg.(vars_order(8)) = mean(step_error_avg.(vars_order(8)));
        avg_step_error_avg.(vars_order(9)) = mean(step_error_avg.(vars_order(9)));    

        avg_step_error_var_avg.(vars_order(1)) = mean(step_error_var_avg.(vars_order(1)));
        avg_step_error_var_avg.(vars_order(2)) = mean(step_error_var_avg.(vars_order(2)));
        avg_step_error_var_avg.(vars_order(3)) = mean(step_error_var_avg.(vars_order(3)));
        avg_step_error_var_avg.(vars_order(4)) = mean(step_error_var_avg.(vars_order(4)));
        avg_step_error_var_avg.(vars_order(5)) = mean(step_error_var_avg.(vars_order(5)));
        avg_step_error_var_avg.(vars_order(6)) = mean(step_error_var_avg.(vars_order(6)));
        avg_step_error_var_avg.(vars_order(7)) = mean(step_error_var_avg.(vars_order(7)));
        avg_step_error_var_avg.(vars_order(8)) = mean(step_error_var_avg.(vars_order(8)));
        avg_step_error_var_avg.(vars_order(9)) = mean(step_error_var_avg.(vars_order(9))); 

        avg_stride_time_avg.(vars_order(1)) = mean(stride_time_avg.(vars_order(1)));
        avg_stride_time_avg.(vars_order(2)) = mean(stride_time_avg.(vars_order(2)));
        avg_stride_time_avg.(vars_order(3)) = mean(stride_time_avg.(vars_order(3)));
        avg_stride_time_avg.(vars_order(4)) = mean(stride_time_avg.(vars_order(4)));
        avg_stride_time_avg.(vars_order(5)) = mean(stride_time_avg.(vars_order(5)));
        avg_stride_time_avg.(vars_order(6)) = mean(stride_time_avg.(vars_order(6)));
        avg_stride_time_avg.(vars_order(7)) = mean(stride_time_avg.(vars_order(7)));
        avg_stride_time_avg.(vars_order(8)) = mean(stride_time_avg.(vars_order(8)));
        avg_stride_time_avg.(vars_order(9)) = mean(stride_time_avg.(vars_order(9)));

        avg_stride_time_var_avg.(vars_order(1)) = mean(stride_time_var_avg.(vars_order(1)));
        avg_stride_time_var_avg.(vars_order(2)) = mean(stride_time_var_avg.(vars_order(2)));
        avg_stride_time_var_avg.(vars_order(3)) = mean(stride_time_var_avg.(vars_order(3)));
        avg_stride_time_var_avg.(vars_order(4)) = mean(stride_time_var_avg.(vars_order(4)));
        avg_stride_time_var_avg.(vars_order(5)) = mean(stride_time_var_avg.(vars_order(5)));
        avg_stride_time_var_avg.(vars_order(6)) = mean(stride_time_var_avg.(vars_order(6)));
        avg_stride_time_var_avg.(vars_order(7)) = mean(stride_time_var_avg.(vars_order(7)));
        avg_stride_time_var_avg.(vars_order(8)) = mean(stride_time_var_avg.(vars_order(8)));
        avg_stride_time_var_avg.(vars_order(9)) = mean(stride_time_var_avg.(vars_order(9)));

        avg_head_angle_avg.(vars_order(1)) = mean(head_angle_avg.(vars_order(1)));
        avg_head_angle_avg.(vars_order(2)) = mean(head_angle_avg.(vars_order(2)));
        avg_head_angle_avg.(vars_order(3)) = mean(head_angle_avg.(vars_order(3)));
        avg_head_angle_avg.(vars_order(4)) = mean(head_angle_avg.(vars_order(4)));
        avg_head_angle_avg.(vars_order(5)) = mean(head_angle_avg.(vars_order(5)));
        avg_head_angle_avg.(vars_order(6)) = mean(head_angle_avg.(vars_order(6)));
        avg_head_angle_avg.(vars_order(7)) = mean(head_angle_avg.(vars_order(7)));
        avg_head_angle_avg.(vars_order(8)) = mean(head_angle_avg.(vars_order(8)));
        avg_head_angle_avg.(vars_order(9)) = mean(head_angle_avg.(vars_order(9))); 

        avg_trunk_angle_avg.(vars_order(1)) = mean(trunk_angle_avg.(vars_order(1)));
        avg_trunk_angle_avg.(vars_order(2)) = mean(trunk_angle_avg.(vars_order(2)));
        avg_trunk_angle_avg.(vars_order(3)) = mean(trunk_angle_avg.(vars_order(3)));
        avg_trunk_angle_avg.(vars_order(4)) = mean(trunk_angle_avg.(vars_order(4)));
        avg_trunk_angle_avg.(vars_order(5)) = mean(trunk_angle_avg.(vars_order(5)));
        avg_trunk_angle_avg.(vars_order(6)) = mean(trunk_angle_avg.(vars_order(6)));
        avg_trunk_angle_avg.(vars_order(7)) = mean(trunk_angle_avg.(vars_order(7)));
        avg_trunk_angle_avg.(vars_order(8)) = mean(trunk_angle_avg.(vars_order(8)));
        avg_trunk_angle_avg.(vars_order(9)) = mean(trunk_angle_avg.(vars_order(9))); 
        
        if plot_number == 5
            
            % Means
            scatter(avg_metabolics_avg.(vars_order(1)), avg_step_width_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(avg_metabolics_avg.(vars_order(2)), avg_step_width_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1});% a=0, b=1 -> Square, Green
            scatter(avg_metabolics_avg.(vars_order(3)), avg_step_width_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Triangle, Green
            scatter(avg_metabolics_avg.(vars_order(4)), avg_step_width_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Circle, Red
            scatter(avg_metabolics_avg.(vars_order(5)), avg_step_width_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(avg_metabolics_avg.(vars_order(6)), avg_step_width_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Triangle, Red
            scatter(avg_metabolics_avg.(vars_order(7)), avg_step_width_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Circle, Blue
            scatter(avg_metabolics_avg.(vars_order(8)), avg_step_width_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Square, Blue
            scatter(avg_metabolics_avg.(vars_order(9)), avg_step_width_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
        
            % Individual points (softer)
            var1 = metabolics_avg;
            var2 = step_width_avg;

            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue

            xlabel('Energy Expenditure (W)');
            ylabel('Step Width (cm)');
            title('Energy Expenditure by Step Width');

        elseif plot_number == 6
        
            scatter(avg_metabolics_avg.(vars_order(1)), avg_step_length_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(avg_metabolics_avg.(vars_order(2)), avg_step_length_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(avg_metabolics_avg.(vars_order(3)), avg_step_length_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(avg_metabolics_avg.(vars_order(4)), avg_step_length_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(avg_metabolics_avg.(vars_order(5)), avg_step_length_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(avg_metabolics_avg.(vars_order(6)), avg_step_length_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(avg_metabolics_avg.(vars_order(7)), avg_step_length_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(avg_metabolics_avg.(vars_order(8)), avg_step_length_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(avg_metabolics_avg.(vars_order(9)), avg_step_length_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
        
            % Individual points (softer)
            var1 = metabolics_avg;
            var2 = step_length_avg;

            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue

            xlabel('Energy Expenditure (W)');
            ylabel('Step Length (cm)');
            title('Energy Expenditure by Step Length');

        elseif plot_number == 7
        
            scatter(avg_step_error_avg.(vars_order(1)), avg_metabolics_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(avg_step_error_avg.(vars_order(2)), avg_metabolics_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(avg_step_error_avg.(vars_order(3)), avg_metabolics_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(avg_step_error_avg.(vars_order(4)), avg_metabolics_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(avg_step_error_avg.(vars_order(5)), avg_metabolics_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(avg_step_error_avg.(vars_order(6)), avg_metabolics_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(avg_step_error_avg.(vars_order(7)), avg_metabolics_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(avg_step_error_avg.(vars_order(8)), avg_metabolics_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(avg_step_error_avg.(vars_order(9)), avg_metabolics_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
        
            % Individual points (softer)
            var2 = metabolics_avg;
            var1 = step_error_avg;

            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue

            if normalize_metabolics == 0
                ylabel('Energy Expenditure (W)');
            elseif normalize_metabolics == 1
                ylabel('Energy Expenditure (W/kg)');
            end
            xlabel('Step Error (cm)');
            if subtract_all == 1
                title('Energy Expenditure by Step Error');
            elseif subtract_speeds == 1
                title('Energy Expenditure by Step Error - Mean Error Subtracted for Each Speed');
            else
                title('Energy Expenditure by Step Error');
            end

        elseif plot_number == 8  
        
            scatter(avg_step_width_avg.(vars_order(1)), avg_step_length_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(avg_step_width_avg.(vars_order(2)), avg_step_length_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(avg_step_width_avg.(vars_order(3)), avg_step_length_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(avg_step_width_avg.(vars_order(4)), avg_step_length_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(avg_step_width_avg.(vars_order(5)), avg_step_length_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(avg_step_width_avg.(vars_order(6)), avg_step_length_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(avg_step_width_avg.(vars_order(7)), avg_step_length_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(avg_step_width_avg.(vars_order(8)), avg_step_length_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(avg_step_width_avg.(vars_order(9)), avg_step_length_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue

            % Individual points (softer)
            var1 = step_width_avg;
            var2 = step_length_avg;

            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue

            xlabel('Step Width (cm)');
            ylabel('Step Length (cm)');
            title('Step Width by Step Length');

        elseif plot_number == 9

            scatter(avg_step_error_avg.(vars_order(1)), avg_step_width_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(avg_step_error_avg.(vars_order(2)), avg_step_width_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(avg_step_error_avg.(vars_order(3)), avg_step_width_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(avg_step_error_avg.(vars_order(4)), avg_step_width_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(avg_step_error_avg.(vars_order(5)), avg_step_width_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(avg_step_error_avg.(vars_order(6)), avg_step_width_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(avg_step_error_avg.(vars_order(7)), avg_step_width_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(avg_step_error_avg.(vars_order(8)), avg_step_width_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(avg_step_error_avg.(vars_order(9)), avg_step_width_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
            
            % Individual points (softer)
            var2 = step_width_avg;
            var1 = step_error_avg;

            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue

            ylabel('Step Width (cm)');
            xlabel('Step Error (cm)');
            if subtract_all == 1
                title('Step Width by Step Error - Average Mean Error Subtracted');
            elseif subtract_speeds == 1
                title('Step Width by Step Error - Mean Error Subtracted for Each Speed');
            else
                title('Step Width by Step Error');
            end

        elseif plot_number == 10 

            scatter(avg_step_error_avg.(vars_order(1)), avg_step_length_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(avg_step_error_avg.(vars_order(2)), avg_step_length_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(avg_step_error_avg.(vars_order(3)), avg_step_length_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(avg_step_error_avg.(vars_order(4)), avg_step_length_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(avg_step_error_avg.(vars_order(5)), avg_step_length_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(avg_step_error_avg.(vars_order(6)), avg_step_length_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(avg_step_error_avg.(vars_order(7)), avg_step_length_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(avg_step_error_avg.(vars_order(8)), avg_step_length_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(avg_step_error_avg.(vars_order(9)), avg_step_length_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue

            % Individual points (softer)
            var2 = step_length_avg;
            var1 = step_error_avg;

            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue

            ylabel('Step Length (cm)');
            xlabel('Step Error (cm)');
            if subtract_all == 1
                title('Step Length by Step Error - Average Mean Error Subtracted');
            elseif subtract_speeds == 1
                title('Step Length by Step Error - Mean Error Subtracted for Each Speed');
            else
                title('Step Length by Step Error');
            end

         elseif plot_number == 12 

            scatter(avg_metabolics_avg.(vars_order(1)), avg_step_width_var_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(avg_metabolics_avg.(vars_order(2)), avg_step_width_var_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(avg_metabolics_avg.(vars_order(3)), avg_step_width_var_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(avg_metabolics_avg.(vars_order(4)), avg_step_width_var_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(avg_metabolics_avg.(vars_order(5)), avg_step_width_var_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(avg_metabolics_avg.(vars_order(6)), avg_step_width_var_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(avg_metabolics_avg.(vars_order(7)), avg_step_width_var_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(avg_metabolics_avg.(vars_order(8)), avg_step_width_var_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(avg_metabolics_avg.(vars_order(9)), avg_step_width_var_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue

            % Individual points (softer)
            var1 = metabolics_avg;
            var2 = step_width_var_avg;

            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue

            xlabel('Energy Expenditure (W)');
            ylabel('Step Width Variability (cm)');
            title('Energy Expenditure by Step Width Variability');

        elseif plot_number == 13

            scatter(avg_step_error_avg.(vars_order(1)), avg_step_width_var_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(avg_step_error_avg.(vars_order(2)), avg_step_width_var_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(avg_step_error_avg.(vars_order(3)), avg_step_width_var_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(avg_step_error_avg.(vars_order(4)), avg_step_width_var_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(avg_step_error_avg.(vars_order(5)), avg_step_width_var_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(avg_step_error_avg.(vars_order(6)), avg_step_width_var_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(avg_step_error_avg.(vars_order(7)), avg_step_width_var_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(avg_step_error_avg.(vars_order(8)), avg_step_width_var_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(avg_step_error_avg.(vars_order(9)), avg_step_width_var_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue

            % Individual points (softer)
            var1 = step_error_avg;
            var2 = step_width_var_avg;

            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue

            xlabel('Step Error (cm)');
            ylabel('Step Width Variability (cm)');
            if subtract_all == 1
                title('Step Error by Step Width Variability - Average Mean Error Subtracted');
            elseif subtract_speeds == 1
                title('Step Error by Step Width Variability - Mean Error Subtracted for Each Speed');
            else
                title('Step Error by Step Width Variability');
            end

        elseif plot_number == 14

            scatter(avg_step_error_avg.(vars_order(1)), avg_head_angle_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(avg_step_error_avg.(vars_order(2)), avg_head_angle_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(avg_step_error_avg.(vars_order(3)), avg_head_angle_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(avg_step_error_avg.(vars_order(4)), avg_head_angle_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(avg_step_error_avg.(vars_order(5)), avg_head_angle_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(avg_step_error_avg.(vars_order(6)), avg_head_angle_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(avg_step_error_avg.(vars_order(7)), avg_head_angle_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(avg_step_error_avg.(vars_order(8)), avg_head_angle_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(avg_step_error_avg.(vars_order(9)), avg_head_angle_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue

            % Individual points (softer)
            var1 = step_error_avg;
            var2 = head_angle_avg;

            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue

            xlabel('Step Error (cm)');
            ylabel('Head Angle (deg)');
            if subtract_all == 1
                title('Step Error by Head Angle - Average Mean Error Subtracted');
            elseif subtract_speeds == 1
                title('Step Error by Head Angle - Mean Error Subtracted for Each Speed');
            else
                title('Step Error by Head Angle');
            end

       elseif plot_number == 15 % All points for all participants (circles, and color by accuracy)
            
            hold on
            scatter(data_labels.s0a0b0.step_error, data_labels.s0a0b0.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s0a0b1.step_error, data_labels.s0a0b1.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s0a0b2.step_error, data_labels.s0a0b2.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s0a1b0.step_error, data_labels.s0a1b0.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s0a1b1.step_error, data_labels.s0a1b1.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s0a1b2.step_error, data_labels.s0a1b2.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s0a2b0.step_error, data_labels.s0a2b0.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s0a2b1.step_error, data_labels.s0a2b1.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s0a2b2.step_error, data_labels.s0a2b2.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s1a0b0.step_error, data_labels.s1a0b0.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s1a0b1.step_error, data_labels.s1a0b1.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s1a0b2.step_error, data_labels.s1a0b2.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s1a1b0.step_error, data_labels.s1a1b0.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s1a1b1.step_error, data_labels.s1a1b1.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s1a1b2.step_error, data_labels.s1a1b2.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s1a2b0.step_error, data_labels.s1a2b0.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s1a2b1.step_error, data_labels.s1a2b1.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s1a2b2.step_error, data_labels.s1a2b2.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s2a0b0.step_error, data_labels.s2a0b0.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s2a0b1.step_error, data_labels.s2a0b1.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s2a0b2.step_error, data_labels.s2a0b2.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s2a1b0.step_error, data_labels.s2a1b0.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s2a1b1.step_error, data_labels.s2a1b1.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s2a1b2.step_error, data_labels.s2a1b2.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s2a2b0.step_error, data_labels.s2a2b0.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s2a2b1.step_error, data_labels.s2a2b1.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s2a2b2.step_error, data_labels.s2a2b2.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}, 'HandleVisibility', 'off'); 

            scatter(nan, nan, 200, 'o', 'MarkerEdgeColor', 'g', 'MarkerFaceColor', 'g', 'DisplayName', 'Null Accuracy');
            scatter(nan, nan, 200, 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b', 'DisplayName', 'Medium Accuracy');
            scatter(nan, nan, 200, 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'DisplayName', 'High Accuracy');
            legend('show', 'Location', 'best');

            xlabel('Step Error (cm)');
            ylabel('Head Angle (deg)');
            if subtract_all == 1
                title('Step Error by Head Angle - Average Mean Error Subtracted');
            elseif subtract_speeds == 1
                title('Step Error by Head Angle - Mean Error Subtracted for Each Speed');
            else
                title('Step Error by Head Angle');
            end

            grid on;
            ax = gca;               % Get current axis
            ax.FontSize = 20;       % Set axis number font size (e.g., 12)
            return

        elseif plot_number == 16 % All points for all participants (different shapes, and color by accuracy)
            
            hold on
            scatter(data_labels.s0a0b0.step_error, data_labels.s0a0b0.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s0a0b1.step_error, data_labels.s0a0b1.head_angle, 80, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s0a0b2.step_error, data_labels.s0a0b2.head_angle, 80, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s0a1b0.step_error, data_labels.s0a1b0.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s0a1b1.step_error, data_labels.s0a1b1.head_angle, 80, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s0a1b2.step_error, data_labels.s0a1b2.head_angle, 80, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s0a2b0.step_error, data_labels.s0a2b0.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s0a2b1.step_error, data_labels.s0a2b1.head_angle, 80, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s0a2b2.step_error, data_labels.s0a2b2.head_angle, 80, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s1a0b0.step_error, data_labels.s1a0b0.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s1a0b1.step_error, data_labels.s1a0b1.head_angle, 80, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s1a0b2.step_error, data_labels.s1a0b2.head_angle, 80, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s1a1b0.step_error, data_labels.s1a1b0.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s1a1b1.step_error, data_labels.s1a1b1.head_angle, 80, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s1a1b2.step_error, data_labels.s1a1b2.head_angle, 80, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s1a2b0.step_error, data_labels.s1a2b0.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s1a2b1.step_error, data_labels.s1a2b1.head_angle, 80, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s1a2b2.step_error, data_labels.s1a2b2.head_angle, 80, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s2a0b0.step_error, data_labels.s2a0b0.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s2a0b1.step_error, data_labels.s2a0b1.head_angle, 80, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s2a0b2.step_error, data_labels.s2a0b2.head_angle, 80, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s2a1b0.step_error, data_labels.s2a1b0.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s2a1b1.step_error, data_labels.s2a1b1.head_angle, 80, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s2a1b2.step_error, data_labels.s2a1b2.head_angle, 80, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s2a2b0.step_error, data_labels.s2a2b0.head_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s2a2b1.step_error, data_labels.s2a2b1.head_angle, 80, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s2a2b2.step_error, data_labels.s2a2b2.head_angle, 80, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}, 'HandleVisibility', 'off'); 

            xlabel('Step Error (cm)');
            ylabel('Head Angle (deg)');
            xlim([0 30])
            if subtract_all == 1
                title('Step Error by Head Angle - Average Mean Error Subtracted');
            elseif subtract_speeds == 1
                title('Step Error by Head Angle - Mean Error Subtracted for Each Speed');
            else
                title('Step Error by Head Angle');
            end

        elseif plot_number == 17

            scatter(avg_step_error_avg.(vars_order(1)), avg_trunk_angle_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(avg_step_error_avg.(vars_order(2)), avg_trunk_angle_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(avg_step_error_avg.(vars_order(3)), avg_trunk_angle_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(avg_step_error_avg.(vars_order(4)), avg_trunk_angle_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(avg_step_error_avg.(vars_order(5)), avg_trunk_angle_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(avg_step_error_avg.(vars_order(6)), avg_trunk_angle_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(avg_step_error_avg.(vars_order(7)), avg_trunk_angle_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(avg_step_error_avg.(vars_order(8)), avg_trunk_angle_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(avg_step_error_avg.(vars_order(9)), avg_trunk_angle_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue


            % Individual points (softer)
            var1 = step_error_avg;
            var2 = trunk_angle_avg;

            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue


            xlabel('Step Error (cm)');
            ylabel('Trunk Angle (deg)');
            if subtract_all == 1
                title('Step Error by Trunk Angle - Average Mean Error Subtracted');
            elseif subtract_speeds == 1
                title('Step Error by Trunk Angle - Mean Error Subtracted for Each Speed');
            else
                title('Step Error by Trunk Angle');
            end

        elseif plot_number == 18

            hold on
            scatter(data_labels.s0a0b0.step_error, data_labels.s0a0b0.trunk_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s0a0b1.step_error, data_labels.s0a0b1.trunk_angle, 80, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s0a0b2.step_error, data_labels.s0a0b2.trunk_angle, 80, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s0a1b0.step_error, data_labels.s0a1b0.trunk_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s0a1b1.step_error, data_labels.s0a1b1.trunk_angle, 80, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s0a1b2.step_error, data_labels.s0a1b2.trunk_angle, 80, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s0a2b0.step_error, data_labels.s0a2b0.trunk_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s0a2b1.step_error, data_labels.s0a2b1.trunk_angle, 80, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s0a2b2.step_error, data_labels.s0a2b2.trunk_angle, 80, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s1a0b0.step_error, data_labels.s1a0b0.trunk_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s1a0b1.step_error, data_labels.s1a0b1.trunk_angle, 80, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s1a0b2.step_error, data_labels.s1a0b2.trunk_angle, 80, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s1a1b0.step_error, data_labels.s1a1b0.trunk_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s1a1b1.step_error, data_labels.s1a1b1.trunk_angle, 80, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s1a1b2.step_error, data_labels.s1a1b2.trunk_angle, 80, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s1a2b0.step_error, data_labels.s1a2b0.trunk_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s1a2b1.step_error, data_labels.s1a2b1.trunk_angle, 80, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s1a2b2.step_error, data_labels.s1a2b2.trunk_angle, 80, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s2a0b0.step_error, data_labels.s2a0b0.trunk_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s2a0b1.step_error, data_labels.s2a0b1.trunk_angle, 80, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s2a0b2.step_error, data_labels.s2a0b2.trunk_angle, 80, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s2a1b0.step_error, data_labels.s2a1b0.trunk_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s2a1b1.step_error, data_labels.s2a1b1.trunk_angle, 80, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s2a1b2.step_error, data_labels.s2a1b2.trunk_angle, 80, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s2a2b0.step_error, data_labels.s2a2b0.trunk_angle, 80, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s2a2b1.step_error, data_labels.s2a2b1.trunk_angle, 80, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}, 'HandleVisibility', 'off'); 
            scatter(data_labels.s2a2b2.step_error, data_labels.s2a2b2.trunk_angle, 80, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}, 'HandleVisibility', 'off'); 

            xlabel('Step Error (cm)');
            ylabel('Trunk Angle (deg)');
            xlim([0 30])
            if subtract_all == 1
                title('Step Error by Trunk Angle - Average Mean Error Subtracted');
            elseif subtract_speeds == 1
                title('Step Error by Trunk Angle - Mean Error Subtracted for Each Speed');
            else
                title('Step Error by Trunk Angle');
            end

        elseif plot_number == 20

            var1 = step_width_avg;
            var2 = trunk_angle_avg;
            var1_avg = avg_step_width_avg;
            var2_avg = avg_trunk_angle_avg;

            scatter(var1_avg.(vars_order(1)), var2_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(var1_avg.(vars_order(2)), var2_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(var1_avg.(vars_order(3)), var2_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(var1_avg.(vars_order(4)), var2_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(var1_avg.(vars_order(5)), var2_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(var1_avg.(vars_order(6)), var2_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(var1_avg.(vars_order(7)), var2_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(var1_avg.(vars_order(8)), var2_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(var1_avg.(vars_order(9)), var2_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue

            % Individual points (softer)
            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue

            xlabel('Step Width (cm)');
            ylabel('Trunk Angle (deg)');
            title('Step Width by Trunk Angle');

        elseif plot_number == 21

            var1 = step_width_var_avg;
            var2 = trunk_angle_avg;
            var1_avg = avg_step_width_var_avg;
            var2_avg = avg_trunk_angle_avg;

            scatter(var1_avg.(vars_order(1)), var2_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(var1_avg.(vars_order(2)), var2_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(var1_avg.(vars_order(3)), var2_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(var1_avg.(vars_order(4)), var2_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(var1_avg.(vars_order(5)), var2_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(var1_avg.(vars_order(6)), var2_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(var1_avg.(vars_order(7)), var2_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(var1_avg.(vars_order(8)), var2_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(var1_avg.(vars_order(9)), var2_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue

            % Individual points (softer)
            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue

            xlabel('Step Width Variability (cm)');
            ylabel('Trunk Angle (deg)');
            title('Step Width Variability by Trunk Angle');

        elseif plot_number == 22

            var1 = step_length_var_avg;
            var2 = stride_time_avg;
            var1_avg = avg_step_length_var_avg;
            var2_avg = avg_stride_time_avg;
    
            scatter(var1_avg.(vars_order(1)), var2_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(var1_avg.(vars_order(2)), var2_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(var1_avg.(vars_order(3)), var2_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(var1_avg.(vars_order(4)), var2_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(var1_avg.(vars_order(5)), var2_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(var1_avg.(vars_order(6)), var2_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(var1_avg.(vars_order(7)), var2_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(var1_avg.(vars_order(8)), var2_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(var1_avg.(vars_order(9)), var2_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
    
            % Individual points (softer)
            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue
    
            xlabel('Step Length Variability (cm)');
            ylabel('Stride Duration (s)');
            title('Step Length Variability by Stride Duration');

        elseif plot_number == 23

            var1 = step_length_var_avg;
            var2 = stride_time_var_avg;
            var1_avg = avg_step_length_var_avg;
            var2_avg = avg_stride_time_var_avg;
    
            scatter(var1_avg.(vars_order(1)), var2_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(var1_avg.(vars_order(2)), var2_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(var1_avg.(vars_order(3)), var2_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(var1_avg.(vars_order(4)), var2_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(var1_avg.(vars_order(5)), var2_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(var1_avg.(vars_order(6)), var2_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(var1_avg.(vars_order(7)), var2_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(var1_avg.(vars_order(8)), var2_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(var1_avg.(vars_order(9)), var2_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
    
            % Individual points (softer)
            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue
    
            xlabel('Step Length Variability (cm)');
            ylabel('Stride Duration Variability (s)');
            title('Step Length Variability by Stride Duration Variability');

        elseif plot_number == 24

            var2 = stride_time_var_avg;
            var1 = step_error_avg;
            var2_avg = avg_stride_time_var_avg;
            var1_avg = avg_step_error_avg;
    
            scatter(var1_avg.(vars_order(1)), var2_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(var1_avg.(vars_order(2)), var2_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(var1_avg.(vars_order(3)), var2_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(var1_avg.(vars_order(4)), var2_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(var1_avg.(vars_order(5)), var2_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(var1_avg.(vars_order(6)), var2_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(var1_avg.(vars_order(7)), var2_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(var1_avg.(vars_order(8)), var2_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(var1_avg.(vars_order(9)), var2_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
    
            % Individual points (softer)
            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue
    
            xlabel('Step Error (cm)');
            ylabel('Stride Duration Variability (s)');
            title('Step Error by Stride Duration Variability');

        elseif plot_number == 25
            
            var2 = step_length_avg;
            var1 = step_width_var_avg;
            var2_avg = avg_step_length_avg;
            var1_avg = avg_step_width_var_avg;
    
            scatter(var1_avg.(vars_order(1)), var2_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(var1_avg.(vars_order(2)), var2_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(var1_avg.(vars_order(3)), var2_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(var1_avg.(vars_order(4)), var2_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(var1_avg.(vars_order(5)), var2_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(var1_avg.(vars_order(6)), var2_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(var1_avg.(vars_order(7)), var2_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(var1_avg.(vars_order(8)), var2_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(var1_avg.(vars_order(9)), var2_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
    
            % Individual points (softer)
            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue
    
            xlabel('Step Width Variability (cm)');
            ylabel('Step Length (cm)');
            title('Step Length by Step Width Variability');
                
        elseif plot_number == 26
            
            var2 = step_length_var_avg;
            var1 = step_width_var_avg;
            var2_avg = avg_step_length_var_avg;
            var1_avg = avg_step_width_var_avg;
    
            scatter(var1_avg.(vars_order(1)), var2_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(var1_avg.(vars_order(2)), var2_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(var1_avg.(vars_order(3)), var2_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(var1_avg.(vars_order(4)), var2_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(var1_avg.(vars_order(5)), var2_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(var1_avg.(vars_order(6)), var2_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(var1_avg.(vars_order(7)), var2_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(var1_avg.(vars_order(8)), var2_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(var1_avg.(vars_order(9)), var2_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
    
            % Individual points (softer)
            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue
    
            xlabel('Step Width Variability (cm)');
            ylabel('Stride Length Variability (cm)');
            title('Stride Length Variability by Step Width Variability');

        elseif plot_number == 27
            
            var2 = step_width_avg;
            var1 = step_width_var_avg;
            var2_avg = avg_step_width_avg;
            var1_avg = avg_step_width_var_avg;
    
            scatter(var1_avg.(vars_order(1)), var2_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(var1_avg.(vars_order(2)), var2_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(var1_avg.(vars_order(3)), var2_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(var1_avg.(vars_order(4)), var2_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(var1_avg.(vars_order(5)), var2_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(var1_avg.(vars_order(6)), var2_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(var1_avg.(vars_order(7)), var2_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(var1_avg.(vars_order(8)), var2_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(var1_avg.(vars_order(9)), var2_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
    
            % Individual points (softer)
            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue
    
            xlabel('Step Width Variability (cm)');
            ylabel('Stride Width (cm)');
            title('Stride Width by Step Width Variability');

            elseif plot_number == 28
            
            var2 = metabolics_avg;
            var1 = step_width_var_avg;
            var2_avg = avg_metabolics_avg;
            var1_avg = avg_step_width_var_avg;
    
            scatter(var1_avg.(vars_order(1)), var2_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(var1_avg.(vars_order(2)), var2_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(var1_avg.(vars_order(3)), var2_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(var1_avg.(vars_order(4)), var2_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(var1_avg.(vars_order(5)), var2_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(var1_avg.(vars_order(6)), var2_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(var1_avg.(vars_order(7)), var2_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(var1_avg.(vars_order(8)), var2_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(var1_avg.(vars_order(9)), var2_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
    
            % Individual points (softer)
            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue
    
            xlabel('Step Width Variability (cm)');
            if normalize_metabolics == 0
                ylabel('Energy Expenditure (W)');
            elseif normalize_metabolics == 1
                ylabel('Energy Expenditure (W/kg)');
            end
            title('Energy Expenditure by Step Width Variability');

            elseif plot_number == 29
            
            var2 = metabolics_avg;
            var1 = step_length_avg;
            var2_avg = avg_metabolics_avg;
            var1_avg = avg_step_length_avg;
    
            scatter(var1_avg.(vars_order(1)), var2_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(var1_avg.(vars_order(2)), var2_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(var1_avg.(vars_order(3)), var2_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(var1_avg.(vars_order(4)), var2_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(var1_avg.(vars_order(5)), var2_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(var1_avg.(vars_order(6)), var2_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(var1_avg.(vars_order(7)), var2_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(var1_avg.(vars_order(8)), var2_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(var1_avg.(vars_order(9)), var2_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
    
            % Individual points (softer)
            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue
    
            xlabel('Step Length (cm)');
            if normalize_metabolics == 0
                ylabel('Energy Expenditure (W)');
            elseif normalize_metabolics == 1
                ylabel('Energy Expenditure (W/kg)');
            end
            title('Energy Expenditure by Step Length');

        elseif plot_number == 30
            
            var2 = metabolics_avg;
            var1 = step_length_var_avg;
            var2_avg = avg_metabolics_avg;
            var1_avg = avg_step_length_var_avg;
    
            scatter(var1_avg.(vars_order(1)), var2_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(var1_avg.(vars_order(2)), var2_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(var1_avg.(vars_order(3)), var2_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(var1_avg.(vars_order(4)), var2_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(var1_avg.(vars_order(5)), var2_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(var1_avg.(vars_order(6)), var2_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(var1_avg.(vars_order(7)), var2_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(var1_avg.(vars_order(8)), var2_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(var1_avg.(vars_order(9)), var2_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
    
            % Individual points (softer)
            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue
    
            xlabel('Step Length Variability (cm)');
            if normalize_metabolics == 0
                ylabel('Energy Expenditure (W)');
            elseif normalize_metabolics == 1
                ylabel('Energy Expenditure (W/kg)');
            end
            title('Energy Expenditure by Step Length Variability');

            elseif plot_number == 31
            
            var2 = metabolics_avg;
            var1 = step_width_avg;
            var2_avg = avg_metabolics_avg;
            var1_avg = avg_step_width_avg;
    
            scatter(var1_avg.(vars_order(1)), var2_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(var1_avg.(vars_order(2)), var2_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(var1_avg.(vars_order(3)), var2_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(var1_avg.(vars_order(4)), var2_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(var1_avg.(vars_order(5)), var2_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(var1_avg.(vars_order(6)), var2_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(var1_avg.(vars_order(7)), var2_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(var1_avg.(vars_order(8)), var2_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(var1_avg.(vars_order(9)), var2_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
    
            % Individual points (softer)
            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue
    
            xlabel('Step Width (cm)');
            if normalize_metabolics == 0
                ylabel('Energy Expenditure (W)');
            elseif normalize_metabolics == 1
                ylabel('Energy Expenditure (W/kg)');
            end
            title('Energy Expenditure by Step Width');

        elseif plot_number == 32
            
            var2 = step_length_var_avg;
            var1 = walking_speed_avg;
            var2_avg = avg_step_length_var_avg;
            var1_avg = avg_walking_speed_avg;
    
            scatter(var1_avg.(vars_order(1)), var2_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(var1_avg.(vars_order(2)), var2_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(var1_avg.(vars_order(3)), var2_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(var1_avg.(vars_order(4)), var2_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(var1_avg.(vars_order(5)), var2_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(var1_avg.(vars_order(6)), var2_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(var1_avg.(vars_order(7)), var2_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(var1_avg.(vars_order(8)), var2_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(var1_avg.(vars_order(9)), var2_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
    
            % Individual points (softer)
            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue
    
            xlabel('Walking Speed (m/s)');
            ylabel('Step Length Variability (cm)');
            title('Step Length Variability by Walking Speed');

        elseif plot_number == 33
            
            var2 = metabolics_avg;
            var1 = walking_speed_var_avg;
            var2_avg = avg_metabolics_avg;
            var1_avg = avg_walking_speed_var_avg;
    
            scatter(var1_avg.(vars_order(1)), var2_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(var1_avg.(vars_order(2)), var2_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(var1_avg.(vars_order(3)), var2_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(var1_avg.(vars_order(4)), var2_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(var1_avg.(vars_order(5)), var2_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(var1_avg.(vars_order(6)), var2_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(var1_avg.(vars_order(7)), var2_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(var1_avg.(vars_order(8)), var2_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(var1_avg.(vars_order(9)), var2_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
    
            % Individual points (softer)
            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue
    
            xlabel('Walking Speed Variability (m/s)');
            ylabel('Energy Expenditure (W/kg)');
            title('Energy Expenditure by Walking Speed Variability');

        elseif plot_number == 34
            
            var2 = trunk_angle_avg;
            var1 = walking_speed_var_avg;
            var2_avg = avg_trunk_angle_avg;
            var1_avg = avg_walking_speed_var_avg;
    
            scatter(var1_avg.(vars_order(1)), var2_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(var1_avg.(vars_order(2)), var2_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(var1_avg.(vars_order(3)), var2_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(var1_avg.(vars_order(4)), var2_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(var1_avg.(vars_order(5)), var2_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(var1_avg.(vars_order(6)), var2_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(var1_avg.(vars_order(7)), var2_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(var1_avg.(vars_order(8)), var2_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(var1_avg.(vars_order(9)), var2_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
    
            % Individual points (softer)
            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue
    
            xlabel('Walking Speed Variability (m/s)');
            ylabel('Trunk Angle (deg)');
            title('Trunk Angle by Walking Speed Variability');

        elseif plot_number == 35
            
            var1 = trunk_angle_avg;
            var2 = metabolics_avg;
            var1_avg = avg_trunk_angle_avg;
            var2_avg = avg_metabolics_avg;
    
            scatter(var1_avg.(vars_order(1)), var2_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(var1_avg.(vars_order(2)), var2_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(var1_avg.(vars_order(3)), var2_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(var1_avg.(vars_order(4)), var2_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(var1_avg.(vars_order(5)), var2_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(var1_avg.(vars_order(6)), var2_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(var1_avg.(vars_order(7)), var2_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(var1_avg.(vars_order(8)), var2_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(var1_avg.(vars_order(9)), var2_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
    
            % Individual points (softer)
            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue
    
            xlabel('Trunk Angle (deg)');
            ylabel('Energy Expenditure (W/kg)');
            title('Energy Expenditure by Trunk Angle');

        elseif plot_number == 36
            
            var1 = step_width_avg;
            var2 = step_length_var_avg;
            var1_avg = avg_step_width_avg;
            var2_avg = avg_step_length_var_avg;
    
            scatter(var1_avg.(vars_order(1)), var2_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(var1_avg.(vars_order(2)), var2_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(var1_avg.(vars_order(3)), var2_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(var1_avg.(vars_order(4)), var2_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(var1_avg.(vars_order(5)), var2_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(var1_avg.(vars_order(6)), var2_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(var1_avg.(vars_order(7)), var2_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(var1_avg.(vars_order(8)), var2_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(var1_avg.(vars_order(9)), var2_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
    
            % Individual points (softer)
            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue
    
            xlabel('Step Width (cm)');
            ylabel('Step Length Variability (cm)');
            title('Step Length Variability by Step Width');

        elseif plot_number == 37
            
            var1 = step_width_avg;
            var2 = stride_time_avg;
            var1_avg = avg_step_width_avg;
            var2_avg = avg_stride_time_avg;
    
            scatter(var1_avg.(vars_order(1)), var2_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(var1_avg.(vars_order(2)), var2_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(var1_avg.(vars_order(3)), var2_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(var1_avg.(vars_order(4)), var2_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(var1_avg.(vars_order(5)), var2_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(var1_avg.(vars_order(6)), var2_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(var1_avg.(vars_order(7)), var2_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(var1_avg.(vars_order(8)), var2_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(var1_avg.(vars_order(9)), var2_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
    
            % Individual points (softer)
            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue
    
            xlabel('Step Width (cm)');
            ylabel('Stride Duration (s)');
            title('Stride Duration by Step Width');

        elseif plot_number == 38
            
            var1 = step_width_avg;
            var2 = stride_time_var_avg;
            var1_avg = avg_step_width_avg;
            var2_avg = avg_stride_time_var_avg;
    
            scatter(var1_avg.(vars_order(1)), var2_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(var1_avg.(vars_order(2)), var2_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(var1_avg.(vars_order(3)), var2_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(var1_avg.(vars_order(4)), var2_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(var1_avg.(vars_order(5)), var2_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(var1_avg.(vars_order(6)), var2_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(var1_avg.(vars_order(7)), var2_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(var1_avg.(vars_order(8)), var2_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(var1_avg.(vars_order(9)), var2_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
    
            % Individual points (softer)
            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue
    
            xlabel('Step Width (cm)');
            ylabel('Stride Duration Variability (s)');
            title('Stride Duration Variability by Step Width');

        elseif plot_number == 39
            
            var1 = head_angle_avg;
            var2 = step_length_var_avg;
            var1_avg = avg_head_angle_avg;
            var2_avg = avg_step_length_var_avg;
    
            scatter(var1_avg.(vars_order(1)), var2_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(var1_avg.(vars_order(2)), var2_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(var1_avg.(vars_order(3)), var2_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(var1_avg.(vars_order(4)), var2_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(var1_avg.(vars_order(5)), var2_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(var1_avg.(vars_order(6)), var2_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(var1_avg.(vars_order(7)), var2_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(var1_avg.(vars_order(8)), var2_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(var1_avg.(vars_order(9)), var2_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
    
            % Individual points (softer)
            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue
    
            xlabel('Head Angle (deg)');
            ylabel('Step Length Variability (s)');
            title('Step Length Variability by Head Angle');

        elseif plot_number == 40
            
            var1 = trunk_angle_avg;
            var2 = step_length_var_avg;
            var1_avg = avg_trunk_angle_avg;
            var2_avg = avg_step_length_var_avg;
    
            scatter(var1_avg.(vars_order(1)), var2_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(var1_avg.(vars_order(2)), var2_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(var1_avg.(vars_order(3)), var2_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(var1_avg.(vars_order(4)), var2_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(var1_avg.(vars_order(5)), var2_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(var1_avg.(vars_order(6)), var2_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(var1_avg.(vars_order(7)), var2_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(var1_avg.(vars_order(8)), var2_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(var1_avg.(vars_order(9)), var2_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
    
            % Individual points (softer)
            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue
    
            xlabel('Trunk Angle (deg)');
            ylabel('Step Length Variability (s)');
            title('Step Length Variability by Trunk Angle');

        elseif plot_number == 41
            
            var1 = stride_time_avg;
            var2 = stride_time_var_avg;
            var1_avg = avg_stride_time_avg;
            var2_avg = avg_stride_time_var_avg;
    
            scatter(var1_avg.(vars_order(1)), var2_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(var1_avg.(vars_order(2)), var2_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(var1_avg.(vars_order(3)), var2_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(var1_avg.(vars_order(4)), var2_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(var1_avg.(vars_order(5)), var2_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(var1_avg.(vars_order(6)), var2_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(var1_avg.(vars_order(7)), var2_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(var1_avg.(vars_order(8)), var2_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(var1_avg.(vars_order(9)), var2_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
    
            % Individual points (softer)
            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue
    
            xlabel('Stride Duration (s)');
            ylabel('Stride Duration Variability (s)');
            title('Stride Duration Variability by Stride Duration');

        elseif plot_number == 42
            
            var1 = trunk_angle_avg;
            var2 = head_angle_avg;
            var1_avg = avg_trunk_angle_avg;
            var2_avg = avg_head_angle_avg;
    
            scatter(var1_avg.(vars_order(1)), var2_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(var1_avg.(vars_order(2)), var2_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(var1_avg.(vars_order(3)), var2_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(var1_avg.(vars_order(4)), var2_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(var1_avg.(vars_order(5)), var2_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(var1_avg.(vars_order(6)), var2_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(var1_avg.(vars_order(7)), var2_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(var1_avg.(vars_order(8)), var2_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(var1_avg.(vars_order(9)), var2_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
    
            % Individual points (softer)
            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue
    
            xlabel('Trunk Angle (deg)');
            ylabel('Head Angle (deg)');
            title('Head Angle by Trunk Angle');

        elseif plot_number == 43
            
            var2 = head_angle_avg;
            var1 = step_error_var_avg;
            var2_avg = avg_head_angle_avg;
            var1_avg = avg_step_error_var_avg;
    
            scatter(var1_avg.(vars_order(1)), var2_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(var1_avg.(vars_order(2)), var2_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(var1_avg.(vars_order(3)), var2_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(var1_avg.(vars_order(4)), var2_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(var1_avg.(vars_order(5)), var2_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(var1_avg.(vars_order(6)), var2_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(var1_avg.(vars_order(7)), var2_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(var1_avg.(vars_order(8)), var2_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(var1_avg.(vars_order(9)), var2_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
    
            % Individual points (softer)
            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue
    
            ylabel('Head Angle (deg)');
            xlabel('Step Error Variability (cm)');
            title('Step Error Variability by Head Angle');

        elseif plot_number == 44
            
            var1 = step_error_avg;
            var2 = step_error_var_avg;
            var1_avg = avg_step_error_avg;
            var2_avg = avg_step_error_var_avg;
    
            scatter(var1_avg.(vars_order(1)), var2_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(var1_avg.(vars_order(2)), var2_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(var1_avg.(vars_order(3)), var2_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(var1_avg.(vars_order(4)), var2_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(var1_avg.(vars_order(5)), var2_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(var1_avg.(vars_order(6)), var2_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(var1_avg.(vars_order(7)), var2_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(var1_avg.(vars_order(8)), var2_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(var1_avg.(vars_order(9)), var2_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
    
            % Individual points (softer)
            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue
    
            xlabel('Step Error (cm)');
            ylabel('Step Error Variability (cm)');
            title('Step Error Variability by Step Error');

        elseif plot_number == 45
            
            var1 = step_width_avg;
            var2 = head_angle_avg;
            var1_avg = avg_step_width_avg;
            var2_avg = avg_head_angle_avg;
    
            scatter(var1_avg.(vars_order(1)), var2_avg.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=0 -> Circle, Green
            scatter(var1_avg.(vars_order(2)), var2_avg.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=1 -> Circle, Red
            scatter(var1_avg.(vars_order(3)), var2_avg.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1}); % a=0, b=2 -> Circle, Blue
            scatter(var1_avg.(vars_order(4)), var2_avg.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=0 -> Square, Green
            scatter(var1_avg.(vars_order(5)), var2_avg.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=1 -> Square, Red
            scatter(var1_avg.(vars_order(6)), var2_avg.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2}); % a=1, b=2 -> Square, Blue
            scatter(var1_avg.(vars_order(7)), var2_avg.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=0 -> Triangle, Green
            scatter(var1_avg.(vars_order(8)), var2_avg.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=1 -> Triangle, Red
            scatter(var1_avg.(vars_order(9)), var2_avg.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3}); % a=2, b=2 -> Triangle, Blue
    
            % Individual points (softer)
            scatter(var1.(vars_order(1)), var2.(vars_order(1)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=0 -> Circle, Green
            scatter(var1.(vars_order(2)), var2.(vars_order(2)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1);% a=0, b=1 -> Square, Green
            scatter(var1.(vars_order(3)), var2.(vars_order(3)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{1}, 'MarkerFaceColor', colors{1},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=0, b=2 -> Triangle, Green
            scatter(var1.(vars_order(4)), var2.(vars_order(4)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=0 -> Circle, Red
            scatter(var1.(vars_order(5)), var2.(vars_order(5)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=1 -> Square, Red
            scatter(var1.(vars_order(6)), var2.(vars_order(6)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{2}, 'MarkerFaceColor', colors{2},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=1, b=2 -> Triangle, Red
            scatter(var1.(vars_order(7)), var2.(vars_order(7)), 200, 'Marker', shapes{1}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=0 -> Circle, Blue
            scatter(var1.(vars_order(8)), var2.(vars_order(8)), 200, 'Marker', shapes{2}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=1 -> Square, Blue
            scatter(var1.(vars_order(9)), var2.(vars_order(9)), 200, 'Marker', shapes{3}, 'MarkerEdgeColor', colors{3}, 'MarkerFaceColor', colors{3},'MarkerFaceAlpha', 0.12,'MarkerEdgeAlpha', 0.1); % a=2, b=2 -> Triangle, Blue
    
            xlabel('Step Width (cm)');
            ylabel('Head Angle (deg)');
            title('Head Angle by Step Width');

        end % end of the big loop

        % Create dummy scatter points for legend
        legend_handles = gobjects(9, 1);  % Preallocate graphics objects for legend entries
                
        grid on;
        ax = gca;               % Get current axis
        ax.FontSize = 20;       % Set axis number font size (e.g., 12)

        % Generate each combination for the legend
        if plot_number ~=15
            index = 1;
            for color = colors
                for shape = shapes
                    legend_handles(index) = scatter(NaN, NaN, 200, 'Marker', shape{1}, 'MarkerEdgeColor', color{1}, 'MarkerFaceColor', color{1});
                    index = index + 1;
                end
            end
        else
            return % dont do cluster analysis 
        end

        try % change this
            var1_all = reshape(struct2array(var1), 1, []);
            var2_all = reshape(struct2array(var2), 1, []);
    
            xlim([min(var1_all) - range(var1_all) * 0.04, max(var1_all) + range(var1_all) * 0.04]);  
            ylim([min(var2_all) - range(var2_all) * 0.04, max(var2_all) + range(var2_all) * 0.04]);  
        catch
            return
        end

        if plot_number == 16
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'southwest');

        elseif plot_number == 5
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'northwest');    
        elseif plot_number == 7
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'southwest');
        elseif plot_number == 9
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'northeast');
        elseif plot_number == 10
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'southwest');
        elseif plot_number == 13
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'northwest');
        elseif plot_number == 14
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'southwest');
        elseif plot_number == 17
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'northeast');
        elseif plot_number == 18
            ylim([0 30])
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'northwest');
        elseif plot_number == 19
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'northeast');
        elseif plot_number == 23
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'northwest');
        elseif plot_number == 24
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'northeast');
        elseif plot_number == 26
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'northwest');
        elseif plot_number == 27
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'northwest');
        elseif plot_number == 29
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'northwest');
        elseif plot_number == 31
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'southeast');
        elseif plot_number == 33
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'northwest');
        elseif plot_number == 34
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'northwest');
        elseif plot_number == 35
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'southeast');
        elseif plot_number == 36
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'northwest');
        elseif plot_number == 37
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'northwest');
        elseif plot_number == 38
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'northwest');
        elseif plot_number == 39
            ylim([2 27])
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'northwest');
        elseif plot_number == 40
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'northwest');
        elseif plot_number == 41
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'northwest');
        elseif plot_number == 42
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'southeast');
        elseif plot_number == 43
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'southwest');
        elseif plot_number == 44
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'northwest');
        elseif plot_number == 45
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'southeast');
        else
            legend_handle = legend(legend_handles, legend_labels, 'Location', 'best');
        end
        
        legend_handle.FontSize = 14; % size of the legend

        if show_line == 1
            % Perform linear regression
            coefficients = polyfit(var1_all, var2_all, 1); % Fit line: y = mx + c
            var2_fit = polyval(coefficients, var1_all);   % Predicted y values
            m = coefficients(1); % Slope
            c = coefficients(2); % Intercept
            
            % Compute R^2
            ss_total = sum((var2_all - mean(var2_all)).^2); % Total sum of squares
            ss_residual = sum((var2_all - var2_fit).^2);   % Residual sum of squares
            r_squared = 1 - (ss_residual / ss_total);
            
            % Plot the regression line
            %var1_range = linspace(min(var1_all)-0.01, max(var1_all)+0.01, 100); % Create range for the line
            var1_range = linspace(min(var1_all), max(var1_all), 100);
            var2_range = polyval(coefficients, var1_range);
            plot(var1_range, var2_range, 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
            
            % Display the formula and R^2 in a location away from the line
            if show_formula == 1
                %formula_text = sprintf('y = %.2fx + %.2f\nR^2 = %.2f', m, c, r_squared);
                if c >= 0
                    formula_text = sprintf('y = %.2fx + %.2f\nR^2 = %.2f', m, c, r_squared);
                else
                    formula_text = sprintf('y = %.2fx - %.2f\nR^2 = %.2f', m, abs(c), r_squared);
                end
                if plot_number == 5
                    text(min(var1_all)+0.002, max(var2_all)-4, formula_text, 'FontSize', 16);
                elseif plot_number == 6
                    text(min(var1_all)+0.002, max(var2_all)-4, formula_text, 'FontSize', 16);
                elseif plot_number == 7
                    text(min(var1_all)+0.002, max(var2_all)-4, formula_text, 'FontSize', 16);
                elseif plot_number == 8
                    text(min(var1_all)+0.002, max(var2_all)-4, formula_text, 'FontSize', 16);                   
                elseif plot_number == 9
                    text(min(var1_all), min(var2_all)+0.2, formula_text, 'FontSize', 16);
                elseif plot_number == 10
                    text(min(var1_all)+0.002, max(var2_all)-4, formula_text, 'FontSize', 16);
                elseif plot_number == 12
                    text(min(var1_all)+0.002, max(var2_all)-4, formula_text, 'FontSize', 16);
                elseif plot_number == 13
                    text(min(var1_all)+0.002, max(var2_all)-4, formula_text, 'FontSize', 16);
                elseif plot_number == 14
                    text(max(var1_all)-1.6, max(var2_all)-1, formula_text, 'FontSize', 16);
                elseif plot_number == 15
                    text(min(var1_all)+0.002, max(var2_all)-4, formula_text, 'FontSize', 16);
                elseif plot_number == 16
                    text(min(var1_all)+0.002, max(var2_all)-4, formula_text, 'FontSize', 16);
                elseif plot_number == 17
                    text(min(var1_all)+0.002, max(var2_all)-4, formula_text, 'FontSize', 16);
                elseif plot_number == 18
                    text(min(var1_all)+0.002, max(var2_all)-4, formula_text, 'FontSize', 16);
                elseif plot_number == 19
                    text(min(var1_all)+0.002, max(var2_all)-4, formula_text, 'FontSize', 16);
                elseif plot_number == 20
                    text(min(var1_all)+0.002, max(var2_all)-4, formula_text, 'FontSize', 16);
                elseif plot_number == 21
                    text(min(var1_all)+0.002, max(var2_all)-4, formula_text, 'FontSize', 16);
                elseif plot_number == 22
                    text(min(var1_all)+0.002, max(var2_all)-4, formula_text, 'FontSize', 16);
                elseif plot_number == 23
                    text(max(var1_all)-3, min(var2_all)+0.01, formula_text, 'FontSize', 16);
                elseif plot_number == 24
                    text(min(var1_all)+0.002, max(var2_all)-4, formula_text, 'FontSize', 16);
                elseif plot_number == 25
                    text(min(var1_all)+0.002, max(var2_all)-4, formula_text, 'FontSize', 16);
                elseif plot_number == 26
                    text(min(var1_all)+0.002, max(var2_all)-4, formula_text, 'FontSize', 16);
                elseif plot_number == 27
                    text(min(var1_all)+0.002, max(var2_all)-4, formula_text, 'FontSize', 16);
                elseif plot_number == 28
                    text(min(var1_all)+0.002, max(var2_all)-4, formula_text, 'FontSize', 16);
                elseif plot_number == 29
                    text(min(var1_all)+0.002, max(var2_all)-4, formula_text, 'FontSize', 16);
                elseif plot_number == 30
                    text(min(var1_all)+0.002, max(var2_all)-4, formula_text, 'FontSize', 16);
                elseif plot_number == 31
                    text(min(var1_all)+0.002, max(var2_all)-4, formula_text, 'FontSize', 16);
                elseif plot_number == 33
                    text(max(var1_all)-0.014, min(var2_all)+0.1, formula_text, 'FontSize', 16);   
                elseif plot_number == 34
                    text(max(var1_all)-0.015, min(var2_all)+0.2, formula_text, 'FontSize', 16);    
                elseif plot_number == 35
                    text(min(var1_all)-0.15, max(var2_all)-5, formula_text, 'FontSize', 16);  
                elseif plot_number == 36
                    text(max(var1_all)-1.2, min(var2_all)+0.5, formula_text, 'FontSize', 16);  
                elseif plot_number == 39
                    text(max(var1_all)-6, min(var2_all)+1, formula_text, 'FontSize', 16);
                elseif plot_number == 41
                    text(max(var1_all)-0.15, min(var2_all)+0.02, formula_text, 'FontSize', 16);
                elseif plot_number == 42
                    text(min(var1_all), max(var2_all)-1, formula_text, 'FontSize', 16);
                elseif plot_number == 43
                    text(max(var1_all)-1, max(var2_all)-1, formula_text, 'FontSize', 16);
                elseif plot_number == 44
                    text(max(var1_all)-1.5, min(var2_all)+0.3, formula_text, 'FontSize', 16);
                elseif plot_number == 45
                    text(min(var1_all)-0.1, max(var2_all)-0.5, formula_text, 'FontSize', 16);
                end
            end
        end

        if plot_number == 16 ||  plot_number == 18 % data structure not prepared for this clustering 
            return
        end

        hold off;
    
        %% Cluster analysis
        if show_cluster_analysis == 1
            points_a0b0(:,1) = var1.(vars_order(1));
            points_a0b0(:,2) = var2.(vars_order(1));
            points_a0b1(:,1) = var1.(vars_order(2));
            points_a0b1(:,2) = var2.(vars_order(2));
            points_a0b2(:,1) = var1.(vars_order(3));
            points_a0b2(:,2) = var2.(vars_order(3));
            points_a1b0(:,1) = var1.(vars_order(4));
            points_a1b0(:,2) = var2.(vars_order(4));
            points_a1b1(:,1) = var1.(vars_order(5));
            points_a1b1(:,2) = var2.(vars_order(5));
            points_a1b2(:,1) = var1.(vars_order(6));
            points_a1b2(:,2) = var2.(vars_order(6));
            points_a2b0(:,1) = var1.(vars_order(7));
            points_a2b0(:,2) = var2.(vars_order(7));
            points_a2b1(:,1) = var1.(vars_order(8));
            points_a2b1(:,2) = var2.(vars_order(8));
            points_a2b2(:,1) = var1.(vars_order(9));
            points_a2b2(:,2) = var2.(vars_order(9));
        
            % Silhouette - By accuracy across all conditions (color)

            clusters = [[points_a0b0; points_a0b1; points_a0b2]; 
                        [points_a1b0; points_a1b1; points_a1b2]; 
                        [points_a2b0; points_a2b1; points_a2b2]];
            
            labels_silhouette = [0*ones(9,1); 1*ones(9,1); 2*ones(9,1)];
            silhouetteValues = silhouette(clusters, labels_silhouette);

            uniqueClusters = flip(unique(labels_silhouette)); 
    
            sortedSilhouette = [];
            sortedLabels = [];
            barPositions = [];
            separation = 1.25; % Distance between clusters
            currentY = 1; 
            
            for i = uniqueClusters'
                clusterIdx = find(labels_silhouette == i); % Índices de los puntos del cluster
                silValues = silhouetteValues(clusterIdx); % Valores silhouette del cluster
                
                [sortedSilValues, sortIdx] = sort(silValues, 'descend'); 
                
                sortedSilhouette = [sortedSilhouette; sortedSilValues];
                sortedLabels = [sortedLabels; i * ones(length(sortedSilValues), 1)];
                
                clusterYPositions = (currentY + length(sortedSilValues) - 1):-1:currentY;
                barPositions = [barPositions, clusterYPositions]; 
                
                currentY = currentY + length(sortedSilValues) + separation;
            end
            
            clusterSizes = histcounts(sortedLabels, 'BinMethod', 'integers');
            clusterPositions = cumsum([1, clusterSizes(1:end-1) + separation]) + clusterSizes / 2 - 0.5;
            
            % Horizontal bars
            figure;
            hold on;
            for i = 1:length(sortedSilhouette)
                clusterID = sortedLabels(i);
                color = colors{clusterID + 1}; % Obtener color del cluster
                barh(barPositions(i), sortedSilhouette(i), 'FaceColor', color, 'EdgeColor', 'none');
            end
            hold off;
            
            xlabel('Silhouette Value');
            %ylabel('Cluster');
            title('Silhouette Analysis');
            yticks(clusterPositions);
            yticklabels({clusters_color_silhouette{3}, clusters_color_silhouette{2}, clusters_color_silhouette{1}});
            
            ax = gca;
            ax.FontSize = 20;
            ylim([0 31.3]);
            xlim([-1 1]);
            
            avgSilhouette = mean(silhouetteValues);
            disp(['Average Silhouette Score by accuracy (color): ', num2str(avgSilhouette)]);
            text(0.99, 0.98, sprintf('Avg Silhouette: %.2f', avgSilhouette),'Units', 'normalized', 'FontSize', 16, 'HorizontalAlignment', 'right');
        
            % Silhouette - By balance across all conditions (shape)
            clusters = [[points_a0b0; points_a1b0; points_a2b0]; [points_a0b1; points_a1b1; points_a2b1]; [points_a0b2; points_a1b2; points_a2b2]];
            labels_silhouette = [0*ones(9,1); 1*ones(9,1); 2*ones(9,1)];
            
            figure;
            silhouette(clusters, labels_silhouette);
            
            xlabel('Silhouette Value');
            %ylabel('Cluster');
            title('Silhouette Plot by Balance (Shape)');
            set(gca, 'FontSize', 20);
        
            silhouetteValues = silhouette(clusters, labels_silhouette);
            avgSilhouette = mean(silhouetteValues);
            disp(['Average Silhouette Score by balance (shape): ', num2str(avgSilhouette)]);
    
            avgSilhouette = mean(silhouetteValues);
            text(0.99, 0.965, sprintf('Avg Silhouette: %.2f', avgSilhouette),'Units', 'normalized', 'FontSize', 14, 'HorizontalAlignment', 'right');
            ylim([0 35])
            xlim([-1 1])
        
            % Silhouette - By speed for that condition (shape + color)
            clusters = [points_a0b0; points_a1b0; points_a2b0; points_a0b1; points_a1b1; points_a2b1; points_a0b2; points_a1b2; points_a2b2];
            labels_silhouette = [1*ones(3,1); 2*ones(3,1);3*ones(3,1); 4*ones(3,1); 5*ones(3,1);6*ones(3,1); 7*ones(3,1); 8*ones(3,1); 9*ones(3,1)];
            
            figure;
            silhouette(clusters, labels_silhouette);
            
            xlabel('Silhouette Value');
            %ylabel('Cluster');
            title('Silhouette Plot by Speed for that Condition (Color + Shape)');
            set(gca, 'FontSize', 20);
        
            silhouetteValues = silhouette(clusters, labels_silhouette);
            avgSilhouette = mean(silhouetteValues);
            disp(['Average Silhouette Score by speed for that condition (color + shape): ', num2str(avgSilhouette)]);
    
            avgSilhouette = mean(silhouetteValues);
            text(0.99, 0.965, sprintf('Avg Silhouette: %.2f', avgSilhouette),'Units', 'normalized', 'FontSize', 14, 'HorizontalAlignment', 'right');
            ylim([0 35])
            xlim([-1 1])
        
            % Mean by accuracy (color)
            clusters = [[mean(points_a0b0); mean(points_a0b1); mean(points_a0b2)]; [mean(points_a1b0); mean(points_a1b1); mean(points_a1b2)]; [mean(points_a2b0); mean(points_a2b1); mean(points_a2b2)]];
            labels_silhouette = [0*ones(3,1); 1*ones(3,1); 2*ones(3,1)];
            
            figure;
            silhouette(clusters, labels_silhouette);
            
            xlabel('Silhouette Value');
            %ylabel('Cluster');
            title('Silhouette Plot of the Means by Accuracy (Color)');
            set(gca, 'FontSize', 20);
        
            silhouetteValues = silhouette(clusters, labels_silhouette);
            avgSilhouette = mean(silhouetteValues);
            disp(['Average Silhouette Score of the means by accuracy (color): ', num2str(avgSilhouette)]);
    
            avgSilhouette = mean(silhouetteValues);
            text(0.99, 0.965, sprintf('Avg Silhouette: %.2f', avgSilhouette),'Units', 'normalized', 'FontSize', 14, 'HorizontalAlignment', 'right');
            ylim([1 17])
            xlim([-1 1])
        
            % Mean by shape (balance)
            clusters = [[mean(points_a0b0); mean(points_a1b0); mean(points_a2b0)]; [mean(points_a0b1); mean(points_a1b1); mean(points_a2b1)]; [mean(points_a0b2); mean(points_a1b2); mean(points_a2b2)]];
            labels_silhouette = [0*ones(3,1); 1*ones(3,1); 2*ones(3,1)];
            
            figure;
            silhouette(clusters, labels_silhouette);
            
            xlabel('Silhouette Value');
            %ylabel('Cluster');
            title('Silhouette Plot of the Means by Balance (Shape)');
            set(gca, 'FontSize', 20);
        
            silhouetteValues = silhouette(clusters, labels_silhouette);
            avgSilhouette = mean(silhouetteValues);
            disp(['Average Silhouette Score of the means by balance (shape): ', num2str(avgSilhouette)]);
    
            avgSilhouette = mean(silhouetteValues);
            text(0.99, 0.965, sprintf('Avg Silhouette: %.2f', avgSilhouette),'Units', 'normalized', 'FontSize', 14, 'HorizontalAlignment', 'right');
            ylim([1 17])
            xlim([-1 1])
        end
    end
end