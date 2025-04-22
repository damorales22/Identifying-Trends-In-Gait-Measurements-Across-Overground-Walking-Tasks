function[] = Accuracy_bars_error(imported_data,error,n,section,foot,colors,subtract_speeds_2,subtract_all_2,accuracy_categories,speed_input_2, balance_input_2,speed_s0a2b0,speed_s1a2b0,speed_s2a2b0,bmh)    

% Check if bmh was provided; if not, set it to 0
    if nargin < 15
        bmh_text = " - Across all Participants";
        participants = "multiple";
    else
        participants = "single";
        bmh_text = string(bmh);
    end

%% Detect common features

    common_labels = [];
    datasets = fieldnames(imported_data);

    for i=1:size(datasets)
        bmh = string(datasets{1});
        numtrials = size(imported_data.(bmh).walkingspeed,2); %Number of trials

        for j=1:numtrials %check when every condition is fulfilled
        
            if j==4    %trialname == "Trial0004"
                continue
            end
            if balance_input_2 == str2double(imported_data.(bmh).prompt.balance(j, 2))
                labels_balance = fieldnames(error.bybalance.(imported_data.(bmh).prompt.balance(j, 1)));
            else
                labels_balance = {};
            end
            if speed_input_2 == str2double(imported_data.(bmh).prompt.speed(j, 2))
                labels_speed = fieldnames(error.byspeed.(imported_data.(bmh).prompt.speed(j, 1)));
            else
                labels_speed = {};
            end
        
            if balance_input_2 == 3
                labels_balance = fieldnames(imported_data.(bmh).markerdata);
            end
            if speed_input_2 == 3
                labels_speed = fieldnames(imported_data.(bmh).markerdata);
            end
        
            if (balance_input_2 == str2double(imported_data.(bmh).prompt.balance(j, 2)) && speed_input_2 == str2double(imported_data.(bmh).prompt.speed(j, 2))) || ...
                    (balance_input_2 == str2double(imported_data.(bmh).prompt.balance(j, 2)) && speed_input_2==3) || ...
                    (balance_input_2 == 3 && speed_input_2 == str2double(imported_data.(bmh).prompt.speed(j, 2))) || ...
                    (balance_input_2 == 3 && speed_input_2 == 3)
                
                common_labels = [common_labels; string(sprintf('Trial%04d', j))];
            end   
        end
    end    
    for i = length(common_labels):-1:1  % Iterate in reverse to avoid indexing issues
        trial_name = common_labels{i};
        trial_num = str2double(trial_name(end-1:end)); % Extract the number 'aa' from 'Trial00aa'
        if trial_num < 3 || trial_num > 29 % Remove trials which are not between 3 and 29
            common_labels(i) = [];  % Remove the trial from the cell array
        end
    end

%% Processing
    % Pre-allocate arrays for means and standard deviations
    accuracy_means = zeros(1, length(accuracy_categories)); % Adjust to number of accuracy categories
    accuracy_stds = zeros(1, length(accuracy_categories));
    accuracy_low = [];
    accuracy_medium = [];
    accuracy_high = [];

    for i = 1:length(common_labels)

        trialname = common_labels{i};
        trial_num = str2double(strrep(trialname, 'Trial00', '')); % Extract trial number

        if speed_input_2 == 0
            speed = "Slow";
        elseif speed_input_2 == 1
            speed = "Medium";
        elseif speed_input_2 == 2
            speed = "Fast";
        elseif speed_input_2 == 3
            num = str2double(extractAfter(trialname, "Trial"));
            speed = imported_data.(bmh).prompt.speed(num,1);
        end
    
        if balance_input_2 == 0
            balance = "None";
        elseif balance_input_2 == 1
            balance = "Medium";
        elseif balance_input_2 == 2
            balance = "High";
        elseif balance_input_2 == 3
            num = str2double(extractAfter(trialname, "Trial"));
            balance = imported_data.(bmh).prompt.balance(num,1);
        end

        accuracy = imported_data.(bmh).prompt.accuracy(trial_num,1);
        
        if section ==1
            accuracy_trial_r = imported_data.(bmh).matches.(trialname).right.match;
            accuracy_trial_l = imported_data.(bmh).matches.(trialname).left.match;
        elseif section == 2
            accuracy_trial_r = imported_data.(bmh).matches.(trialname).right.match_straight;
            accuracy_trial_l = imported_data.(bmh).matches.(trialname).left.match_straight;
        elseif section == 3
            accuracy_trial_r = imported_data.(bmh).matches.(trialname).right.match_curve;
            accuracy_trial_l = imported_data.(bmh).matches.(trialname).left.match_curve;
        end

        accuracy_trial_r = accuracy_trial_r(:,n);
        accuracy_trial_l = accuracy_trial_l(:,n);

        if foot==1
            accuracy_trial = [accuracy_trial_r; accuracy_trial_l];
        elseif foot==2
            accuracy_trial = accuracy_trial_r;
        elseif foot==3
            accuracy_trial = accuracy_trial_l;
        end

         %Separate depending on the accuracy
        if accuracy == "Low"
            accuracy_low = [accuracy_low; accuracy_trial];
        elseif accuracy == "Medium"
            accuracy_medium = [accuracy_medium; accuracy_trial];
        elseif accuracy == "High"
            accuracy_high = [accuracy_high; accuracy_trial];
        end
    end

%% Mean substracction

    subtract_text = "";
    accuracy_means = [mean(accuracy_low), mean(accuracy_medium), mean(accuracy_high)]/10; %to convert it to cm
    
    if subtract_all_2 == 1 && (n == 3 || n == 4)
        accuracy_means = accuracy_means - ((speed_s0a2b0+speed_s1a2b0+speed_s2a2b0)/30); %3 to calculate the mean * 10 for cm
        subtract_text = " - Mean Subtracted";     
    elseif subtract_speeds_2 == 1 && (n == 3 || n == 4)
        accuracy_means(1) = accuracy_means(1) - speed_s0a2b0/10;
        accuracy_means(2) = accuracy_means(2) - speed_s1a2b0/10;
        accuracy_means(3) = accuracy_means(3) - speed_s2a2b0/10;
        subtract_text = " - Mean Subtracted According to Speed";   
    end
    
    accuracy_stds = [std(accuracy_low), std(accuracy_medium), std(accuracy_high)]/10;     %to convert it to cm

%% Plot
    % Create the bar graph with error bars for accuracy
    figure;
    b_acc = bar(accuracy_means); % Bar graph of means;
    b_acc.FaceColor = 'flat'; % Enable color customization
    for k = 1:length(accuracy_categories)
        b_acc.CData(k, :) = colors(k, :); % Assign a different color to each bar
    end
    hold on;
    errorbar(1:length(accuracy_categories), accuracy_means, accuracy_stds, 'k', 'linestyle', 'none'); % Add error bars
    
    % Customize plot
    set(gca, 'XTick', 1:length(accuracy_categories), 'XTickLabel', accuracy_categories, 'FontSize', 15);
    ylabel('Step Error (cm)', 'FontSize', 15);
    xlabel('Accuracy Condition', 'FontSize', 15);
    
    if n==3
        n_text = "Longitudinal Error";
    elseif n==4
        n_text = "Longitudinal Error (Absolute Value)";
    elseif n==5
        n_text = "Total Error";
    elseif n==6
        n_text = "Lateral Error";
    end

    if section==1
        section_text = "All the Circuit";
    elseif section==2
        section_text = "Straight Sections";
    elseif section==3
        section_text = "Curved Sections";
    end

    if foot==1
        foot_text = "Both Feet";
    elseif foot==2
        foot_text = "Right Foot";
    elseif foot==3
        foot_text = "Left Foot";
    end
    
    if speed_input_2==3
        speed_text = 'All';
    elseif speed_input_2==0
        speed_text = "Slow";
    elseif speed_input_2==1
        speed_text = "Medium";
    elseif speed_input_2==2
        speed_text = "Fast";
    end

    if balance_input_2==3
        balance_text = 'All';
    elseif balance_input_2==0
        balance_text = "None";
    elseif balance_input_2==1
        balance_text = "Medium";
    elseif balance_input_2==2
        balance_text = "High";
    end

    title({'Step Error Grouped by Accuracy Condition' + bmh_text + ' - ' + n_text + ' - ' + section_text; foot_text + ' - Speed: ' + speed_text + ' - Balance: ' + balance_text + subtract_text}, 'FontSize', 20);

end %end of the function
