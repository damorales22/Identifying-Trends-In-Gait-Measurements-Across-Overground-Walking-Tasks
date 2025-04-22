function [speed_s0a2b0, speed_s1a2b0, speed_s2a2b0] = BaselineErrors(imported_data,bmh,n,section)

    speed_s0a2b0 = [];
    speed_s1a2b0 = [];
    speed_s2a2b0 = [];
   
    for i = 3:29 % To take only Trials 3-29
        
        trialname = sprintf('Trial%04d', i);

        speed_val = str2double(imported_data.(bmh).prompt.speed(i, 2));
        accuracy_val = str2double(imported_data.(bmh).prompt.accuracy(i, 2));
        balance_val = str2double(imported_data.(bmh).prompt.balance(i, 2));
    
        if speed_val == 0 && accuracy_val == 2 && balance_val == 0
            label_s0a2b0 = trialname;
            if n == 3 % longitudinal - absolute; just straight sections
                if section == 1
                    speed_s0a2b0 = [speed_s0a2b0, imported_data.(bmh).matches.(trialname).overall.overall_mean_absolute];
                elseif section == 2
                    speed_s0a2b0 = [speed_s0a2b0, imported_data.(bmh).matches.(trialname).overall.straight_mean_absolute];
                elseif section == 3
                    speed_s0a2b0 = [speed_s0a2b0, imported_data.(bmh).matches.(trialname).overall.curve_mean_absolute];
                end
            elseif n == 4 % longitudinal - signed; just straight sections
                if section == 1
                    speed_s0a2b0 = [speed_s0a2b0, imported_data.(bmh).matches.(trialname).overall.overall_mean_signed];
                elseif section == 2
                    speed_s0a2b0 = [speed_s0a2b0, imported_data.(bmh).matches.(trialname).overall.straight_mean_signed];
                elseif section == 3
                    speed_s0a2b0 = [speed_s0a2b0, imported_data.(bmh).matches.(trialname).overall.curve_mean_signed];
                end
            end
        elseif speed_val == 1 && accuracy_val == 2 && balance_val == 0
            label_s1a2b0 = trialname;
            if n == 3 % longitudinal - absolute; just straight sections
                if section == 1
                    speed_s1a2b0 = [speed_s1a2b0, imported_data.(bmh).matches.(trialname).overall.overall_mean_absolute];
                elseif section == 2
                    speed_s1a2b0 = [speed_s1a2b0, imported_data.(bmh).matches.(trialname).overall.straight_mean_absolute];
                elseif section == 3
                    speed_s1a2b0 = [speed_s1a2b0, imported_data.(bmh).matches.(trialname).overall.curve_mean_absolute];
                end
            elseif n == 4 % longitudinal - signed; just straight sections
                if section == 1
                    speed_s1a2b0 = [speed_s1a2b0, imported_data.(bmh).matches.(trialname).overall.overall_mean_signed];
                elseif section == 2
                    speed_s1a2b0 = [speed_s1a2b0, imported_data.(bmh).matches.(trialname).overall.straight_mean_signed];
                elseif section == 3
                    speed_s1a2b0 = [speed_s1a2b0, imported_data.(bmh).matches.(trialname).overall.curve_mean_signed];
                end
            end
        elseif speed_val == 2 && accuracy_val == 2 && balance_val == 0
            label_s2a2b0 = trialname;
            if n == 3 % longitudinal - absolute; just straight sections
                if section == 1
                    speed_s2a2b0 = [speed_s2a2b0, imported_data.(bmh).matches.(trialname).overall.overall_mean_absolute];
                elseif section == 2
                    speed_s2a2b0 = [speed_s2a2b0, imported_data.(bmh).matches.(trialname).overall.straight_mean_absolute];
                elseif section == 3
                    speed_s2a2b0 = [speed_s2a2b0, imported_data.(bmh).matches.(trialname).overall.curve_mean_absolute];
                end
            elseif n == 4 % longitudinal - signed; just straight sections
                if section == 1
                    speed_s2a2b0 = [speed_s2a2b0, imported_data.(bmh).matches.(trialname).overall.overall_mean_signed];
                elseif section == 2
                    speed_s2a2b0 = [speed_s2a2b0, imported_data.(bmh).matches.(trialname).overall.straight_mean_signed];
                elseif section == 3
                    speed_s2a2b0 = [speed_s2a2b0, imported_data.(bmh).matches.(trialname).overall.curve_mean_signed];
                end
            end
        end
    end
    
    % Mean across all participants
    speed_s0a2b0 = mean(speed_s0a2b0);
    speed_s1a2b0 = mean(speed_s1a2b0);
    speed_s2a2b0 = mean(speed_s2a2b0);

end