%% New Method for Calculating Step Width 
function [stepparameters] = StepParametersHuxham(markerdata, steps, pos_limit, neg_limit)    
    %% Get Step Locations + Times
    left = zeros(length(steps.left.frames(:,1)),4);
    right = zeros(length(steps.right.frames(:,1)),4);

    % left
    for i = 1:length(left(:,1))
        timestamp = steps.left.frames(i);
        step = markerdata.LHEEL(:, 1:3, timestamp); %just takes heel marker
        left(i,1) = timestamp; %when the foot is horizontally on the ground
        left(i, 2:4) = step;  %x,y,z
    end 
    % right
    for i = 1:length(right(:,1))
        timestamp = steps.right.frames(i);
        step = markerdata.RHEEL(:, 1:3, timestamp);
        right(i,1) = timestamp;
        right(i, 2:4) = step;
    end
    
    % Declares empty array for steps from both feet in successive (cronological) order
    combined_steps = NaN(length(left(:,1)) + length(right(:,1)), 5);
    
    % Variables to keep track of indices/skipped steps during running
    ind_r = 1;
    ind_l = 1;
    problem = 0;
    
    % Combines steps based on time stamp
    for i=1:length(combined_steps(:,1))-2
        if ind_r < length(right(:,1)) && right(ind_r,1) <= left(ind_l,1)
            combined_steps(i,1:4) = right(ind_r,:);
            combined_steps(i,5) = 0; % Right foot
            ind_r = ind_r + 1;
        elseif ind_l < length(left(:,1)) && left(ind_l,1) < right(ind_r,1)
            combined_steps(i,1:4) = left(ind_l, :);
            combined_steps(i,5) = 1; % Left foot
            ind_l = ind_l + 1;
        else
            problem = problem + 1; % Increments if any steps are skipped (two consecutive right or left)
        end
    end

    %disp(["Number of combined_steps failed: " problem])

    combined_steps(:,4) = 0; % makes the z's zero
    combined_steps = combined_steps(~isnan(combined_steps(:, 1)), :); % remove unused rows

    %% Calculate Step Width on Straights
    % Declares empty arrays to save steps
    stepwidthstraights = NaN(length(combined_steps(:,1))-2,1);
    steplengthstraights = NaN(length(combined_steps(:,1))-2,1);
    stepwidthcurves = NaN(length(combined_steps(:,1))-2,1);
    steplengthcurves = NaN(length(combined_steps(:,1))-2,1);
    
    % Counts how many skipped steps for successive foots of same foot
    samefoot = 0;
    % Counts how many negated crossover steps
    negated = 0;
    crossovers = [];
    % coordinates for point (to check crossover step)
    crossover_point = [-1650 2150]; % x, y of designated center point based on current track

    for i=1:length(stepwidthstraights (:,1))
        % Checks whether successive steps are from same foot
        if combined_steps(i,5) == combined_steps(i+1,5) || combined_steps(i+1,5) == combined_steps(i+2,5)
            samefoot = samefoot + 1; %If two consecutive steps are from the same foot, increment samefoot and skip 
            continue
        end

        % Sets variables for parameters calculations
        oppstep = norm(combined_steps(i,2:4)-combined_steps(i+1,2:4));
        step = norm(combined_steps(i+1,2:4)-combined_steps(i+2,2:4));
        stride = norm(combined_steps(i,2:4)-combined_steps(i+2,2:4));
        
        if (combined_steps(i+2,3) < pos_limit && combined_steps(i,3) > neg_limit && combined_steps(i,2) > -1300) ... % Straights
                || combined_steps(i,3) < pos_limit && combined_steps(i+2,3) > neg_limit && combined_steps(i,2) < -1300
            
            % Steps on straights
            steplengthstraights(i) = (step^2 + stride^2 - oppstep^2)/(2*stride);
            stepwidthstraights(i) = sqrt(step^2 - steplengthstraights (i)^2);

            % Check whether crossover step
            fststepdist = abs(combined_steps(i,2) - crossover_point(1));
            sndstepdist = abs(combined_steps(i+1,2) - crossover_point(1));
            if combined_steps(i,5) == 0 && fststepdist < sndstepdist ...
                    || combined_steps(i,5) == 1 && fststepdist > sndstepdist
                % Negate crossover steps
                stepwidthstraights(i) = -1*stepwidthstraights(i);
                crossovers(end+1) = combined_steps(i,1);
                negated = negated + 1;
            end

        else 
            % Steps on curves
            steplengthcurves(i) = (step^2 + stride^2 - oppstep^2)/(2*stride); %These formulas because we are in a 3D environment; stride: dist between strikes of the same foot
            stepwidthcurves(i) = sqrt(step^2 - steplengthcurves (i)^2);

            % Check whether crossover step
            fststepdist = norm(combined_steps(i,2:3) - crossover_point);
            sndstepdist = norm(combined_steps(i,2:3) - crossover_point);
            if combined_steps(i,5) == 0 && fststepdist < sndstepdist ...
                    || combined_steps(i,5) == 1 && fststepdist > sndstepdist
                % Negate crossover steps
                stepwidthcurves(i) = -1*stepwidthcurves(i);
                crossovers(end+1) = combined_steps(i,1);
                negated = negated + 1;
            end
        end
    end

    % Display how many removed because of consecutive steps with the same foot
    %disp(["Number removed bc same foot: " samefoot]); 
    %disp(["Number of crossover steps: " negated]); %step in which one foot is crossed over another
    %disp(["Crossover steps: " crossovers]);    
    
    % Remove unused rows
    stepwidthstraights = stepwidthstraights(~isnan(stepwidthstraights(:, 1)), :);
    steplengthstraights = steplengthstraights(~isnan(steplengthstraights (:, 1)), :);

    stepwidthcurves = stepwidthcurves(~isnan(stepwidthcurves(:, 1)), :);
    steplengthcurves = steplengthcurves(~isnan(steplengthcurves (:, 1)), :);
    
    % Creates struct for output, saves arrays + calculates mean/variability; there are straight and curve sections in thebcircuit
    stepparameters = struct();

    stepparameters.stepwidth.step_width_straights = stepwidthstraights;
    stepparameters.stepwidth.straights_variability = std(stepwidthstraights);
    stepparameters.stepwidth.straights_mean = mean(stepwidthstraights);
    stepparameters.stepwidth.straights_median = median(stepwidthstraights);

    stepparameters.stepwidth.step_width_curves = stepwidthcurves;
    stepparameters.stepwidth.curves_variability = std(stepwidthcurves);
    stepparameters.stepwidth.curves_mean = mean(stepwidthcurves);

    stepparameters.steplength.step_length_straights = steplengthstraights;
    stepparameters.steplength.straights_variability = std(steplengthstraights);
    stepparameters.steplength.straights_mean = mean(steplengthstraights);
    stepparameters.steplength.straights_median = median(steplengthstraights);

    stepparameters.steplength.step_length_curves = steplengthcurves;
    stepparameters.steplength.curves_variability = std(steplengthcurves);
    stepparameters.steplength.curves_mean = mean(steplengthcurves);
end