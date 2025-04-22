function [dist_gaze] = Dist_gaze(trial, markerdata, trialname, static, head_angles, trunk_angles, height, bmh)

    pos_C7 = find(strcmp(static.Labels, 'C7'));
    coords_C7 = permute(static.Data(pos_C7,1:3,:), [3, 2, 1]);
    
    pos_LPSIS = find(strcmp(static.Labels, 'LPSIS'));
    coords_LPSIS = permute(static.Data(pos_LPSIS,1:3,:), [3, 2, 1]);

    pos_RPSIS = find(strcmp(static.Labels, 'RPSIS'));
    coords_RPSIS = permute(static.Data(pos_RPSIS,1:3,:), [3, 2, 1]);

    dist = (coords_C7(100,3) - ((coords_LPSIS(100,3) + coords_RPSIS(100,3))/2))/10; % to cm
    
    head_angles = mean(head_angles);
    head_angles = head_angles(~any(isnan(head_angles), 2), :); % Filter out rows with NaN

    d = abs(dist * cos(deg2rad(trunk_angles)));
    %r = abs(d * sin(trunk_angles));
    r = abs(dist * sin(deg2rad(trunk_angles)));

    delta_h = dist-d;

    b = abs((height - delta_h) / tan(deg2rad(head_angles)));
    dist_gaze = b - r;
end