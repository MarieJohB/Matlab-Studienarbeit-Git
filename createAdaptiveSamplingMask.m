% Helper function to create adaptive sampling mask for efficient plotting
function mask = createAdaptiveSamplingMask(y, t, target_points)
    % Create a mask that preferentially keeps points in regions of high signal change
    
    % Start with a base mask for important points
    mask = false(size(t));
    mask(1) = true;  % First point
    mask(end) = true;  % Last point
    
    % Calculate signal changes
    dy = abs(diff(y));
    
    % Normalize changes to 0-1 range
    if max(dy) > 0
        norm_dy = dy / max(dy);
    else
        norm_dy = zeros(size(dy));
    end
    
    % Find regions with significant changes (transient regions)
    significant_changes = norm_dy > 0.01;
    
    % Always include points with significant changes
    for i = 1:length(significant_changes)
        if significant_changes(i)
            mask(i:i+1) = true;  % Keep point and next point
        end
    end
    
    % Count how many points we have so far
    current_points = sum(mask);
    
    % If we need more points to reach target, add them adaptively
    if current_points < target_points
        % Sort remaining candidates by normalized change
        candidates = find(~mask);
        if ~isempty(candidates) && length(candidates) > 1
            % For first point, use the first difference
            if candidates(1) == 1
                candidate_changes = [norm_dy(1); norm_dy(candidates(2:end)-1)];
            else
                candidate_changes = norm_dy(candidates-1);
            end
            
            % Sort candidates by their change values
            [~, sorted_idx] = sort(candidate_changes, 'descend');
            sorted_candidates = candidates(sorted_idx);
            
            % Add points up to target or until we run out
            points_to_add = min(target_points - current_points, length(sorted_candidates));
            mask(sorted_candidates(1:points_to_add)) = true;
        end
    end
    
    % If we still need more points, add evenly spaced ones
    current_points = sum(mask);
    if current_points < target_points
        remaining_candidates = find(~mask);
        spacing = max(1, floor(length(remaining_candidates) / (target_points - current_points)));
        mask(remaining_candidates(1:spacing:end)) = true;
    end
end