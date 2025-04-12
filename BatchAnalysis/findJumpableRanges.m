function ranges = findJumpableRanges(paramValues, isJumpable)
% Find continuous ranges of parameter values where system is jumpable

% Initialize
ranges = [];

% No jumpable values
if ~any(isJumpable)
    return;
end

% Find transitions between jumpable and non-jumpable
transitions = find(diff([0 isJumpable 0]));

% Convert to ranges
for i = 1:2:length(transitions)
    if i+1 <= length(transitions)
        startIdx = transitions(i);
        endIdx = transitions(i+1) - 1;
        ranges = [ranges; paramValues(startIdx) paramValues(endIdx)];
    end
end
end