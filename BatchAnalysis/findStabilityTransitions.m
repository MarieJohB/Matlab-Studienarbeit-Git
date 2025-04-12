function transitions = findStabilityTransitions(batchResults)
% Find indices where stability changes
stability = batchResults.stability;
transitions = find(diff(stability) ~= 0);
end