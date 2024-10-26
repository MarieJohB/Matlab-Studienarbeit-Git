function valuesVector = generateFunctionValues(r, t)
    % Validate inputs
    if ~isa(r, 'function_handle')
        error('The first input must be a function handle.');
    end
    if ~isnumeric(t) || ~isvector(t)
        error('The second input must be a numeric vector.');
    end
    
    % Generate values of r at each step in t
    valuesVector = arrayfun(r, t);
end