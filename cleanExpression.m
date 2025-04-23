function cleaned_str = cleanExpression(str)
    % Clean up expressions:
    % 1. Remove "+ 0*t" or "+ 0 * t" patterns
    % 2. Remove "- 0*t" or "- 0 * t" patterns
    % 3. Clean up multiplication signs with t for better display

    % Remove "+ 0*t" patterns with various spacing
    str = regexprep(str, '\+\s*0\s*\*\s*t', '');
    
    % Remove "- 0*t" patterns with various spacing
    str = regexprep(str, '\-\s*0\s*\*\s*t', '');
    
    % Improve formatting for multiplication with t (x*t → x·t)
    str = regexprep(str, '(\d+)\s*\*\s*t', '$1 · t');
    
    % Trim any leading/trailing whitespace
    cleaned_str = strtrim(str);
    
    % If the result is empty, return a '0'
    if isempty(cleaned_str)
        cleaned_str = '0';
    end
end