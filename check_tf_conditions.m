function [isStable, isProper, isStrictlyProper] = check_tf_conditions(tf_sys)
    % CHECK_TF_CONDITIONS Check if a transfer function satisfies various conditions
    %
    % This function checks if the transfer function system tf_sys:
    % - is stable (all poles have negative real parts, verified using Hurwitz criterion),
    % - is proper (numerator degree ≤ denominator degree), and
    % - is strictly proper (numerator degree < denominator degree).
    %
    % Input:
    %   tf_sys - A transfer function object
    %
    % Output:
    %   isStable - Boolean indicating stability (true if stable)
    %   isProper - Boolean indicating if transfer function is proper
    %   isStrictlyProper - Boolean indicating if transfer function is strictly proper
    
    % Extract the coefficients of the numerator and denominator polynomials
    [num, den] = tfdata(tf_sys, 'v');
    
    % Remove leading zeros for correct degree calculation
    num = removeLeadingZeros(num);
    den = removeLeadingZeros(den);
    
    % Determine the degrees (number of coefficients minus 1)
    degNum = length(num) - 1;
    degDen = length(den) - 1;
    
    % Check if proper: numerator degree ≤ denominator degree
    isProper = (degNum <= degDen);
    
    % Check if strictly proper: numerator degree < denominator degree
    isStrictlyProper = (degNum < degDen);
    
    % Check stability using Hurwitz criterion
    isStable = check_stability_hurwitz(tf_sys);
end