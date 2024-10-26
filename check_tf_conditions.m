function [isStable, isProper, isSternlyProper] = check_tf_conditions(tf_sys)
    % Check if the transfer function is sternly proper
    if sternlyproper_tf(tf_sys)
        isSternlyProper = '\surd';
        isProper = '\surd';
    elseif proper_tf(tf_sys)
        isSternlyProper = '\times';
        isProper = '\surd';
    else
        isSternlyProper = '\times';
        isProper = '\times';
    end
    
    % Check stability using negative real parts of roots of the denominator
    if checkRootsNegativeReal(tf_sys)
        isStable = '\surd';
    else
        isStable = '\times';
    end
end