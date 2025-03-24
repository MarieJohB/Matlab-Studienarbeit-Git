function value = generateRandomParameter(paramSetting)
    % GENERATERANDOMPARAMETER - Generate random parameter values based on distribution
    
    switch paramSetting.dist
        case 'Normal'
            % Generate normal distribution with mean param1 and std dev param2
            value = paramSetting.param1 + paramSetting.param2 * randn();
            
        case 'Uniform'
            % Generate uniform distribution between param1 and param2
            value = paramSetting.param1 + (paramSetting.param2 - paramSetting.param1) * rand();
            
        case 'Fixed'
            % Return fixed value param1
            value = paramSetting.param1;
            
        otherwise
            % Default to fixed value
            value = paramSetting.param1;
    end
end