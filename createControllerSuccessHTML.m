function successMsg = createControllerSuccessHTML(K, num, den)
% CREATECONTROLLERSUCCESSHTML - Generate HTML success message for controller design
%
% This function creates an HTML success message with details about the
% designed controller and transfer function.
%
% Parameters:
%   K - The designed state feedback gain matrix
%   num - Numerator coefficients of the transfer function
%   den - Denominator coefficients of the transfer function
%
% Returns:
%   successMsg - HTML formatted success message

    % Create success message with controller info
    successMsg = '<html><div style="color: green; text-align: center; margin-top: 10px; font-weight: bold; padding: 10px; border: 1px solid green; border-radius: 5px;">';
    successMsg = [successMsg, 'State-space controller successfully designed!<br><br>'];
    
    % Add controller gain information
    successMsg = [successMsg, 'Controller gain K = ['];
    for i = 1:size(K, 1)
        if i > 1
            successMsg = [successMsg, '; '];
        end
        for j = 1:size(K, 2)
            if j > 1
                successMsg = [successMsg, ', '];
            end
            successMsg = [successMsg, num2str(K(i,j), '%.4g')];
        end
    end
    successMsg = [successMsg, ']<br><br>'];
    
    % Add transfer function information
    successMsg = [successMsg, 'Closed-loop transfer function G(s) created<br>'];
    successMsg = [successMsg, 'Numerator order: ', num2str(length(num)-1), '<br>'];
    successMsg = [successMsg, 'Denominator order: ', num2str(length(den)-1), '<br><br>'];
    
    % Add next steps
    successMsg = [successMsg, 'The controller is incorporated into G(s), and K(s) is set to 1.<br>'];
    successMsg = [successMsg, 'You can now click "Analyze" to analyze the closed-loop system.'];
    successMsg = [successMsg, '</div></html>'];
end