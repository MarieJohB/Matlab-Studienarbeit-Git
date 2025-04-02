% Local function for creating success message with controller details
function htmlMsg = createSuccessHTML(K, num, den)
    htmlMsg = '<html><div style="color: green; text-align: center; margin-top: 10px; font-weight: bold; padding: 10px; border: 1px solid green; border-radius: 5px;">';
    htmlMsg = [htmlMsg, 'Controller designed and state-space model converted to transfer function G(s)!<br><br>'];
    
    % Add model information
    htmlMsg = [htmlMsg, 'Continuous-time transfer function<br>'];
    
    htmlMsg = [htmlMsg, 'Numerator order: ', num2str(length(num)-1), '<br>'];
    htmlMsg = [htmlMsg, 'Denominator order: ', num2str(length(den)-1), '<br><br>'];
    
    % Print coefficient information
    htmlMsg = [htmlMsg, 'Numerator coefficients: ['];
    for i = 1:length(num)
        if i > 1
            htmlMsg = [htmlMsg, ', '];
        end
        htmlMsg = [htmlMsg, num2str(num(i), '%.4g')];
    end
    htmlMsg = [htmlMsg, ']<br>'];
    
    htmlMsg = [htmlMsg, 'Denominator coefficients: ['];
    for i = 1:length(den)
        if i > 1
            htmlMsg = [htmlMsg, ', '];
        end
        htmlMsg = [htmlMsg, num2str(den(i), '%.4g')];
    end
    htmlMsg = [htmlMsg, ']<br><br>'];
    
    % Add controller gain matrix
    htmlMsg = [htmlMsg, 'Controller Gain Matrix K:<br>'];
    for i = 1:size(K, 1)
        htmlMsg = [htmlMsg, '['];
        for j = 1:size(K, 2)
            if j > 1
                htmlMsg = [htmlMsg, ', '];
            end
            htmlMsg = [htmlMsg, num2str(K(i,j), '%.4g')];
        end
        htmlMsg = [htmlMsg, ']<br>'];
    end
    
    htmlMsg = [htmlMsg, '<br>The controller implements state feedback u = -Kx.<br>'];
    htmlMsg = [htmlMsg, 'The transfer function G(s) represents the closed-loop system with this controller.<br><br>'];
    
    htmlMsg = [htmlMsg, 'You can now analyze the control system.'];
    htmlMsg = [htmlMsg, '</div></html>'];
end