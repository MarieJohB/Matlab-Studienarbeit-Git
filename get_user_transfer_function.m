function G = get_user_transfer_function()
    % Create input dialog for numerator and denominator coefficients
    prompt = {'Enter the numerator coefficients as a vector [b0 b1 b2 ...]:', ...
              'Enter the denominator coefficients as a vector [a0 a1 a2 ...]:'};
    dlgtitle = 'Transfer Function Input';
    dims = [1 50];

    % Display input dialog
    answer = inputdlg(prompt, dlgtitle, dims);

    % Validate input and create the transfer function
    try
        num = str2num(answer{1}); %#ok<ST2NM>
        den = str2num(answer{2}); %#ok<ST2NM>
        G = tf(num, den);
    catch
        uiwait(msgbox('Invalid input. Please enter the coefficients correctly.', 'Error','error'));
        G = [];
    end
    
    % Display the transfer function
    if ~isempty(G)
        disp('The entered transfer function is:');
        disp(G);
    end
end