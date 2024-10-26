function G = get_user_transfer_function()
    % Create input dialog for numerator and denominator coefficients
    prompt = {'Enter the numerator coefficients as a vector [b0 b1 b2 ...]:', ...
              'Enter the denominator coefficients as a vector [a0 a1 a2 ...]:'};
    dlgtitle = 'Transfer Function Input';
    dims = [1 50];

    % Initialize G to be empty
    G = [];

    % Loop until valid coefficients are entered or user cancels the dialog
    while isempty(G)
        % Display input dialog
        answer = inputdlg(prompt, dlgtitle, dims);

        % Check if user cancelled the dialog
        if isempty(answer)
            disp('Operation cancelled by user.');
            return;
        end

        % Validate input and create the transfer function
        try
            % Replace commas with periods
            num_str = strrep(answer{1}, ',', '.');
            den_str = strrep(answer{2}, ',', '.');
            % Split by space to handle zeroes correctly
            num = str2double(strsplit(num_str, ' '));
            den = str2double(strsplit(den_str, ' '));

            if any(isnan(num)) || any(isnan(den))
                error('Invalid input. Please enter non-empty coefficient vectors.');
            end

            G = tf(num, den);
        catch
            uiwait(msgbox('Invalid input. Please enter the coefficients correctly.', 'Error', 'error'));
            G = [];
        end
    end

    % Display the transfer function
    disp('The entered transfer function is:');
    disp(G);
end
