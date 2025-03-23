function target_cancel_input(fig)
% TARGET_CANCEL_INPUT - Handle cancellation
%
% Parameters:
%   fig - The figure to close

uiresume(fig);
delete(fig);
end