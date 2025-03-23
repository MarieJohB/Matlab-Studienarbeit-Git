function disturbance_cancel_input(fig)
% DISTURBANCE_CANCEL_INPUT - Handle cancellation
%
% Parameters:
%   fig - The figure to close

uiresume(fig);
delete(fig);
end