function [op, args] = extractOperationAndArgs(expr)
    % Helper function to extract operation and arguments from symbolic expression
    
    if isa(expr, 'sym') && ~isempty(symvar(expr))
        % Try to identify the operation
        if isa(expr, 'symfun')
            % Handle function applications
            op = 'function';
            args = {expr.func, expr.vars};
        elseif isEquivalent(expr, sin(symvar(expr)))
            op = 'sin';
            % Extract the argument of sin (which may be more complex than just a variable)
            args = {extractFunctionArg(expr, 'sin')};
        elseif isEquivalent(expr, cos(symvar(expr)))
            op = 'cos';
            args = {extractFunctionArg(expr, 'cos')};
        elseif isEquivalent(expr, tan(symvar(expr)))
            op = 'tan';
            args = {extractFunctionArg(expr, 'tan')};
        elseif isEquivalent(expr, exp(symvar(expr)))
            op = 'exp';
            args = {extractFunctionArg(expr, 'exp')};
        elseif isEquivalent(expr, log(symvar(expr)))
            op = 'log';
            args = {extractFunctionArg(expr, 'log')};
        else
            try
                % For addition, multiplication, etc.
                C = children(expr);
                if length(C) == 2
                    if isEquivalent(expr, C{1} + C{2})
                        op = 'add';
                        args = {C{1}, C{2}};
                    elseif isEquivalent(expr, C{1} - C{2})
                        op = 'subtract';
                        args = {C{1}, C{2}};
                    elseif isEquivalent(expr, C{1} * C{2})
                        op = 'multiply';
                        args = {C{1}, C{2}};
                    elseif isEquivalent(expr, C{1} / C{2})
                        op = 'divide';
                        args = {C{1}, C{2}};
                    elseif isEquivalent(expr, C{1} ^ C{2})
                        op = 'power';
                        args = {C{1}, C{2}};
                    else
                        op = '';
                        args = {};
                    end
                else
                    op = '';
                    args = {};
                end
            catch
                op = '';
                args = {};
            end
        end
    else
        % Simple expression or number
        op = '';
        args = {};
    end
end