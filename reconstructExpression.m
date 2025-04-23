function expr = reconstructExpression(op, args)
    % Helper function to reconstruct symbolic expression from operation and arguments
    
    switch op
        case 'add'
            expr = args{1} + args{2};
        case 'subtract'
            expr = args{1} - args{2};
        case 'multiply'
            expr = args{1} * args{2};
        case 'divide'
            expr = args{1} / args{2};
        case 'power'
            expr = args{1} ^ args{2};
        case 'sin'
            expr = sin(args{1});
        case 'cos'
            expr = cos(args{1});
        case 'tan'
            expr = tan(args{1});
        case 'exp'
            expr = exp(args{1});
        case 'log'
            expr = log(args{1});
        case 'function'
            expr = symfun(args{1}, args{2});
        otherwise
            expr = args{1};
    end
end