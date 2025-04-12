function tfText = getTransferFunctionText(tf_obj, name)
% Get a text representation of a transfer function
[num, den] = tfdata(tf_obj, 'v');

numStr = coeffsToString(num);
denStr = coeffsToString(den);

tfText = sprintf('%s(s) = (%s)/(%s)', name, numStr, denStr);
end