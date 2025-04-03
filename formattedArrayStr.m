% Format array to string with spacing
function str = formattedArrayStr(arr)
    str = '';
    for i = 1:length(arr)
        if i > 1
            str = [str, ' '];
        end
        str = [str, sprintf('%.4g', arr(i))];
    end
end