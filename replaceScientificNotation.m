function formatted_str = replaceScientificNotation(input_str)
    % replaceScientificNotation
    % ---------------------------
    % Diese Funktion ersetzt in einem gegebenen String die wissenschaftliche
    % Notation (z. B. "2.5e+3") durch eine HTML-kompatible Darstellung, die
    % Exponenten als Hochzahlen (mittels <sup>...</sup>) darstellt.
    %
    % Eingabeparameter:
    %   input_str - Der Eingabestring, der wissenschaftliche Notation enthält.
    %
    % Rückgabewert:
    %   formatted_str - Der String mit ersetzter Notation.
    
    formatted_str = regexprep(input_str, 'e([-+]?\d+)', '&times;10<sup>$1</sup>');
end