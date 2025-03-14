function formatted_poly = formatPolynomial(coeffs)
    % formatPolynomial
    % ----------------
    % Diese Funktion formatiert einen Koeffizientenvektor (Polynom) in einen
    % HTML-kompatiblen String. Es werden Terme mit Variablen (s) und Hochzahlen
    % für Potenzen generiert.
    %
    % Eingabeparameter:
    %   coeffs - Vektor mit Koeffizienten des Polynoms.
    %
    % Rückgabewert:
    %   formatted_poly - Der formatierte HTML-String des Polynoms.
    
    % Falls das Polynom ausschließlich Nullen enthält, gebe "0" zurück.
    if isempty(coeffs) || all(coeffs == 0)
        formatted_poly = '0';
        return;
    end

    deg = length(coeffs) - 1;  % Bestimme den Grad des Polynoms
    terms = {};                % Zell-Array zum Speichern einzelner Terme

    % Iteriere über alle Koeffizienten
    for i = 1:length(coeffs)
        coeff = coeffs(i);
        exp = deg - (i - 1);   % Bestimme den Exponenten für den aktuellen Term

        % Überspringe Nullen, da sie keinen Einfluss auf den Ausdruck haben
        if coeff == 0
            continue;
        end

        % Vorzeichenbehandlung: Füge ein " + " hinzu, wenn es nicht der erste Term ist;
        % Bei negativen Koeffizienten wird ein " - " gesetzt und der Betrag genutzt.
        if coeff > 0 && ~isempty(terms)
            term = ' + ';
        elseif coeff < 0
            term = ' - ';
            coeff = abs(coeff);  % Verwende den absoluten Wert für die Anzeige
        else
            term = '';
        end

        % Zeige den Koeffizienten an, außer er ist 1 und es handelt sich nicht um den
        % konstanten Term (exp==0)
        if coeff ~= 1 || exp == 0
            term = strcat(term, num2str(coeff));
        end

        % Hänge das "s" und ggf. eine Hochzahl hinzu, falls der Exponent größer als 0 ist.
        if exp > 1
            term = strcat(term, 's<sup>', num2str(exp), '</sup>');
        elseif exp == 1
            term = strcat(term, 's');
        end

        % Füge den formatierten Term zum Zell-Array hinzu.
        terms{end + 1} = term;
    end

    % Verbinde alle Terme zu einem einzelnen String
    formatted_poly = strjoin(terms, '');
end