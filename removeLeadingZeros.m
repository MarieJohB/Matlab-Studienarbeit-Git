function polyOut = removeLeadingZeros(polyIn)
    % removeLeadingZeros
    % ------------------
    % Diese Funktion entfernt führende Nullen aus einem Koeffizientenvektor,
    % damit die korrekte Gradbestimmung des Polynoms möglich ist.
    %
    % Eingabeparameter:
    %   polyIn - Koeffizientenvektor des Polynoms.
    %
    % Rückgabewert:
    %   polyOut - Der modifizierte Koeffizientenvektor ohne führende Nullen.
    
    idx = find(polyIn, 1, 'first');  % Finde den Index des ersten ungleich Null-Elements
    if isempty(idx)
        polyOut = 0;  % Falls alle Werte null sind, gebe 0 zurück
    else
        polyOut = polyIn(idx:end);  % Andernfalls den Vektor ab dem ersten ungleich Null-Wert zurückgeben
    end
end