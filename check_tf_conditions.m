function [isStable, isProper, isStrictlyProper] = check_tf_conditions(tf_sys)
    % check_tf_conditions überprüft, ob das Übertragungssystem tf_sys:
    % - stabil ist (alle Pole haben negative Realteile),
    % - proper ist (Grad des Zählers <= Grad des Nenners) und
    % - strictly proper (Grad des Zählers < Grad des Nenners).
    %
    % Es werden boolesche Werte zurückgegeben, die später in der UI (z.B. in einer
    % HTML-Tabelle) in entsprechende Symbole (z.B. Häkchen oder Kreuze) umgewandelt werden können.
    
    % Extrahiere die Koeffizienten der Zähler- und Nennerpolynome
    [num, den] = tfdata(tf_sys, 'v');
    num = removeLeadingZeros(num);
    den = removeLeadingZeros(den);
    
    % Bestimme den Grad (Anzahl der Koeffizienten minus 1)
    degNum = length(num) - 1;
    degDen = length(den) - 1;
    
    % Properness: Zählergrad <= Nennergrad
    isProper = (degNum <= degDen);
    
    % Strictly Proper: Zählergrad < Nennergrad
    isStrictlyProper = (degNum < degDen);
    
    % Stabilität: Alle Pole haben negative Realteile
    p = pole(tf_sys);
    isStable = all(real(p) < 0);
end