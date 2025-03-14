function html_str = createHtmlStringfortransfer(tf_sys, name)
    % createHtmlStringfortransfer
    % ---------------------------
    % Diese Funktion erstellt einen HTML-String, der eine Transferfunktion formatiert darstellt.
    % Dabei werden Zähler und Nenner in einen Bruch umgewandelt – mithilfe von CSS – sodass nur der
    % Bruchstrich angezeigt wird. Ein äußerer Container mit fester Breite wird verwendet, um
    % die Gesamtbreite zu definieren, während der Bruch selbst automatisch an seinen Inhalt angepasst und zentriert wird.
    %
    % Eingabeparameter:
    %   tf_sys - Das Transferfunktionssystem (tf-Objekt).
    %   name   - Der Name der Transferfunktion (z. B. 'G' oder 'K').
    %
    % Rückgabewert:
    %   html_str - Der HTML-String, der die formatierte Transferfunktion enthält.
    
    % Feste Breite für den äußeren Container (hier anpassbar)
    containerWidth = '300px'; 
    
    % Extrahiere die Koeffizienten der Transferfunktion (Zähler und Nenner)
    [num, den] = tfdata(tf_sys, 'v');
    
    % Entferne eventuelle führende Nullen, damit die Grade korrekt bestimmt werden.
    num = removeLeadingZeros(num);
    den = removeLeadingZeros(den);
    
    % Formatiere den Zähler und Nenner als HTML-kompatible Polynome
    num_str = formatPolynomial(num);
    den_str = formatPolynomial(den);
    
    % Ersetze wissenschaftliche Notation (z. B. "2.5e+3") durch HTML-kompatibles Format (Hochzahlen)
    num_str = replaceScientificNotation(num_str);
    den_str = replaceScientificNotation(den_str);
    
    % Feste Schriftgröße für die Darstellung des Bruchs
    fixedFontSize = '14px';
    
    % Erstelle den HTML-String:
    % 1. Ein äußerer Container (div) mit fester Breite containerWidth,
    %    der den gesamten Bruch zentriert.
    % 2. Innerhalb dieses Containers wird der Bruch als inline-flex Element
    %    dargestellt, sodass Zähler und Nenner (als separate Zeilen) zentriert sind.
    % 3. Der Name (z. B. "G(s) =") wird links davon angezeigt und vertikal mittig ausgerichtet.
    html_str = strcat( ...
        '<div style="display: flex; align-items: center;">', ...
            '<span style="font-size: 14px; margin-right: 5px; vertical-align: left;">', name, '(s)= </span>', ...
            '<div style="width: ', containerWidth, '; display: flex; justify-content: center; align-items: center;">', ...
                '<div style="display: inline-flex; flex-direction: column; text-align: center; font-size: ', fixedFontSize, ';">', ...
                    '<div style="border-bottom: 1px solid black; padding-bottom: 2px;">', num_str, '</div>', ...
                    '<div style="padding-top: 2px;">', den_str, '</div>', ...
                '</div>', ...
            '</div>', ...
        '</div>' ...
    );
end