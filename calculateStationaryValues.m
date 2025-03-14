function [y_inf_1, e_inf_1, y_inf_2, e_inf_2, e_inf_3, cansysjump] = calculateStationaryValues(T_s, S_s, R_s, D1_s, D2_s, G_s)
    % calculateStationaryValues berechnet die stationären Werte des Regelkreises.
    % Diese Funktion unterstützt sowohl kontinuierliche als auch abschnittsweise
    % definierte Funktionen für R, D1 und D2. Liegen Eingaben als double vor,
    % werden sie in tf-Objekte konvertiert. Dabei wird sichergestellt, dass
    % numerische Eingaben als Zeilenvektoren vorliegen, um Fehler bei tf zu vermeiden.
    %
    % \ac{Robustheit} und \ac{Flexibilität} werden durch die Unterstützung
    % unterschiedlicher Eingabetypen erreicht \cite{MatlabBestPractices2019}.
    
    % Sicherstellen, dass T_s, S_s und G_s als tf-Objekte vorliegen
    if ~isa(T_s, 'tf')
        T_s = tf(T_s, 1);
    end
    if ~isa(S_s, 'tf')
        S_s = tf(S_s, 1);
    end
    if ~isa(G_s, 'tf')
        G_s = tf(G_s, 1);
    end
    
    % Konvertiere T_s, S_s und G_s in symbolische Übertragungsfunktionen.
    % Verwende s = sym('s') anstelle von syms s, um Probleme in statischen Arbeitsbereichen zu vermeiden.
    s = sym('s');
    [T_num, T_den] = tfdata(T_s, 'v');
    [S_num, S_den] = tfdata(S_s, 'v');
    [G_num, G_den] = tfdata(G_s, 'v');
    
    T_sys = poly2sym(T_num, s) / poly2sym(T_den, s);
    S_sys = poly2sym(S_num, s) / poly2sym(S_den, s);
    G_sys = poly2sym(G_num, s) / poly2sym(G_den, s);
    
    % Hilfsfunktion: Sicherstellen, dass die Eingabe als Zell-Array vorliegt
    function cellOut = ensureCell(x)
        if iscell(x)
            cellOut = x;
        else
            cellOut = {x};
        end
    end

    % Sicherstellen, dass R_s, D1_s und D2_s als Zell-Arrays vorliegen,
    % um abschnittsweise und kontinuierliche Eingaben zu unterstützen.
    R_cell  = ensureCell(R_s);
    D1_cell = ensureCell(D1_s);
    D2_cell = ensureCell(D2_s);
    
    % Bestimme die Anzahl der Segmente anhand der längsten Zell-Array-Länge
    nSeg = max([numel(R_cell), numel(D1_cell), numel(D2_cell)]);
    
    % Initialisiere Ergebnisvektoren
    y_inf_1   = zeros(1, nSeg);
    e_inf_1   = zeros(1, nSeg);
    y_inf_2   = zeros(1, nSeg);
    e_inf_2   = zeros(1, nSeg);
    e_inf_3   = zeros(1, nSeg);
    cansysjump = zeros(1, nSeg);
    
    % Hilfsfunktion zur Auswahl des Elements: Falls das i-te Segment nicht existiert,
    % wird das erste Element (kontinuierliche Eingabe) verwendet.
    function out = selectSegment(cellArray, idx)
        if numel(cellArray) >= idx
            out = cellArray{idx};
        else
            out = cellArray{1};
        end
    end

    % Schleife über alle Segmente bzw. bei kontinuierlicher Eingabe über ein einzelnes Segment
    for i = 1:nSeg
        % Auswahl der entsprechenden Segmente (oder kontinuierliche Werte, falls kein Array)
        R_i  = selectSegment(R_cell, i);
        D1_i = selectSegment(D1_cell, i);
        D2_i = selectSegment(D2_cell, i);
        
        % Konvertiere die Eingaben in tf-Objekte, falls notwendig. Dabei wird sichergestellt,
        % dass numerische Eingaben als Zeilenvektoren vorliegen.
        if ~isa(R_i, 'tf')
            if isnumeric(R_i)
                R_i = tf(R_i(:)', 1);
            else
                error('R_i muss entweder numerisch oder ein tf-Objekt sein.');
            end
        end
        if ~isa(D1_i, 'tf')
            if isnumeric(D1_i)
                D1_i = tf(D1_i(:)', 1);
            else
                error('D1_i muss entweder numerisch oder ein tf-Objekt sein.');
            end
        end
        if ~isa(D2_i, 'tf')
            if isnumeric(D2_i)
                D2_i = tf(D2_i(:)', 1);
            else
                error('D2_i muss entweder numerisch oder ein tf-Objekt sein.');
            end
        end
        
        % Umwandlung in symbolische Übertragungsfunktionen
        [R_num, R_den]   = tfdata(R_i, 'v');
        [D1_num, D1_den] = tfdata(D1_i, 'v');
        [D2_num, D2_den] = tfdata(D2_i, 'v');
        
        R_sys  = poly2sym(R_num, s) / poly2sym(R_den, s);
        D1_sys = poly2sym(D1_num, s) / poly2sym(D1_den, s);
        D2_sys = poly2sym(D2_num, s) / poly2sym(D2_den, s);
        
        % Berechnung der stationären Werte für das aktuelle Segment
        y_inf_1(i)    = double(limit(s * T_sys * R_sys, s, 0));
        e_inf_1(i)    = double(limit(s * S_sys * R_sys, s, 0));
        y_inf_2(i)    = double(limit(s * S_sys * D1_sys, s, 0));
        e_inf_2(i)    = double(limit(-s * S_sys * D1_sys, s, 0));
        e_inf_3(i)    = double(limit(s * (S_sys * R_sys - G_sys * S_sys * D2_sys), s, 0));
        cansysjump(i) = double(limit(s * T_sys * R_sys, s, inf));
    end

    % Falls nur ein Segment vorliegt, werden skalare Werte zurückgegeben
    if nSeg == 1
        y_inf_1    = y_inf_1(1);
        e_inf_1    = e_inf_1(1);
        y_inf_2    = y_inf_2(1);
        e_inf_2    = e_inf_2(1);
        e_inf_3    = e_inf_3(1);
        cansysjump = cansysjump(1);
    end
end
