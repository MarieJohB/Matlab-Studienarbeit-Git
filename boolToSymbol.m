function symbol = boolToSymbol(x)
    % boolToSymbol
    % --------------
    % Diese Funktion wandelt einen booleschen Wert in ein HTML-Symbol um.
    % Bei true wird ein Häkchen (✓) mit grünem Hintergrund angezeigt,
    % bei false ein Kreuz (✗) mit rotem Hintergrund.
    
    if x
        % Häkchen (✓) mit grünem Hintergrund
        symbol = '<span style="font-size: 14px; background-color: #90EE90; padding: 2px; border-radius: 2px;">&#10003;</span>';
    else
        % Kreuz (✗) mit rotem Hintergrund
        symbol = '<span style="font-size: 14px; background-color: #FFCCCC; padding: 2px; border-radius: 2px;">&#10007;</span>';
    end
    
    % Alternative: Einfache Textversion für nicht-HTML-Kontexte
    % symbol = x ? "✓" : "✗";
end