function [values_vector] = generate_Function_Values_Version2_with_time(function_vector, time_vector)
    % Anzahl der Abschnitte prüfen
    num_sections = length(function_vector);
    
    % Zeitvektor verwenden
    t = time_vector;
    
    % Wertevektor initialisieren
    values_vector = zeros(size(t));
    
    if num_sections == 1 % Funktion ist kontinuierlich
        u = function_vector;
        
        if isa(u, 'function_handle')
            values_vector = arrayfun(u, t);
        else
            values_vector = u * ones(size(t));
        end
    elseif num_sections > 1 % Funktion ist stückweise definiert
        % Abrufen der Zeitparameter für die stückweise Funktion
        [start_time, end_time, ~] = get_time_vector_disturbance(num_sections);
        
        % Für jeden Zeitpunkt den korrekten Funktionswert bestimmen
        for i = 1:length(t)
            current_time = t(i);
            
            % Bestimmen des aktiven Abschnitts für diesen Zeitpunkt
            active_section = 0;
            for j = 1:num_sections
                if current_time >= start_time(j) && current_time <= end_time(j)
                    active_section = j;
                    break;
                end
            end
            
            % Wenn ein aktiver Abschnitt gefunden wurde, Wert berechnen
            if active_section > 0
                if isa(function_vector{active_section}, 'function_handle')
                    values_vector(i) = function_vector{active_section}(current_time);
                else
                    values_vector(i) = function_vector{active_section};
                end
            end
        end
    else
        disp('Something went wrong. Number of sections is smaller than 1.');
        values_vector = [];
        return;
    end
    
    % Debug-Info
    disp('Debugging Info:');
    disp(['Length of t: ', num2str(length(t))]);
    disp(['Length of values_vector: ', num2str(length(values_vector))]);
end