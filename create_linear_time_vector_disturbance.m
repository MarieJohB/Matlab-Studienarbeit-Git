function [t] = create_linear_time_vector_disturbance(num_sections, disturbance_label)
    % CREATE_LINEAR_TIME_VECTOR_DISTURBANCE - Erzeugt einen linearen Zeitvektor für eine Störgröße
    % Diese Funktion ist an create_linear_time_vector angelehnt, berücksichtigt aber
    % die spezifischen Zeitparameter der Störgröße.
    
    % Zeitvektor initialisieren
    t = [];
    
    % Prüfen, ob num_sections gültig ist
    if isempty(num_sections)
        disp('Operation cancelled: No sections defined.');
        return;
    end
    
    % Zeitparameter vom get_time_vector_disturbance abrufen
    [start_time, end_time, time_steps] = get_time_vector_disturbance(num_sections, [], [], [], disturbance_label);
    
    % Prüfen, ob Zeitparameter zurückgegeben wurden
    if isempty(start_time) || isempty(end_time) || isempty(time_steps)
        disp('Operation cancelled: Time parameters not provided.');
        return;
    end
    
    % Zeitvektor für jeden Abschnitt erstellen
    for i = 1:num_sections
        % Anzahl der Punkte in diesem Abschnitt berechnen
        num_points = round((end_time(i) - start_time(i)) / time_steps(i)) + 1;
        
        % Zeitvektor für diesen Abschnitt erstellen
        section_t = linspace(start_time(i), end_time(i), num_points);
        
        % Wenn dies nicht der erste Abschnitt ist, den ersten Punkt entfernen, um Duplizierung zu vermeiden
        if i > 1 && ~isempty(t) && ~isempty(section_t)
            section_t = section_t(2:end);
        end
        
        % Dem gesamten Zeitvektor hinzufügen
        t = [t, section_t];
    end
    
    % Informationen zum Zeitvektor anzeigen (nur wenn wir einen gültigen Zeitvektor haben)
    if ~isempty(t)
        disp(['Created ', disturbance_label, ' time vector with ', num2str(length(t)), ' points from ', ...
            num2str(t(1)), ' to ', num2str(t(end)), ' across ', num2str(num_sections), ' section(s).']);
    else
        disp(['Warning: Created an empty time vector for ', disturbance_label, '.']);
    end
end