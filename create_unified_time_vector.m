%% 2. Funktion zum Erstellen eines gemeinsamen Zeitvektors
function [unified_time_vector] = create_unified_time_vector(r_times, r_start_times, r_end_times, d1_times, d1_start_times, d1_end_times, d2_times, d2_start_times, d2_end_times)
    % CREATE_UNIFIED_TIME_VECTOR - Erzeugt einen einheitlichen Zeitvektor für Simulationen
    % Diese Funktion kombiniert die Zeitvektoren von Sollwert und Störgrößen in einem
    % einheitlichen Zeitvektor für die Simulation.
    %
    % Eingaben (alle optional):
    %   r_times - Zeitvektor für Sollwert
    %   r_start_times, r_end_times - Start- und Endzeiten der Sollwert-Abschnitte
    %   d1_times - Zeitvektor für erste Störgröße
    %   d1_start_times, d1_end_times - Start- und Endzeiten der d1-Abschnitte
    %   d2_times - Zeitvektor für zweite Störgröße
    %   d2_start_times, d2_end_times - Start- und Endzeiten der d2-Abschnitte
    %
    % Ausgabe:
    %   unified_time_vector - Kombinierter Zeitvektor, der alle Übergangspunkte enthält

    % Initialisieren aller Start- und Endzeiten
    all_transition_times = [];
    
    % Hinzufügen der Sollwert-Zeitübergänge
    if nargin >= 3 && ~isempty(r_start_times) && ~isempty(r_end_times)
        all_transition_times = [all_transition_times, r_start_times, r_end_times];
    elseif nargin >= 1 && ~isempty(r_times)
        all_transition_times = [all_transition_times, r_times(1), r_times(end)];
    end
    
    % Hinzufügen der d1-Zeitübergänge
    if nargin >= 6 && ~isempty(d1_start_times) && ~isempty(d1_end_times)
        all_transition_times = [all_transition_times, d1_start_times, d1_end_times];
    elseif nargin >= 4 && ~isempty(d1_times)
        all_transition_times = [all_transition_times, d1_times(1), d1_times(end)];
    end
    
    % Hinzufügen der d2-Zeitübergänge
    if nargin >= 9 && ~isempty(d2_start_times) && ~isempty(d2_end_times)
        all_transition_times = [all_transition_times, d2_start_times, d2_end_times];
    elseif nargin >= 7 && ~isempty(d2_times)
        all_transition_times = [all_transition_times, d2_times(1), d2_times(end)];
    end
    
    % Falls keine Zeiten gegeben sind, Standardzeitvektor zurückgeben
    if isempty(all_transition_times)
        unified_time_vector = linspace(0, 10, 1000);
        return;
    end
    
    % Einzigartige Zeitpunkte extrahieren und sortieren
    all_transition_times = unique(all_transition_times);
    all_transition_times = sort(all_transition_times);
    
    % Mindestens Anfang und Ende haben
    if length(all_transition_times) < 2
        if length(all_transition_times) == 1
            % Nur ein Zeitpunkt, füge einen weiteren hinzu
            all_transition_times = [all_transition_times, all_transition_times(1) + 10];
        else
            % Keine Zeitpunkte, Standardzeitvektor zurückgeben
            unified_time_vector = linspace(0, 10, 1000);
            return;
        end
    end
    
    % Bestimmen der feinsten Schrittweite aus den vorhandenen Zeitvektoren
    min_step = 0.1; % Standardwert
    
    if ~isempty(r_times) && length(r_times) > 1
        r_step = min(diff(r_times));
        min_step = min(min_step, r_step);
    end
    
    if ~isempty(d1_times) && length(d1_times) > 1
        d1_step = min(diff(d1_times));
        min_step = min(min_step, d1_step);
    end
    
    if ~isempty(d2_times) && length(d2_times) > 1
        d2_step = min(diff(d2_times));
        min_step = min(min_step, d2_step);
    end
    
    % Erzeugen des einheitlichen Zeitvektors
    unified_time_vector = [];
    
    % Für jedes Übergangsintervall einen Zeitvektor mit der minimalen Schrittweite erzeugen
    for i = 1:length(all_transition_times)-1
        start_t = all_transition_times(i);
        end_t = all_transition_times(i+1);
        
        % Sicherstellen, dass das Intervall nicht zu klein ist
        if end_t - start_t < min_step
            % Überspringen von zu kleinen Intervallen
            continue;
        end
        
        % Anzahl der Punkte im Intervall berechnen
        num_points = max(2, ceil((end_t - start_t) / min_step) + 1);
        
        % Zeitvektor für dieses Intervall erstellen
        interval_times = linspace(start_t, end_t, num_points);
        
        % Wenn nicht das erste Intervall, den ersten Punkt entfernen, um Duplikate zu vermeiden
        if ~isempty(unified_time_vector) && ~isempty(interval_times)
            interval_times = interval_times(2:end);
        end
        
        % Zum einheitlichen Zeitvektor hinzufügen
        unified_time_vector = [unified_time_vector, interval_times];
    end
    
    % Sicherstellen, dass der Zeitvektor nicht leer ist
    if isempty(unified_time_vector)
        unified_time_vector = linspace(all_transition_times(1), all_transition_times(end), 1000);
    end
end