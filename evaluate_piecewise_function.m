% Dieses Skript enthält verbesserte Funktionen für das Management von Zeitvektoren 
% und stückweisen Funktionen in der Regelkreis-App

%% 1. Verbesserte stückweise Funktionsauswertung
function [values_vector] = evaluate_piecewise_function(function_vector, time_vector, start_times, end_times)
    % EVALUATE_PIECEWISE_FUNCTION - Wertet eine stückweise Funktion für einen gegebenen Zeitvektor aus
    % Diese Funktion behandelt sowohl kontinuierliche als auch stückweise Funktionen.
    %
    % Eingaben:
    %   function_vector - Entweder ein function_handle (kontinuierlich) oder eine Zelle von function_handles
    %   time_vector - Der Zeitvektor, für den die Funktion ausgewertet werden soll
    %   start_times - Startzeiten für jeden Abschnitt (nur für stückweise Funktionen)
    %   end_times - Endzeiten für jeden Abschnitt (nur für stückweise Funktionen)
    %
    % Ausgabe:
    %   values_vector - Ausgewertete Funktionswerte für jeden Zeitpunkt

    % Initialisieren des Ergebnisvektors
    values_vector = zeros(size(time_vector));
    
    % Überprüfen, ob es sich um eine kontinuierliche oder stückweise Funktion handelt
    if ~iscell(function_vector)
        % Kontinuierlicher Fall - einfache Auswertung
        if isa(function_vector, 'function_handle')
            values_vector = arrayfun(function_vector, time_vector);
        else
            % Konstanter Wert
            values_vector = function_vector * ones(size(time_vector));
        end
    else
        % Stückweise Funktion - Abschnitt für Abschnitt auswerten
        num_sections = length(function_vector);
        
        % Für jeden Zeitpunkt bestimmen, welcher Abschnitt aktiv ist
        for i = 1:length(time_vector)
            current_time = time_vector(i);
            
            % Finden des aktiven Abschnitts für diesen Zeitpunkt
            section_idx = find(current_time >= start_times & current_time <= end_times, 1, 'first');
            
            % Wenn kein Abschnitt gefunden wurde, den nächstgelegenen wählen
            if isempty(section_idx)
                % Finden des nächstgelegenen Abschnitts (extrapolieren)
                if current_time < start_times(1)
                    % Zeit liegt vor dem ersten Abschnitt
                    section_idx = 1;
                elseif current_time > end_times(end)
                    % Zeit liegt nach dem letzten Abschnitt
                    section_idx = num_sections;
                end
            end
            
            % Auswerten der Funktion für den aktuellen Zeitpunkt
            if ~isempty(section_idx) && section_idx <= num_sections
                if isa(function_vector{section_idx}, 'function_handle')
                    values_vector(i) = function_vector{section_idx}(current_time);
                else
                    % Konstanter Wert für diesen Abschnitt
                    values_vector(i) = function_vector{section_idx};
                end
            end
        end
    end
end



