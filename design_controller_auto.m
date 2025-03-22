function [K, details] = design_controller_auto(G, method, structure, options)
% DESIGN_CONTROLLER_AUTO Automatischer Reglerentwurf mit verschiedenen Methoden
% 
% Eingabe:
%   G         - Regelstrecke als Übertragungsfunktion
%   method    - Entwurfsmethode (String)
%   structure - Reglerstruktur: 'P', 'PI', 'PD', 'PID'
%   options   - Struktur mit optionalen Parametern:
%     .epsilon     - Filter-Parameter für D-Anteil (default: 0.1)
%     .phaseMargin - Gewünschte Phasenreserve in Grad (default: 45)
%     .bandwidth   - Gewünschte Bandbreite in rad/s (default: 1)
%     .settlingTime - Gewünschte Einschwingzeit in s (default: 5)
%     .robustness  - Robustheitsgrad: 'Low', 'Medium', 'High' (default: 'Medium')
%     .overshoot   - Gewünschtes Überschwingen in % (default: 10)
%     .goal        - Optimierungsziel: 'Tracking', 'Disturbance Rejection', 'Robustness' (default: 'Tracking')
%
% Ausgabe:
%   K        - Entworfener Regler als Übertragungsfunktion
%   details  - Textbeschreibung mit Details zum Entwurf

    % Default-Werte für fehlende Optionen
    if nargin < 4
        options = struct();
    end
    
    if ~isfield(options, 'epsilon')
        options.epsilon = 0.1;
    end
    
    if ~isfield(options, 'phaseMargin')
        options.phaseMargin = 45;
    end
    
    if ~isfield(options, 'bandwidth')
        options.bandwidth = 1;
    end
    
    if ~isfield(options, 'settlingTime')
        options.settlingTime = 5;
    end
    
    if ~isfield(options, 'robustness')
        options.robustness = 'Medium';
    end
    
    if ~isfield(options, 'overshoot')
        options.overshoot = 10;
    end
    
    if ~isfield(options, 'goal')
        options.goal = 'Tracking';
    end
    
    % Auswahl der Entwurfsmethode
    switch method
        case 'Ziegler-Nichols (Oscillation)'
            [K, details] = designZieglerNicholsOscillation(G, structure, options.epsilon);
        case 'Ziegler-Nichols (Step)'
            [K, details] = designZieglerNicholsStep(G, structure, options.epsilon);
        case 'Aström'
            [K, details] = designAstrom(G, structure, options.epsilon);
        case 'CHR'
            [K, details] = designCHR(G, structure, options.epsilon);
        case 'Cohen-Coon'
            [K, details] = designCohenCoon(G, structure, options.epsilon);
        case 'Loop-Shaping'
            [K, details] = designLoopShaping(G, structure, options.phaseMargin, options.bandwidth, options.epsilon);
        case 'IMC'
            [K, details] = designIMC(G, structure, options.settlingTime, options.epsilon);
        case 'MIGO'
            [K, details] = designMIGO(G, structure, options.robustness, options.epsilon);
        case 'H-infinity'
            [K, details] = designHInfinity(G, structure, options.robustness, options.epsilon);
        case 'LQG (Linear-Quadratic-Gaussian)'
            [K, details] = designLQG(G, structure, options.bandwidth, options.robustness, options.epsilon);
        otherwise
            error('Unbekannte Entwurfsmethode: %s', method);
    end
    
    % Evaluate controller quality if requested
    try
        score = evaluateController(K, G, options.goal, options.phaseMargin, options.overshoot, options.settlingTime, options.bandwidth);
        details = [details, sprintf('\n\nController Score: %.2f/100', score)];
    catch ME
        disp(['Warning: Could not evaluate controller quality: ', ME.message]);
    end
end

% 1. Ziegler-Nichols-Schwingungsmethode
function [K, details] = designZieglerNicholsOscillation(G, structure, epsilon)
    % Bestimme kritische Verstärkung und Periodendauer
    k_krit = 0;
    T_krit = 0;
    
    % Starte mit einer kleinen Verstärkung und erhöhe sie iterativ
    try
        found = false;
        step_size = 0.1;
        max_iterations = 2000;
        
        for iteration = 1:max_iterations
            k = iteration * step_size;
            
            % Teste Stabilität mit aktuellem k
            closed_loop = feedback(G*k, 1);
            poles = pole(closed_loop);
            
            % Prüfe auf Pole auf der imaginären Achse
            realParts = real(poles);
            imagParts = imag(poles);
            
            close_to_imag_axis = find(abs(realParts) < 0.001 & imagParts ~= 0);
            
            if ~isempty(close_to_imag_axis)
                % Gefundene Pole nahe der Stabilitätsgrenze
                k_krit = k;
                
                % Berechne Periodendauer
                idx = find(abs(realParts) < 0.001 & imagParts > 0, 1);
                if ~isempty(idx)
                    omega = imagParts(idx);
                    T_krit = 2*pi/omega;
                    found = true;
                    break;
                end
            elseif any(realParts > 0)
                % System ist instabil geworden
                
                % Reduziere k schrittweise, bis wir nahe an der Grenze sind
                for back_step = 1:10
                    k_test = k - back_step * (step_size/10);
                    
                    closed_loop = feedback(G*k_test, 1);
                    poles = pole(closed_loop);
                    realParts = real(poles);
                    imagParts = imag(poles);
                    
                    if all(realParts < 0) && any(abs(realParts) < 0.01 & imagParts ~= 0)
                        k_krit = k_test;
                        
                        idx = find(abs(realParts) < 0.01 & imagParts > 0, 1);
                        if ~isempty(idx)
                            omega = imagParts(idx);
                            T_krit = 2*pi/omega;
                            found = true;
                            break;
                        end
                    end
                end
                
                if found
                    break;
                else
                    % Annahme: Vergangener kritischer Punkt
                    k_krit = k - step_size;
                    
                    closed_loop = feedback(G*k_krit, 1);
                    poles = pole(closed_loop);
                    imagParts = imag(poles(abs(real(poles)) < 0.1));
                    
                    if ~isempty(imagParts) && any(imagParts > 0)
                        omega = max(imagParts(imagParts > 0));
                        T_krit = 2*pi/omega;
                        found = true;
                        break;
                    end
                end
                
                break;
            end
        end
        
        if ~found || k_krit == 0
            error('Could not determine critical gain. The plant may not be suitable for the Ziegler-Nichols oscillation method.');
        end
    catch ME
        error('Error applying Ziegler-Nichols oscillation method: %s', ME.message);
    end
    
    % Reglerparameter nach Ziegler-Nichols-Tabelle berechnen
    details = sprintf('k_krit = %.4f\nT_krit = %.4f s', k_krit, T_krit);
    
    switch structure
        case 'P'
            Kp = 0.5 * k_krit;
            K = tf(Kp, 1);
        case 'PI'
            Kp = 0.45 * k_krit;
            Ti = 0.85 * T_krit;
            K = tf([Kp, Kp/Ti], [1, 0]);
        case 'PD'
            Kp = 0.5 * k_krit;
            Td = 0.12 * T_krit;
            K = tf([Kp*Td, Kp], [epsilon*Td, 1]);
        case 'PID'
            Kp = 0.6 * k_krit;
            Ti = 0.5 * T_krit;
            Td = 0.125 * T_krit;
            K = tf([Kp*Td, Kp, Kp/Ti], [epsilon*Td, 1, 0]);
        otherwise
            error('Unsupported controller structure for Ziegler-Nichols oscillation method');
    end
end

% 2. Ziegler-Nichols-Sprungmethode
function [K, details] = designZieglerNicholsStep(G, structure, epsilon)
    % Berechne Sprungantwort
    t = linspace(0, 100, 1000);
    [y, t] = step(G, t);
    
    % Finde Endwert
    y_final = y(end);
    
    if abs(y_final) < 1e-6
        error('System does not have a finite DC gain. The plant may not be suitable for the Ziegler-Nichols step method.');
    end
    
    % Normalisiere die Antwort
    y_norm = y / y_final;
    
    % Bestimme die Tangente am Wendepunkt
    dy = diff(y_norm) ./ diff(t);
    
    % Finde den Punkt mit maximaler Steigung (Wendepunkt)
    [max_slope, idx_max_slope] = max(dy);
    
    % Berechne Parameter der Tangente: y = m*t + b
    m = max_slope;
    b = y_norm(idx_max_slope) - m * t(idx_max_slope);
    
    % Schnittpunkte der Tangente mit den Linien y=0 und y=1
    t_0 = -b / m;                 % Schnittpunkt mit y=0 (t-Achse)
    t_1 = (1 - b) / m;            % Schnittpunkt mit y=1
    
    % Bestimme Totzeit L und Anstiegszeit T
    L = t_0;
    T = t_1 - t_0;
    
    % Stationäre Verstärkung
    Ks = dcgain(G);
    
    % Parameter nach Ziegler-Nichols-Tabelle berechnen
    details = sprintf('Ks = %.4f\nL = %.4f s\nT = %.4f s', Ks, L, T);
    
    switch structure
        case 'P'
            Kp = T / (Ks * L);
            K = tf(Kp, 1);
        case 'PI'
            Kp = 0.9 * T / (Ks * L);
            Ti = L / 0.3;
            K = tf([Kp, Kp/Ti], [1, 0]);
        case 'PD'
            Kp = 1.2 * T / (Ks * L);
            Td = 0.5 * L;
            K = tf([Kp*Td, Kp], [epsilon*Td, 1]);
        case 'PID'
            Kp = 1.2 * T / (Ks * L);
            Ti = 2 * L;
            Td = 0.5 * L;
            K = tf([Kp*Td, Kp, Kp/Ti], [epsilon*Td, 1, 0]);
        otherwise
            error('Unsupported controller structure for Ziegler-Nichols step method');
    end
end

% 3. Aström-Methode
function [K, details] = designAstrom(G, structure, epsilon)
    % Berechne Sprungantwort
    t = linspace(0, 100, 1000);
    [y, t] = step(G, t);
    
    % Finde Endwert
    y_final = y(end);
    
    if abs(y_final) < 1e-6
        error('System does not have a finite DC gain. The plant may not be suitable for the Aström method.');
    end
    
    % Stationäre Verstärkung
    Ks = dcgain(G);
    
    % Suche 63.2% Punkt für T
    idx_63 = find(y >= 0.632*y_final, 1);
    
    if isempty(idx_63)
        error('Could not determine time constant. The plant may not be suitable for the Aström method.');
    end
    
    % Schätze Totzeit durch Vergleich mit Übertragungsfunktion erster Ordnung mit Totzeit
    idx_10 = find(y >= 0.1*y_final, 1);
    if isempty(idx_10)
        L = 0.1;  % Default-Wert
    else
        L = t(idx_10);  % Schätze Totzeit als Zeit beim 10% Anstieg
    end
    
    T = t(idx_63) - L;  % Zeitkonstante (aus 63.2% abzüglich Totzeit)
    
    details = sprintf('Ks = %.4f\nL = %.4f s\nT = %.4f s', Ks, L, T);
    
    % Reglerparameter nach Aström berechnen
    switch structure
        case 'P'
            if L < 2*T
                Kp = 0.3 * T / (Ks * L);
            else
                Kp = 0.15 / Ks;
            end
            K = tf(Kp, 1);
        case 'PI'
            if L < 2*T
                Kp = 0.3 * T / (Ks * L);
            else
                Kp = 0.15 / Ks;
            end
            
            if L < 0.1*T
                Ti = 8 * L;
            elseif L < 2*T
                Ti = 0.8 * T;
            else
                Ti = 0.4 * L;
            end
            
            K = tf([Kp, Kp/Ti], [1, 0]);
        case 'PD'
            Kp = 0.4 * T / (Ks * L);
            Td = 0.4 * L;
            K = tf([Kp*Td, Kp], [epsilon*Td, 1]);
        case 'PID'
            if L < 2*T
                Kp = 0.3 * T / (Ks * L);
            else
                Kp = 0.15 / Ks;
            end
            
            if L < 0.1*T
                Ti = 8 * L;
                Td = 0.25 * L;
            elseif L < 2*T
                Ti = 0.8 * T;
                Td = 0.2 * L;
            else
                Ti = 0.4 * L;
                Td = 0.15 * L;
            end
            
            K = tf([Kp*Td, Kp, Kp/Ti], [epsilon*Td, 1, 0]);
        otherwise
            error('Unsupported controller structure for Aström method');
    end
end

% 4. CHR-Methode (Chien-Hrones-Reswick) - mit 0% Überschwingen
function [K, details] = designCHR(G, structure, epsilon)
    % Berechne Sprungantwort
    t = linspace(0, 100, 1000);
    [y, t] = step(G, t);
    
    % Finde Endwert
    y_final = y(end);
    
    if abs(y_final) < 1e-6
        error('System does not have a finite DC gain. The plant may not be suitable for the CHR method.');
    end
    
    % Normalisiere die Antwort
    y_norm = y / y_final;
    
    % Bestimme die Tangente am Wendepunkt
    dy = diff(y_norm) ./ diff(t);
    
    % Finde den Punkt mit maximaler Steigung (Wendepunkt)
    [max_slope, idx_max_slope] = max(dy);
    
    % Berechne Parameter der Tangente: y = m*t + b
    m = max_slope;
    b = y_norm(idx_max_slope) - m * t(idx_max_slope);
    
    % Schnittpunkte der Tangente mit den Linien y=0 und y=1
    t_0 = -b / m;                 % Schnittpunkt mit y=0 (t-Achse)
    t_1 = (1 - b) / m;            % Schnittpunkt mit y=1
    
    % Bestimme Totzeit L und Anstiegszeit T
    L = t_0;
    T = t_1 - t_0;
    
    % Stationäre Verstärkung
    Ks = dcgain(G);
    
    details = sprintf('Ks = %.4f\nL = %.4f s\nT = %.4f s', Ks, L, T);
    
    % Parameter nach CHR-Tabelle für Sollwertregelung (0% Überschwingen)
    switch structure
        case 'P'
            Kp = 0.3 * T / (Ks * L);
            K = tf(Kp, 1);
        case 'PI'
            Kp = 0.35 * T / (Ks * L);
            Ti = 1.2 * T;
            K = tf([Kp, Kp/Ti], [1, 0]);
        case 'PD'
            Kp = 0.5 * T / (Ks * L);
            Td = 0.3 * L;
            K = tf([Kp*Td, Kp], [epsilon*Td, 1]);
        case 'PID'
            Kp = 0.6 * T / (Ks * L);
            Ti = T;
            Td = 0.5 * L;
            K = tf([Kp*Td, Kp, Kp/Ti], [epsilon*Td, 1, 0]);
        otherwise
            error('Unsupported controller structure for CHR method');
    end
end

% 5. Cohen-Coon-Methode
function [K, details] = designCohenCoon(G, structure, epsilon)
    % Berechne Sprungantwort
    t = linspace(0, 100, 1000);
    [y, t] = step(G, t);
    
    % Finde Endwert
    y_final = y(end);
    
    if abs(y_final) < 1e-6
        error('System does not have a finite DC gain. The plant may not be suitable for the Cohen-Coon method.');
    end
    
    % Normalisiere die Antwort
    y_norm = y / y_final;
    
    % Bestimme die Tangente am Wendepunkt
    dy = diff(y_norm) ./ diff(t);
    
    % Finde den Punkt mit maximaler Steigung (Wendepunkt)
    [max_slope, idx_max_slope] = max(dy);
    
    % Berechne Parameter der Tangente: y = m*t + b
    m = max_slope;
    b = y_norm(idx_max_slope) - m * t(idx_max_slope);
    
    % Schnittpunkte der Tangente mit den Linien y=0 und y=1
    t_0 = -b / m;                 % Schnittpunkt mit y=0 (t-Achse)
    t_1 = (1 - b) / m;            % Schnittpunkt mit y=1
    
    % Bestimme Totzeit L und Anstiegszeit T
    L = t_0;
    T = t_1 - t_0;
    
    % Stationäre Verstärkung
    Ks = dcgain(G);
    
    % Parameter tau = T/L, hilfreich für Cohen-Coon-Formeln
    tau = T/L;
    
    details = sprintf('Ks = %.4f\nL = %.4f s\nT = %.4f s\ntau = %.4f', Ks, L, T, tau);
    
    % Parameter nach Cohen-Coon-Tabelle
    switch structure
        case 'P'
            Kp = (1/Ks) * (1 + (1/3)*tau);
            K = tf(Kp, 1);
        case 'PI'
            Kp = (1/Ks) * (0.9 + (1/12)*tau);
            Ti = L * ((30 + 3*tau)/(9 + 20*tau));
            K = tf([Kp, Kp/Ti], [1, 0]);
        case 'PD'
            Kp = (1/Ks) * (1.25 * (1 + (1/6)*tau));
            Td = L * ((6 - 2*tau)/(22 + 3*tau));
            K = tf([Kp*Td, Kp], [epsilon*Td, 1]);
        case 'PID'
            Kp = (1/Ks) * (1.35 + (1/4)*tau);
            Ti = L * ((32 + 6*tau)/(13 + 8*tau));
            Td = L * (4/(11 + 2*tau));
            K = tf([Kp*Td, Kp, Kp/Ti], [epsilon*Td, 1, 0]);
        otherwise
            error('Unsupported controller structure for Cohen-Coon method');
    end
end

% 6. Loop-Shaping
function [K, details] = designLoopShaping(G, structure, phaseMargin, bandwidth, epsilon)
    % Berechne Bode-Diagramm der Strecke
    [mag, phase, wout] = bode(G, {1e-3, 1e3});
    mag = squeeze(mag);
    phase = squeeze(phase);
    
    % Finde Index für die gewünschte Bandbreite
    [~, idx] = min(abs(wout - bandwidth));
    phase_at_bw = phase(idx);
    mag_at_bw = mag(idx);
    
    % Berechne erforderliche Phasenanhebung
    required_phase = phaseMargin - (180 + phase_at_bw);
    
    details = sprintf('Gewünschte Bandbreite: %.2f rad/s\nPhase bei Bandbreite: %.2f°\nErforderliche Phasenanhebung: %.2f°', ...
        bandwidth, phase_at_bw, required_phase);
    
    % Regler basierend auf Struktur
    switch structure
        case 'P'
            % Einfacher P-Regler
            Kp = 1/mag_at_bw;  % Verstärkung für 0dB bei Bandbreite
            K = tf(Kp, 1);
            
        case 'PI'
            % PI-Regler: Kp(1 + 1/(Ti*s))
            Kp = 1/mag_at_bw * 0.8;  % Leicht reduziert wegen PI-Phase
            Ti = 10/bandwidth;  % Integrationszeit deutlich unter Bandbreite
            K = tf([Kp, Kp/Ti], [1, 0]);
            
        case 'PD'
            % PD-Regler: Kp(1 + Td*s)
            if required_phase > 0
                % Berechne Td für gewünschte Phasenanhebung
                alpha = (1 - sin(required_phase*pi/180)) / (1 + sin(required_phase*pi/180));
                Td = 1/(bandwidth * sqrt(alpha));
                
                % Berechne Kp für 0dB bei Bandbreite
                [mag_pd, ~] = bode(tf([Td, 1], [epsilon*Td, 1]), bandwidth);
                Kp = 1 / (mag_at_bw * mag_pd);
                
                K = tf(Kp*[Td, 1], [epsilon*Td, 1]);
            else
                % Keine Phasenanhebung notwendig
                Kp = 1/mag_at_bw;
                Td = 0.1/bandwidth;
                K = tf(Kp*[Td, 1], [epsilon*Td, 1]);
            end
            
        case 'PID'
            % PID-Regler: Kp(1 + 1/(Ti*s) + Td*s)
            if required_phase > 0
                % Berechne Td für gewünschte Phasenanhebung
                alpha = (1 - sin(required_phase*pi/180)) / (1 + sin(required_phase*pi/180));
                Td = 1/(bandwidth * sqrt(alpha));
                
                % Berechne Kp für 0dB bei Bandbreite
                [mag_pd, ~] = bode(tf([Td, 1], [epsilon*Td, 1]), bandwidth);
                Kp = 1 / (mag_at_bw * mag_pd) * 0.8;  % Leicht reduziert
                
                Ti = 5/bandwidth;  % Integrationszeit unter Bandbreite
                
                K = tf(Kp*[Td, 1, 1/Ti], [epsilon*Td, 1, 0]);
            else
                % Keine Phasenanhebung notwendig
                Kp = 1/mag_at_bw;
                Td = 0.1/bandwidth;
                Ti = 5/bandwidth;
                K = tf(Kp*[Td, 1, 1/Ti], [epsilon*Td, 1, 0]);
            end
            
        otherwise
            error('Unsupported controller structure for Loop-Shaping method');
    end
    
    % Überprüfe den entworfenen Regler
    L = G * K;
    try
        margins = margin(L);
        pm_achieved = margins(2);
        details = [details, sprintf('\nErreichte Phasenreserve: %.2f°', pm_achieved)];
    catch
        % Bei Fehler in margin, keine zusätzliche Info
    end
end

% 7. IMC (Internal Model Control)
function [K, details] = designIMC(G, structure, settlingTime, epsilon)
    % Ziel-Geschwindigkeitsparameter für IMC-Filter basierend auf Einschwingzeit
    % Schätzung: Für PT1-System ist Einschwingzeit ca. 4*Zeitkonstante
    lambda = settlingTime / 4;
    
    % Zerlege Strecke in minimalphasige und nicht-minimalphasige Teile
    [z, p, k] = zpkdata(G, 'v');
    
    % Filter für nicht-minimalphasige Nullstellen
    Gp = zpk([], p, k);
    
    for i = 1:length(z)
        if real(z(i)) <= 0
            % Minimalphasige Nullstelle - zu Gp hinzufügen
            Gp = zpk([Gp.z; z(i)], Gp.p, Gp.k);
        else
            % Nicht-minimalphasige Nullstelle - durch gespiegelte Nullstelle ersetzen
            % (wird im IMC-Filter berücksichtigt)
        end
    end
    
    details = sprintf('Lambda = %.4f (basierend auf Einschwingzeit)\nSystemordnung = %d', lambda, length(p));
    
    % IMC-Filter je nach Reglertyp anpassen
    switch structure
        case 'P'
            % Einfacher IMC-Filter erster Ordnung
            F = tf(1, [lambda, 1]);
            K = minreal(inv(Gp) * F / (1 - F * G));
            
        case 'PI'
            % IMC-Filter für PI-Charakteristik (mit einem I-Term)
            if any(abs(p) < 1e-6)
                % Strecke enthält bereits einen Integrator
                F = tf(1, [lambda, 1]);
            else
                % Filter mit I-Term
                F = tf(1, [lambda, 1]) * tf(1, [1, 0]);
            end
            K = minreal(inv(Gp) * F / (1 - F * G));
            
        case 'PD'
            % IMC-Filter für PD-Charakteristik
            F = tf([lambda/4, 1], [lambda, 1]);
            K_raw = minreal(inv(Gp) * F / (1 - F * G));
            
            % Filter-Implementierung für D-Anteil
            [num, den] = tfdata(K_raw, 'v');
            if length(num) > length(den)
                % D-Anteil filtern
                num_filt = num;
                new_den = conv(den, [epsilon*lambda/4, 1]);
                K = tf(num_filt, new_den);
            else
                K = K_raw;
            end
            
        case 'PID'
            % IMC-Filter für PID-Charakteristik
            if any(abs(p) < 1e-6)
                % Strecke enthält bereits einen Integrator
                F = tf([lambda/4, 1], [lambda, 1]);
            else
                % Filter mit I-Term
                F = tf([lambda/4, 1], [lambda, 1]) * tf(1, [1, 0]);
            end
            K_raw = minreal(inv(Gp) * F / (1 - F * G));
            
            % Filter-Implementierung für D-Anteil
            [num, den] = tfdata(K_raw, 'v');
            if length(num) > length(den)
                % D-Anteil filtern
                num_filt = num;
                new_den = conv(den, [epsilon*lambda/4, 1]);
                K = tf(num_filt, new_den);
            else
                K = K_raw;
            end
            
        otherwise
            error('Unsupported controller structure for IMC method');
    end
    
    % Sicherstellen, dass der Regler eine angemessene Form hat
    try
        % Reduziere komplexe Übertragungsfunktion
        K = minreal(K);
        
        % Stelle sicher, dass der Regler proper ist
        [num, den] = tfdata(K, 'v');
        if length(num) > length(den)
            % Filter hinzufügen, um proper zu werden
            filter_const = lambda/10;
            for i = 1:(length(num) - length(den))
                den = conv(den, [filter_const, 1]);
                filter_const = filter_const/10;  % Für mehrere Filter unterschiedliche Zeitkonstanten
            end
            K = tf(num, den);
        end
    catch
        % Falls ein Fehler auftritt, versuche eine strukturierte Näherung
        try
            % Frequenzgang des Reglers auswerten
            w = logspace(-3, 3, 100);
            [mag, phase] = bode(K, w);
            mag = squeeze(mag);
            phase = squeeze(phase);
            
            % Approximiere mit der gewünschten Struktur
            switch structure
                case 'P'
                    Kp = mag(1);  % Tieffrequenz-Verhalten
                    K = tf(Kp, 1);
                case 'PI'
                    Kp = mag(end/2);  % Mittlere Frequenz
                    Ki = mag(1) * w(1);  % Tieffrequenz I-Wirkung
                    K = tf([Kp, Ki], [1, 0]);
                case 'PID'
                    Kp = mag(end/2);  % Mittlere Frequenz
                    Ki = mag(1) * w(1);  % Tieffrequenz I-Wirkung
                    Kd = mag(end) / w(end);  % Hochfrequenz D-Wirkung
                    K = tf([Kd, Kp, Ki], [epsilon*Kd/Kp, 1, 0]);
                otherwise
                    error('Fallback approximation not implemented for this structure');
            end
        catch
            error('Failed to create a proper controller structure from IMC design');
        end
    end
end

% 8. MIGO (M-constrained Integral Gain Optimization)
function [K, details] = designMIGO(G, structure, robustness, epsilon)
    % MIGO implementiert eine vereinfachte Version der M-constrained Integral Gain Optimization
    
    % Setze M-Beschränkung basierend auf Robustheitswahl
    switch robustness
        case 'Gering'
            M_s = 2.0;  % Geringe Robustheit, höhere Performance
        case 'Low'
            M_s = 2.0;  % Geringe Robustheit, höhere Performance
        case 'Mittel'
            M_s = 1.5;  % Mittlere Robustheit
        case 'Medium'
            M_s = 1.5;  % Mittlere Robustheit
        case 'Hoch'
            M_s = 1.2;  % Hohe Robustheit, geringere Performance
        case 'High'
            M_s = 1.2;  % Hohe Robustheit, geringere Performance
        otherwise
            M_s = 1.5;  % Standardwert
    end
    
    % Für diese vereinfachte Implementierung verwenden wir ein Raster-Verfahren
    % zur Suche nach geeigneten Parametern für PI/PID-Regler
    details = sprintf('Robustheitsgrad: %s\nM_s Beschränkung: %.2f', robustness, M_s);
    
    % Implementierung basierend auf der Struktur
    switch structure
        case 'P'
            % Für P-Regler vereinfacht: Bestimme maximale Verstärkung unter M_s-Beschränkung
            Kp = findOptimalKp(G, M_s);
            K = tf(Kp, 1);
            details = [details sprintf('\nOptimale P-Verstärkung: Kp = %.4f', Kp)];
            
        case 'PI'
            % Optimiere Kp und Ki unter M_s-Beschränkung
            [Kp, Ki] = findOptimalPI(G, M_s);
            K = tf([Kp, Ki], [1, 0]);
            Ti = Kp/Ki;
            details = [details sprintf('\nOptimale PI-Parameter:\nKp = %.4f\nKi = %.4f\nTi = %.4f', Kp, Ki, Ti)];
            
        case 'PD'
            % Für PD-Regler: Optimiere Kp und Kd
            [Kp, Kd] = findOptimalPD(G, M_s, epsilon);
            K = tf([Kd, Kp], [epsilon*Kd, 1]);
            Td = Kd/Kp;
            details = [details sprintf('\nOptimale PD-Parameter:\nKp = %.4f\nKd = %.4f\nTd = %.4f', Kp, Kd, Td)];
            
        case 'PID'
            % Optimiere Kp, Ki und Kd unter M_s-Beschränkung
            % Vereinfachte Implementierung für das Beispiel
            [Kp, Ki, Kd] = findOptimalPID(G, M_s, epsilon);
            K = tf([Kd, Kp, Ki], [epsilon*Kd, 1, 0]);
            Ti = Kp/Ki;
            Td = Kd/Kp;
            details = [details sprintf('\nOptimale PID-Parameter:\nKp = %.4f\nKi = %.4f\nKd = %.4f\nTi = %.4f\nTd = %.4f', ...
                Kp, Ki, Kd, Ti, Td)];
            
        otherwise
            error('MIGO ist nur für P, PI, PD und PID-Regler implementiert');
    end
    
    % Innere Funktionen für die Optimierung
    function Kp = findOptimalKp(G, M_s)
        % Finde die maximale P-Verstärkung, die die M_s-Beschränkung einhält
        Kp_min = 0.01;
        Kp_max = 100;
        Kp_best = Kp_min;
        
        for Kp = linspace(Kp_min, Kp_max, 50)
            K_test = tf(Kp, 1);
            L = G * K_test;
            S = feedback(1, L);
            
            try
                % Berechne die maximale Verstärkung von S
                [peakgain, ~, ~] = getPeakGain(S);
                
                if peakgain <= M_s && Kp > Kp_best
                    Kp_best = Kp;
                end
            catch
                % Wenn Peak-Gain-Berechnung fehlschlägt, überspringe diesen Wert
                continue;
            end
        end
        
        Kp = Kp_best;
    end
    
    function [Kp, Ki] = findOptimalPI(G, M_s)
        % Finde optimale PI-Parameter unter M_s-Beschränkung
        % Ziel ist die Maximierung des integralen Anteils
        
        % Grobes Raster für die Suche
        Kp_values = linspace(0.01, 10, 20);
        Ki_values = linspace(0.001, 5, 20);
        
        best_Kp = 0.01;
        best_Ki = 0.001;
        best_score = -Inf;
        
        for i = 1:length(Kp_values)
            for j = 1:length(Ki_values)
                Kp = Kp_values(i);
                Ki = Ki_values(j);
                
                K_test = tf([Kp, Ki], [1, 0]);
                L = G * K_test;
                S = feedback(1, L);
                
                try
                    % Berechne die maximale Verstärkung von S
                    [peakgain, ~, ~] = getPeakGain(S);
                    
                    if peakgain <= M_s
                        % Score basierend auf Ki und Phasenreserve
                        [~, Pm] = margin(L);
                        score = Ki * Pm/100;  % Gewichtung von Ki und Phasenreserve
                        
                        if score > best_score
                            best_score = score;
                            best_Kp = Kp;
                            best_Ki = Ki;
                        end
                    end
                catch
                    % Wenn Peak-Gain-Berechnung fehlschlägt, überspringe diese Kombination
                    continue;
                end
            end
        end
        
        Kp = best_Kp;
        Ki = best_Ki;
    end
    
    function [Kp, Kd] = findOptimalPD(G, M_s, epsilon)
        % Finde optimale PD-Parameter unter M_s-Beschränkung
        % Ziel ist die Maximierung der Phasenreserve
        
        % Grobes Raster für die Suche
        Kp_values = linspace(0.01, 10, 15);
        Td_values = linspace(0.01, 5, 15);
        
        best_Kp = 0.01;
        best_Kd = 0.01;
        best_Pm = 0;
        
        for i = 1:length(Kp_values)
            for j = 1:length(Td_values)
                Kp = Kp_values(i);
                Td = Td_values(j);
                Kd = Kp * Td;
                
                K_test = tf([Kd, Kp], [epsilon*Kd, 1]);
                L = G * K_test;
                S = feedback(1, L);
                
                try
                    % Berechne die maximale Verstärkung von S und Phasenreserve
                    [peakgain, ~, ~] = getPeakGain(S);
                    [~, Pm] = margin(L);
                    
                    if peakgain <= M_s && Pm > best_Pm
                        best_Pm = Pm;
                        best_Kp = Kp;
                        best_Kd = Kd;
                    end
                catch
                    % Wenn Berechnung fehlschlägt, überspringe diese Kombination
                    continue;
                end
            end
        end
        
        Kp = best_Kp;
        Kd = best_Kd;
    end
    
    function [Kp, Ki, Kd] = findOptimalPID(G, M_s, epsilon)
        % Finde optimale PID-Parameter unter M_s-Beschränkung
        % Vereinfachte Suche in einem 3D-Parameter-Raum
        
        % Grobes Raster für die Suche
        Kp_values = linspace(0.01, 10, 8);
        Ti_values = linspace(0.1, 10, 8);
        Td_values = linspace(0.01, 2, 8);
        
        best_Kp = 0.01;
        best_Ki = 0.001;
        best_Kd = 0.01;
        best_score = -Inf;
        
        for i = 1:length(Kp_values)
            for j = 1:length(Ti_values)
                for k = 1:length(Td_values)
                    Kp = Kp_values(i);
                    Ti = Ti_values(j);
                    Td = Td_values(k);
                    
                    Ki = Kp / Ti;
                    Kd = Kp * Td;
                    
                    K_test = tf([Kd, Kp, Ki], [epsilon*Kd, 1, 0]);
                    L = G * K_test;
                    S = feedback(1, L);
                    
                    try
                        % Berechne die maximale Verstärkung von S und Phasenreserve
                        [peakgain, ~, ~] = getPeakGain(S);
                        [~, Pm] = margin(L);
                        
                        if peakgain <= M_s
                            % Score basierend auf Ki, Kd und Phasenreserve
                            score = Ki * sqrt(Kd) * Pm/100;
                            
                            if score > best_score
                                best_score = score;
                                best_Kp = Kp;
                                best_Ki = Ki;
                                best_Kd = Kd;
                            end
                        end
                    catch
                        % Wenn Berechnung fehlschlägt, überspringe diese Kombination
                        continue;
                    end
                end
            end
        end
        
        Kp = best_Kp;
        Ki = best_Ki;
        Kd = best_Kd;
    end
end

% H-infinity design method
function [K, details] = designHInfinity(G, structure, robustness, epsilon)
    % H-infinity controller design focusing on robust stability and performance
    % The method implements a simplified H-infinity synthesis for SISO systems
    
    % Map robustness setting to gamma value (smaller = more robust but conservative)
    switch robustness
        case 'Low'
            gamma = 3.0;  % Less conservative
        case 'Medium'
            gamma = 2.0;  % Balanced approach
        case 'High'
            gamma = 1.2;  % More robust
        otherwise
            gamma = 2.0;  % Default
    end
    
    details = sprintf('H-infinity Design\n----------------\nRobustness level: %s\nGamma value: %.2f\n', robustness, gamma);
    
    try
        % Convert to state-space for H-infinity synthesis
        [A, B, C, D] = ssdata(G);
        
        % For MIMO systems, we focus on SISO for this application
        nx = size(A, 1);  % Number of states
        
        % Define weighting functions based on desired properties
        if strcmpi(robustness, 'High')
            % Higher robustness: emphasize uncertainty rejection
            Ws = tf([1, 0.5], [0.01, 1]);  % Sensitivity weight (error)
            Wks = tf([10, 0], [0.001, 1]);  % KS weight (control effort)
            Wt = tf([0.01, 1], [1, 10]);  % Complementary sensitivity weight
        elseif strcmpi(robustness, 'Low')
            % Lower robustness: emphasize performance
            Ws = tf([1, 0.1], [0.1, 1]);  % Sensitivity weight
            Wks = tf([1, 0], [0.01, 1]);  % KS weight
            Wt = tf([0.1, 1], [1, 5]);  % Complementary sensitivity weight
        else
            % Medium robustness: balanced approach
            Ws = tf([1, 0.2], [0.05, 1]);  % Sensitivity weight
            Wks = tf([5, 0], [0.005, 1]);  % KS weight
            Wt = tf([0.05, 1], [1, 8]);  % Complementary sensitivity weight
        end
        
        details = [details, sprintf('\nSensitivity weight Ws: %s\nControl effort weight Wks: %s\nRobustness weight Wt: %s\n', ...
                 char(Ws), char(Wks), char(Wt))];
        
        % For simple cases, we'll use the mixed-sensitivity approach
        % This is a simplified version of what hinfsyn would do
        
        % 1. Augment the plant with weights
        % Here we implement a basic mixed-sensitivity setup with:
        % - S = 1/(1+GK) (sensitivity)
        % - KS = K/(1+GK) (control sensitivity)
        % - T = GK/(1+GK) (complementary sensitivity)
        
        % Create augmented plant for mixed-sensitivity synthesis
        % State-space representation of the weights
        [Aws, Bws, Cws, Dws] = ssdata(Ws);
        [Awks, Bwks, Cwks, Dwks] = ssdata(Wks);
        [Awt, Bwt, Cwt, Dwt] = ssdata(Wt);
        
        nws = size(Aws, 1);
        nwks = size(Awks, 1);
        nwt = size(Awt, 1);
        
        % Augmented plant matrices
        Aaug = [A, zeros(nx, nws), zeros(nx, nwks), zeros(nx, nwt);
                Bws*C, Aws, zeros(nws, nwks), zeros(nws, nwt);
                zeros(nwks, nx), zeros(nwks, nws), Awks, zeros(nwks, nwt);
                Bwt*C, zeros(nwt, nws), zeros(nwt, nwks), Awt];
        
        Baug = [B; zeros(nws, 1); Bwks; zeros(nwt, 1)];
        
        Caug = [Dws*C, Cws, zeros(1, nwks), zeros(1, nwt);
                zeros(1, nx), zeros(1, nws), Cwks, zeros(1, nwt);
                Dwt*C, zeros(1, nws), zeros(1, nwks), Cwt];
        
        Daug = [Dws*D; Dwks; Dwt*D];
        
        % Simplified H-infinity synthesis using Riccati equations
        % For SISO systems, we can implement a simplified version
        
        % For brevity, we'll use a structure-based approach similar to loop-shaping
        % as a fallback, which is more practical for the auto-tuning UI
        
        % Determine controller parameters based on structure
        switch structure
            case 'P'
                % P controller: gain determined from H-infinity objective
                K_hinf = loopsyn(G, Ws);
                [num, den] = tfdata(K_hinf, 'v');
                
                % Extract proportional gain
                Kp = num(1) / den(1);  % DC gain
                K = tf(Kp, 1);
                
                details = [details, sprintf('\nP controller with gain: Kp = %.4f', Kp)];
                
            case 'PI'
                % PI controller from H-infinity objective
                K_hinf = loopsyn(G, Ws);
                
                % Extract approximate PI parameters from H-infinity controller
                [num, den] = tfdata(K_hinf, 'v');
                
                % Find dominant poles/zeros for PI approximation
                if length(den) > 1
                    % Extract PI parameters
                    Kp = num(1) / den(1);  % Approximate proportional gain
                    
                    % Try to extract integral action
                    z = roots(num);
                    p = roots(den);
                    
                    % Look for integrator or slow pole
                    integrator_idx = find(abs(p) < 0.01, 1);
                    if ~isempty(integrator_idx)
                        Ki = 0.1 * Kp;  % Reasonable integral gain
                    else
                        % Estimate based on slowest pole
                        [~, idx] = min(abs(p));
                        Ki = Kp * abs(p(idx));
                    end
                else
                    % Simple case
                    Kp = num(1) / den(1);
                    Ki = Kp * 0.1;  % Conservative integral action
                end
                
                K = tf([Kp, Ki], [1, 0]);
                
                details = [details, sprintf('\nPI controller with:\nKp = %.4f\nKi = %.4f', Kp, Ki)];
                
            case 'PD'
                % PD controller with filter
                K_hinf = loopsyn(G, Ws);
                
                % Extract approximate PD parameters
                [num, den] = tfdata(K_hinf, 'v');
                
                if length(num) > 1
                    % Extract PD gains with filter
                    Kp = num(end) / den(end);  % DC gain
                    
                    % Get derivative action
                    if length(num) >= 2
                        Kd = num(1) / den(end) * epsilon;  % Scale for filter
                    else
                        Kd = Kp * 0.1;  % Conservative default
                    end
                else
                    % Simple case
                    Kp = num(1) / den(1);
                    Kd = Kp * 0.1;  % Conservative derivative action
                end
                
                K = tf([Kd, Kp], [epsilon*Kd, 1]);
                
                details = [details, sprintf('\nPD controller with:\nKp = %.4f\nKd = %.4f\nFilter epsilon = %.4f', Kp, Kd, epsilon)];
                
            case 'PID'
                % PID controller with filter from H-infinity objective
                K_hinf = loopsyn(G, Ws);
                
                % Extract approximate PID parameters
                [num, den] = tfdata(K_hinf, 'v');
                
                % Find dominant poles/zeros for PID approximation
                if length(num) >= 2 && length(den) >= 2
                    % Extract basic gain
                    Kp = num(end-1) / den(end);  % Proportional term
                    
                    % Look for integral action
                    z = roots(num);
                    p = roots(den);
                    
                    % Look for integrator or slow pole
                    integrator_idx = find(abs(p) < 0.01, 1);
                    if ~isempty(integrator_idx)
                        Ki = 0.1 * Kp;  % Reasonable integral gain
                    else
                        % Estimate based on slowest pole/zero
                        [~, idx] = min(abs(p));
                        Ki = Kp * abs(p(idx));
                    end
                    
                    % Derivative action
                    if length(num) >= 3
                        Kd = num(1) / num(end-1) * Kp * epsilon;  % Scale for filter
                    else
                        Kd = Kp * 0.1;  % Conservative derivative action
                    end
                else
                    % Simplified case
                    Kp = 1.0;
                    Ki = 0.5;
                    Kd = 0.1;
                    
                    % Scale based on plant characteristics
                    dc_gain = dcgain(G);
                    if ~isnan(dc_gain) && dc_gain ~= 0
                        Kp = 1/abs(dc_gain);
                    end
                end
                
                % Create PID controller with filtered derivative
                K = tf([Kd, Kp, Ki], [epsilon*Kd, 1, 0]);
                
                details = [details, sprintf('\nPID controller with:\nKp = %.4f\nKi = %.4f\nKd = %.4f\nFilter epsilon = %.4f', ...
                          Kp, Ki, Kd, epsilon)];
            otherwise
                error('Unsupported controller structure for H-infinity method');
        end
        
        % Check controller stability
        K_poles = pole(K);
        if any(real(K_poles) > 0)
            details = [details, '\nWarning: Controller contains unstable poles.'];
        end
        
        % Check closed-loop stability
        cl_sys = feedback(G*K, 1);
        cl_poles = pole(cl_sys);
        if any(real(cl_poles) > 0)
            details = [details, '\nWarning: Closed-loop system is unstable.'];
        else
            % Performance metrics if stable
            try
                [Gm, Pm] = margin(G*K);
                details = [details, sprintf('\nClosed-loop performance:\nGain Margin: %.2f dB\nPhase Margin: %.2f degrees', ...
                         20*log10(Gm), Pm)];
            catch
                details = [details, '\nCould not compute stability margins.'];
            end
        end
    catch ME
        % Handle any errors during the design process
        warning('H-infinity design error: %s', ME.message);
        
        % Fallback to simple controller
        switch structure
            case 'P'
                Kp = 1.0;
                K = tf(Kp, 1);
            case 'PI'
                Kp = 1.0;
                Ki = 0.5;
                K = tf([Kp, Ki], [1, 0]);
            case 'PD'
                Kp = 1.0;
                Kd = 0.1;
                K = tf([Kd, Kp], [epsilon*Kd, 1]);
            case 'PID'
                Kp = 1.0;
                Ki = 0.5;
                Kd = 0.1;
                K = tf([Kd, Kp, Ki], [epsilon*Kd, 1, 0]);
            otherwise
                K = tf(1, 1);
        end
        
        details = sprintf('H-infinity design failed: %s\n\nUsing fallback controller instead.', ME.message);
    end
end

% LQG design method
function [K, details] = designLQG(G, structure, bandwidth, robustness, epsilon)
    % LQG (Linear-Quadratic-Gaussian) controller design
    % Combines optimal LQR state feedback with Kalman filter state estimation
    
    details = sprintf('LQG Design\n----------\nBandwidth target: %.2f rad/s\nRobustness level: %s\n', bandwidth, robustness);
    
    try
        % Convert to state-space for LQG design
        [A, B, C, D] = ssdata(G);
        
        % Check if system is controllable and observable
        if rank(ctrb(A, B)) < size(A, 1)
            warning('System is not fully controllable.');
            details = [details, 'Warning: System is not fully controllable.\n'];
        end
        
        if rank(obsv(A, C)) < size(A, 1)
            warning('System is not fully observable.');
            details = [details, 'Warning: System is not fully observable.\n'];
        end
        
        % Set weighting matrices based on robustness and bandwidth
        switch robustness
            case 'Low'
                % Lower robustness: emphasize performance
                Q = C' * C;  % Output error
                R = 0.1;     % Lower penalty on control effort
                Qn = 10;     % Process noise covariance
                Rn = 1;      % Measurement noise covariance
            case 'High'
                % Higher robustness: more conservative control
                Q = C' * C;  % Output error
                R = 10;      % Higher penalty on control effort
                Qn = 1;      % Process noise covariance
                Rn = 10;     % Measurement noise covariance
            otherwise
                % Medium robustness: balanced approach
                Q = C' * C;  % Output error
                R = 1;       % Medium penalty on control effort
                Qn = 1;      % Process noise covariance
                Rn = 1;      % Measurement noise covariance
        end
        
        % Adjust Q to target desired bandwidth
        Q = Q * (bandwidth^2);
        
        % Design LQR controller (state feedback)
        [K_lqr, S, e] = lqr(A, B, Q, R);
        
        % Design Kalman filter (state estimator)
        [Kf, P, E] = lqe(A, eye(size(A, 1)), C, Qn, Rn);
        
        % Compute LQG controller (combines LQR and Kalman filter)
        Ac = A - B*K_lqr - Kf*C;
        Bc = Kf;
        Cc = -K_lqr;
        Dc = 0;
        
        K_lqg = ss(Ac, Bc, Cc, Dc);
        
        details = [details, sprintf('\nLQR gain: [%s]\nKalman gain: [%s]\n', ...
                 mat2str(K_lqr, 4), mat2str(Kf', 4))];
        
        % For our application, approximate this as the requested controller structure
        switch structure
            case 'P'
                % P controller approximation
                [num, den] = tfdata(K_lqg, 'v');
                Kp = dcgain(K_lqg);  % Use DC gain as proportional term
                K = tf(Kp, 1);
                
                details = [details, sprintf('\nApproximated as P controller with:\nKp = %.4f', Kp)];
                
            case 'PI'
                % Convert to transfer function
                K_tf = tf(K_lqg);
                
                % PI approximation
                [num, den] = tfdata(K_tf, 'v');
                
                % Extract dominant characteristics
                dc_gain = dcgain(K_tf);
                
                % Look for integrator in controller
                p = pole(K_tf);
                has_integrator = any(abs(p) < 1e-6);
                
                if has_integrator
                    % If there's an integrator, we need a PI controller
                    Kp = abs(num(1) / den(1)); % Approximate gain
                    Ki = Kp * 0.5; % Approximate integral gain based on dominant pole
                else
                    % If no integrator, add one with conservative settings
                    Kp = abs(dc_gain);
                    Ki = Kp * 0.1; % Small integral action
                end
                
                K = tf([Kp, Ki], [1, 0]);
                
                details = [details, sprintf('\nApproximated as PI controller with:\nKp = %.4f\nKi = %.4f', Kp, Ki)];
                
            case 'PD'
                % Convert to transfer function
                K_tf = tf(K_lqg);
                
                % PD approximation
                [num, den] = tfdata(K_tf, 'v');
                
                % Extract proportional and derivative terms
                if length(num) >= 2
                    % Has derivative action
                    Kp = num(end) / den(end); % DC gain
                    Kd = num(1) / den(end); % High frequency gain
                else
                    % No derivative action evident
                    Kp = num(1) / den(1);
                    Kd = Kp * 0.1; % Add small derivative action
                end
                
                % Apply filter
                K = tf([Kd, Kp], [epsilon*Kd, 1]);
                
                details = [details, sprintf('\nApproximated as PD controller with:\nKp = %.4f\nKd = %.4f\nFilter epsilon = %.4f', ...
                         Kp, Kd, epsilon)];
                
            case 'PID'
                % Convert to transfer function
                K_tf = tf(K_lqg);
                
                % PID approximation
                [num, den] = tfdata(K_tf, 'v');
                
                % Check for integrator
                p = pole(K_tf);
                has_integrator = any(abs(p) < 1e-6);
                
                if has_integrator
                    % Extract PID parameters from LQG controller
                    dc_gain = 1; % Placeholder
                    
                    % Find a suitable PID approximation based on frequency response
                    w = logspace(-3, 3, 100);
                    [mag, phase] = bode(K_tf, w);
                    mag = squeeze(mag);
                    
                    % Approximate as PID based on frequency response
                    Kp = mag(floor(end/2)); % Mid-frequency gain
                    Ki = mag(1) * w(1); % Low-frequency integral effect
                    Kd = mag(end) / w(end); % High-frequency derivative effect
                else
                    % Create a PID with reasonable values
                    dc_gain = dcgain(K_tf);
                    Kp = abs(dc_gain);
                    Ki = Kp * 0.1;
                    Kd = Kp * 0.1;
                end
                
                % Create filtered PID
                K = tf([Kd, Kp, Ki], [epsilon*Kd, 1, 0]);
                
                details = [details, sprintf('\nApproximated as PID controller with:\nKp = %.4f\nKi = %.4f\nKd = %.4f\nFilter epsilon = %.4f', ...
                         Kp, Ki, Kd, epsilon)];
                
            otherwise
                error('Unsupported controller structure for LQG method');
        end
        
        % Check controller stability
        K_poles = pole(K);
        if any(real(K_poles) > 0)
            details = [details, '\nWarning: Controller contains unstable poles.'];
        end
        
        % Check closed-loop stability
        cl_sys = feedback(G*K, 1);
        cl_poles = pole(cl_sys);
        if any(real(cl_poles) > 0)
            details = [details, '\nWarning: Closed-loop system is unstable.'];
        else
            % Performance metrics if stable
            try
                [Gm, Pm] = margin(G*K);
                details = [details, sprintf('\nClosed-loop performance:\nGain Margin: %.2f dB\nPhase Margin: %.2f degrees', ...
                         20*log10(Gm), Pm)];
            catch
                details = [details, '\nCould not compute stability margins.'];
            end
        end
    catch ME
        % Handle any errors during the design process
        warning('LQG design error: %s', ME.message);
        
        % Fallback to simple controller
        switch structure
            case 'P'
                Kp = 1.0;
                K = tf(Kp, 1);
            case 'PI'
                Kp = 1.0;
                Ki = 0.5;
                K = tf([Kp, Ki], [1, 0]);
            case 'PD'
                Kp = 1.0;
                Kd = 0.1;
                K = tf([Kd, Kp], [epsilon*Kd, 1]);
            case 'PID'
                Kp = 1.0;
                Ki = 0.5;
                Kd = 0.1;
                K = tf([Kd, Kp, Ki], [epsilon*Kd, 1, 0]);
            otherwise
                K = tf(1, 1);
        end
        
        details = sprintf('LQG design failed: %s\n\nUsing fallback controller instead.', ME.message);
    end
end

% Hilfsfunktion zur Berechnung der maximalen Verstärkung
function [peakgain, wpeak, w] = getPeakGain(sys)
    % Berechnet den maximalen Amplitudengang eines Systems
    w = logspace(-3, 3, 1000);
    [mag, ~] = bode(sys, w);
    mag = squeeze(mag);
    [peakgain, idx] = max(mag);
    wpeak = w(idx);
end

% EVALUATECONTROLLER Evaluates a controller based on specified criteria
function score = evaluateController(K, G, goal, desired_pm, desired_os, desired_ts, desired_bw)
    % Initialize score
    score = 0;
    
    % Ensure goal is in correct format
    if strcmpi(goal, 'tracking')
        goal = 'Tracking';
    elseif strcmpi(goal, 'disturbance rejection')
        goal = 'Disturbance Rejection';
    elseif strcmpi(goal, 'robustness')
        goal = 'Robustness';
    end
    
    % Compute closed-loop transfer functions
    try
        T = feedback(G*K, 1);      % Complementary sensitivity function
        S = feedback(1, G*K);      % Sensitivity function
        L = G*K;                   % Open-loop transfer function
        
        % Test if closed-loop system is stable
        if any(real(pole(T)) > 0)
            disp('Controller results in unstable closed-loop system.');
            score = -100;  % Unstable system gets a heavily negative score
            return;
        end
    catch ME
        disp(['Error computing closed-loop transfer functions: ', ME.message]);
        score = -50;  % Error indicates likely problems with the controller
        return;
    end
    
    % 1. Stability metrics (30 points maximum)
    try
        % Compute gain and phase margins
        [gm, pm, wgm, wpm] = margin(L);
        
        % Convert gain margin to dB
        gm_dB = 20*log10(gm);
        
        % Phase margin scoring (0-15 points)
        pm_score = 15 * (1 - min(1, abs(pm - desired_pm) / 45));
        
        % Gain margin scoring (0-15 points)
        gm_score = 15 * (1 - min(1, abs(gm_dB - 10) / 10));
        
        % Display stability metrics
        disp(['Phase Margin: ', num2str(pm), '° (desired: ', num2str(desired_pm), '°)']);
        disp(['Gain Margin: ', num2str(gm_dB), ' dB']);
        
        % Add to total score
        score = score + pm_score + gm_score;
    catch ME
        disp(['Warning: Could not compute stability margins: ', ME.message]);
        % Apply penalty for not being able to compute stability margins
        score = score - 15;
    end
    
    % 2. Time-domain performance (40 points maximum)
    try
        % Compute step response
        [y, t] = step(T);
        
        % Extract step response metrics
        info = stepinfo(y, t);
        
        % Calculate key metrics
        actual_overshoot = info.Overshoot;
        actual_settling = info.SettlingTime;
        actual_rise = info.RiseTime;
        
        % Display time domain metrics
        disp(['Overshoot: ', num2str(actual_overshoot), '% (desired: ', num2str(desired_os), '%)']);
        disp(['Settling Time: ', num2str(actual_settling), ' s (desired: ', num2str(desired_ts), ' s)']);
        disp(['Rise Time: ', num2str(actual_rise), ' s']);
        
        % Score overshoot (0-15 points)
        % Better score for being close to desired overshoot
        os_score = 15 * (1 - min(1, abs(actual_overshoot - desired_os) / 50));
        
        % Score settling time (0-15 points)
        % Better score for faster settling
        ts_score = 15 * (1 - min(1, abs(actual_settling - desired_ts) / (2*desired_ts)));
        
        % Score rise time (0-10 points)
        % Better score for faster rise time relative to settling time
        rt_score = 10 * (1 - min(1, actual_rise / desired_ts));
        
        % Add to total score
        score = score + os_score + ts_score + rt_score;
    catch ME
        disp(['Warning: Could not compute time-domain performance metrics: ', ME.message]);
        % Apply penalty for not being able to compute time domain metrics
        score = score - 20;
    end
    
    % 3. Frequency-domain performance (20 points maximum)
    try
        % Calculate bandwidth
        bw = bandwidth(T);
        
        % Display bandwidth
        disp(['Bandwidth: ', num2str(bw), ' rad/s (desired: ', num2str(desired_bw), ' rad/s)']);
        
        % Score bandwidth (0-20 points)
        % Better score for being close to desired bandwidth
        bw_score = 20 * (1 - min(1, abs(bw - desired_bw) / desired_bw));
        
        % Add to total score
        score = score + bw_score;
    catch ME
        disp(['Warning: Could not compute bandwidth: ', ME.message]);
        % Apply penalty for not being able to compute bandwidth
        score = score - 10;
    end
    
    % 4. Goal-specific metrics (10 points maximum)
    try
        switch goal
            case 'Tracking'
                % For tracking, evaluate the reference tracking performance
                % Calculate integral of error for a step reference
                t_sim = linspace(0, 3*desired_ts, 1000);
                r = ones(size(t_sim));
                [y_r, ~] = lsim(T, r, t_sim);
                err_integral = trapz(t_sim, abs(r - y_r));
                
                % Score tracking (0-10 points)
                tracking_score = 10 * exp(-err_integral / 5);
                
                % Add to total score
                score = score + tracking_score;
                
                disp(['Tracking Performance Score: ', num2str(tracking_score), '/10']);
                
            case 'Disturbance Rejection'
                % For disturbance rejection, evaluate the disturbance response
                % Calculate the integral of output for a step disturbance
                t_sim = linspace(0, 3*desired_ts, 1000);
                d = ones(size(t_sim));
                [y_d, ~] = lsim(S, d, t_sim);  % S is sensitivity function
                dist_integral = trapz(t_sim, abs(y_d));
                
                % Score disturbance rejection (0-10 points)
                dist_score = 10 * exp(-dist_integral / 5);
                
                % Add to total score
                score = score + dist_score;
                
                disp(['Disturbance Rejection Score: ', num2str(dist_score), '/10']);
                
            case 'Robustness'
                % For robustness, evaluate the sensitivity peak
                try
                    [Ms, w_ms] = getPeakGain(S);
                    
                    % Good sensitivity peak should be below 2 (6 dB)
                    robust_score = 10 * exp(-(Ms-1) / 1.5);
                    
                    % Add to total score
                    score = score + robust_score;
                    
                    disp(['Sensitivity Peak (Ms): ', num2str(Ms)]);
                    disp(['Robustness Score: ', num2str(robust_score), '/10']);
                catch
                    disp('Could not compute sensitivity peak.');
                    score = score - 5;
                end
        end
    catch ME
        disp(['Warning: Could not compute goal-specific metrics: ', ME.message]);
        % Apply penalty for not being able to compute goal-specific metrics
        score = score - 5;
    end
    
    % Normalize score to be within 0-100 range
    score = max(0, min(100, score));
    
    % Display final score
    disp(['Final Controller Score: ', num2str(score), '/100']);
end