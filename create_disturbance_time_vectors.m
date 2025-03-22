function create_disturbance_time_vectors(app)
    % Zeitvektor für d_1 erzeugen (wenn definiert)
    if ~isempty(app.d1_num_sections) && ~isempty(app.d1_start_time) && ~isempty(app.d1_end_time)
        app.d1_time_vector = create_linear_time_vector_disturbance(app.d1_num_sections, 'd_1');
    end
    
    % Zeitvektor für d_2 erzeugen (wenn definiert)
    if ~isempty(app.d2_num_sections) && ~isempty(app.d2_start_time) && ~isempty(app.d2_end_time)
        app.d2_time_vector = create_linear_time_vector_disturbance(app.d2_num_sections, 'd_2');
    end
end